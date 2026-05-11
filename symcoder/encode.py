from __future__ import annotations
import importlib.util
import numpy as np
from itertools import groupby
from pathlib import Path
from symatom.atoms import ArgumentSymmetry
from symatom.rep import canonical_pair_flavours
from symatom import repL
from .pairs import eval_pair_orbit, eval_pair_orbit_positive, eval_single_orbit, eval_single_orbit_compressed, _is_self_pair

# Load the Echinocoder zip embedder from the repo root.  It is not a proper
# package, so we locate it relative to this file (repo_root/symcoder/encode.py
# → repo_root/Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py).
_EMBEDDER_PATH = (
    Path(__file__).resolve().parent.parent
    / "Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py"
)
_spec = importlib.util.spec_from_file_location("_echinocoder", _EMBEDDER_PATH)
_echinocoder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_echinocoder)
_zip_embed = _echinocoder.embed


def _force_positive_side_key(op, flavour):
    """
    Canonical key for one side of a pair, collapsing ANTISYMMETRIC/SYMMETRIC
    distinction into the same key (forcePositive: negating all atoms in an
    ANTISYMMETRIC row gives an equivalent z-vector, so both the signed and
    unsigned variants of that row are treated as the same 'side').
    """
    return (op.name, tuple(flavour.counts))


def _encoding_canonical_key(pf):
    """
    Canonical key for a PairFlavour under encoding-level equivalence.

    Two pairs of rows are encoding-equivalent if they produce z-vectors whose
    polynomial embeddings carry the same information, up to:
      - u↔v swap        (z_k → i·conj(z_k)):   handled by sorting the sides
      - negating row u  (z_k → −Re + i·Im):     handled by forcePositive
      - negating row v  (z_k → conj(z_k)):       handled by forcePositive

    The key is (sorted pair of side-keys, overlap), which is invariant under
    all three operations.
    """
    u_key = _force_positive_side_key(pf.op_u, pf.flavour_u)
    v_key = _force_positive_side_key(pf.op_v, pf.flavour_v)
    return (tuple(sorted([u_key, v_key])), tuple(pf.overlap))


def _embed_compressed(z_pos: np.ndarray, pf) -> np.ndarray:
    """
    Embed the sign=(+1,+1) orbit evaluations z_pos into a permutation-invariant
    vector of exactly pf.count(...) complex numbers, exploiting the polynomial
    symmetry structure that ANTISYMMETRIC operations impose.

    Background
    ----------
    The characteristic polynomial of the full orbit multiset {z_k} has forced
    structure whenever either operation is ANTISYMMETRIC, because the orbit
    includes both a z_k and its "partner" obtained by negating that operation's
    sign.  In z-space:

      SYM×ANTISYM  : partner of z = t+ib  is  z̄ = t-ib  (conjugate)
      ANTISYM×SYM  : partner of z = t+ib  is −z̄ = -t+ib (anti-conjugate)
      ANTISYM×ANTISYM: both partners present → {z, z̄, -z, -z̄}

    These structures force polynomial coefficients to be real, or purely
    imaginary, or zero, allowing a factor 2 (one ANTISYM) or 4 (both ANTISYM)
    reduction in the number of independent scalars that must be stored.

    All three non-trivial cases reduce to the same output size as SYM×SYM:
    exactly n = pf.count(group_sizes) complex numbers, where n is the number
    of distinct label combinations before signing.

    Why the sign=(+1,+1) subset
    ---------------------------
    For each base label combination, the Atom constructor sorts labels and
    absorbs permutation parity into the stored sign.  For ANTISYMMETRIC op_u
    this yields two u-atoms with stored signs +s and -s, and similarly for
    op_v.  Exactly one of the four (±u, ±v) pairs has stored (u.sign=+1,
    v.sign=+1) — so this subset always has exactly n elements regardless of
    label ordering (even pathological orderings like "toast" < "apple" <
    "zebra" where most labels sort to sign=-1 internally).

    The z_k values from this subset are z_k = eval(u_canonical, E) + i·eval(v_canonical, E),
    where "canonical" means the Atom constructor's sorted, sign-absorbed form.

    Permutation invariance of the invariant multisets
    -------------------------------------------------
    Under a label permutation that flips the sign of op_v (ANTISYM) atoms:
      z_k → z̄_k.   {z_k, z̄_k} is invariant. ✓
    Under a permutation that flips op_u sign:
      z_k → -z̄_k.  {z_k, -z̄_k} is invariant. ✓
    Under both flips (ANTISYM×ANTISYM):
      z_k → ±z_k or ±z̄_k, so z_k² → z_k² or z̄_k².
      {z_k², z̄_k²} is invariant. ✓
    Under permutations that permute base label combinations (not sign-flipping):
      z_k → z_{π(k)}: the multiset changes, but it is a multiset so order
      does not matter — the polynomial is unchanged. ✓

    Embedder coefficient convention
    --------------------------------
    The Echinocoder zip embedder computes
        polyfromroots(-data)[-2::-1]
    i.e. the polynomial with roots at -data, coefficients in DESCENDING degree
    order with the monic leading term dropped.  For n input values the output
    has n complex entries:
        coeffs[0] = coefficient of degree n-1
        coeffs[1] = coefficient of degree n-2
        ...
        coeffs[n-1] = constant term (degree 0)

    For a degree-2n polynomial (invariant multiset of 2n roots) the output
    has 2n entries; the degree-parity of entry i is (2n-1-i) mod 2:
        i even → odd degree
        i odd  → even degree

    Algorithm per symmetry class
    ----------------------------
    SYM×SYM:
        Embed {z_k} directly (n roots → n complex coefficients). No compression.

    SYM×ANTISYM (op_v ANTISYM):
        Invariant multiset: {z_k, conj(z_k)}   (conjugate-closed, 2n roots)
        Polynomial roots: {-z_k, -conj(z_k)}   (still conjugate-closed)
        → polynomial has real coefficients → coeffs[i].imag ≈ 0 for all i
        Independent reals: coeffs[0].real, coeffs[1].real, ..., coeffs[2n-1].real
        Pack as n complex: coeffs.real[0::2] + 1j * coeffs.real[1::2]

    ANTISYM×SYM (op_u ANTISYM):
        Invariant multiset: {z_k, -conj(z_k)}  (anti-conjugate-closed, 2n roots)
        Polynomial roots: {-z_k, conj(z_k)}
        Each factor: (x + z_k)(x - conj(z_k)) = (x+ib)² - t²
                   = x² + 2ib·x - (t²+b²)
        → even-degree coefficients are real, odd-degree are purely imaginary
        In descending order: even indices (0,2,4,...) hold ODD-degree (pure imag)
                             odd  indices (1,3,5,...) hold EVEN-degree (real)
        Independent reals: coeffs.imag[0::2] and coeffs.real[1::2]
        Pack as n complex: coeffs.imag[0::2] + 1j * coeffs.real[1::2]

    ANTISYM×ANTISYM:
        Invariant multiset: {z_k², conj(z_k²)}  (conjugate-closed, 2n roots)
        Same structure as SYM×ANTISYM but in the squared variable w = z².
        Polynomial has real coefficients.
        Pack as n complex: w_coeffs.real[0::2] + 1j * w_coeffs.real[1::2]

    Validation
    ----------
    Cross-check tests in test_encode.py verify that all four cases produce
    permutation-invariant output, that the output length equals pf.count(),
    and that the full-orbit embedding (eval_pair_orbit + _zip_embed, the brute
    force reference) contains the same information as the compressed form.
    """
    antisym_u = pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
    antisym_v = pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

    if not antisym_u and not antisym_v:
        # SYM×SYM: no forced structure; standard embedding.
        coeffs, _, _ = _zip_embed(z_pos)
        return coeffs

    elif not antisym_u and antisym_v:
        # SYM×ANTISYM: conjugate-closed invariant multiset → real polynomial.
        z_full = np.concatenate([z_pos, np.conj(z_pos)])
        coeffs, _, _ = _zip_embed(z_full)
        r = coeffs.real   # all imaginary parts ≈ 0
        return r[0::2] + 1j * r[1::2]

    elif antisym_u and not antisym_v:
        # ANTISYM×SYM: anti-conjugate-closed → even-real/odd-imag polynomial.
        z_full = np.concatenate([z_pos, -np.conj(z_pos)])
        coeffs, _, _ = _zip_embed(z_full)
        # even indices → odd-degree (purely imaginary); odd indices → even-degree (real)
        return coeffs.imag[0::2] + 1j * coeffs.real[1::2]

    else:
        # ANTISYM×ANTISYM: squared conjugate-closed → real polynomial in w=z².
        w = z_pos ** 2
        w_full = np.concatenate([w, np.conj(w)])
        coeffs, _, _ = _zip_embed(w_full)
        r = coeffs.real   # all imaginary parts ≈ 0
        return r[0::2] + 1j * r[1::2]


def _complex_to_reals(c: np.ndarray) -> np.ndarray:
    """Unpack n complex values to 2n reals: [re0, im0, re1, im1, ...]."""
    out = np.empty(2 * len(c), dtype=float)
    out[0::2] = c.real
    out[1::2] = c.imag
    return out


def _overlap_block_key(pf):
    """Key that identifies an OVERLAP BLOCK: fixed (op_u, flavour_u, op_v, flavour_v)."""
    u = (pf.op_u.name, pf.op_u.rank, pf.op_u.parity,
         pf.op_u.argument_symmetry.value, pf.flavour_u.counts)
    v = (pf.op_v.name, pf.op_v.rank, pf.op_v.parity,
         pf.op_v.argument_symmetry.value, pf.flavour_v.counts)
    return (u, v)


def _sort_encode(values: np.ndarray) -> np.ndarray:
    """
    Sort real eval values into a permutation-invariant float64 array.

    This is the Phase 1 encoding for a single-operator orbit.  Sorting is the
    permutation-invariant operation; the output length equals the number of
    atoms in the orbit (fo.count()), one real per atom.

    Note: for ANTISYMMETRIC operations the orbit contains {a_k, -a_k} pairs,
    so the sorted list is antisymmetric about 0 and only half the values carry
    independent information.  A future optimisation (announced) will exploit
    this to halve Phase 1 storage for ANTISYMMETRIC operators.
    """
    if len(values) == 0:
        return np.array([], dtype=float)
    return np.sort(np.asarray(values, dtype=float))


def encode_brute(plan, event: dict) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector (5a only, no 5b).

    This is the reference implementation kept for cross-checking.  It encodes
    every PairFlavour via _embed_compressed (exploitation 5a: sign symmetry in
    polynomial coefficients) but does NOT apply the 5b overlap-block dropping
    optimisation.  Use encode() for the full optimised encoding.

    Returns a 1-D numpy array of float64.  Each ASSOC contributes 2*pf.count()
    reals (the complex polynomial coefficients viewed as (re, im) pairs).
    If the plan has no registered operations, returns an empty array.
    """
    fo_list = repL(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    if not pf_list:
        return np.array([], dtype=float)

    seen_keys = set()
    group_sizes = tuple(g.size for g in plan.context.groups)
    parts = []
    for pf in pf_list:
        key = _encoding_canonical_key(pf)
        if key in seen_keys:
            print("BEWARE: a deduplication event fired in encode.py despite it being thought that this could/should never happen with the implementation at time of coding.")
            continue
        seen_keys.add(key)
        z_pos = np.array(eval_pair_orbit_positive(pf, plan, event), dtype=complex)
        assert len(z_pos) == pf.count(group_sizes), (
            f"eval_pair_orbit_positive returned {len(z_pos)} values "
            f"but expected pf.count={pf.count(group_sizes)} for {pf!r}"
        )
        parts.append(_complex_to_reals(_embed_compressed(z_pos, pf)))

    return np.concatenate(parts)


def encode(plan, event: dict) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector (optimised: 5a + 5b).

    Two-phase encoding
    ------------------
    Phase 1 — single-orbit sort encoding:
        For each distinct FlavouredOperator fo in repL(plan), sort the real
        values {eval(u, event) : u in fo.atoms()} and pack adjacent pairs as
        complex numbers.  Each fo contributes ceil(fo.count() / 2) complex
        values.  A fo is stored exactly once even if it appears as the u-side
        or v-side of many PairFlavours.

    Phase 2 — compressed pair encoding with largest association dropped:
        PairFlavours from canonical_pair_flavours are grouped into OVERLAP
        BLOCKS: consecutive PairFlavours sharing (op_u, flavour_u, op_v,
        flavour_v).  Within each block the association (PairFlavour) with the
        largest pf.count() is omitted — its values are deducible from the
        Phase 1 orbit headers plus the other stored associations in the block
        (by complementation in the association table).  All remaining
        associations are encoded via _embed_compressed (5a compression).

    Information completeness
    ------------------------
    The dropped association's z-values {eval(u_k)+i*eval(v_k)} can be
    recovered as the complement of the stored associations' (real, imag) pairs
    within the full Cartesian product of Phase-1 u-values × Phase-1 v-values.
    The complementation operates on multisets, so repeated eval values (e.g.
    two particles with identical position vectors) are handled correctly: their
    z-values appear with the right multiplicity in the complement, and the
    polynomial embedding captures that multiplicity exactly.  No information
    about the event is lost for any event — degenerate or otherwise.  (Particles
    that are numerically identical are correctly represented as indistinguishable,
    which is the desired behaviour for a permutation-invariant encoding.)

    Output structure
    ----------------
    [Phase 1 ORBIT segments ... | Phase 2 ASSOC segments ...]
    Use describe_encoding(plan) for a segment-by-segment map.

    Returns a 1-D numpy array of float64.  The output is entirely real:
    ORBIT segments contribute fo.count() reals (sorted eval values);
    ASSOC segments contribute 2*pf.count() reals (complex polynomial
    coefficients viewed as interleaved (re, im) pairs).
    """
    fo_list = repL(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    group_sizes = tuple(g.size for g in plan.context.groups)

    parts = []

    # Phase 1: sort-encode the orbit of each distinct FlavouredOperator.
    # fo is identified by (op.name, op.rank, op.parity, op.argument_symmetry, flavour.counts)
    # so that we don't depend on FlavouredOperator.__hash__ including Context.
    seen_fo_keys: set = set()
    for fo in fo_list:
        fo_key = (fo.operation.name, fo.operation.rank, fo.operation.parity,
                  fo.operation.argument_symmetry.value, fo.flavour.counts)
        if fo_key in seen_fo_keys:
            continue
        seen_fo_keys.add(fo_key)
        if fo.count() == 0:
            continue
        values = np.array(eval_single_orbit_compressed(fo, plan, event), dtype=float)
        parts.append(_sort_encode(values))

    # Phase 2: group PairFlavours into OVERLAP BLOCKS.
    # Drop order (both commute; self-pair first so complementation removes the
    # largest of what remains, not accidentally the self-pair):
    #   NULL_SELF: any self-pairing association (z_k = (1+i)*a_k, Phase 1 redundant)
    #   NULL_COMP: largest remaining association per block (complementation drop)
    for _block_key, block_iter in groupby(pf_list, key=_overlap_block_key):
        block = list(block_iter)
        non_self_idx = [i for i, pf in enumerate(block) if not _is_self_pair(pf)]
        comp_drop = (max(non_self_idx, key=lambda i: block[i].count(group_sizes))
                     if non_self_idx else None)
        for i, pf in enumerate(block):
            if _is_self_pair(pf):
                continue   # NULL_SELF: deducible from Phase 1
            if i == comp_drop:
                continue   # NULL_COMP: deducible by complementation
            z_pos = np.array(eval_pair_orbit_positive(pf, plan, event), dtype=complex)
            assert len(z_pos) == pf.count(group_sizes), (
                f"eval_pair_orbit_positive returned {len(z_pos)} values "
                f"but expected pf.count={pf.count(group_sizes)} for {pf!r}"
            )
            parts.append(_complex_to_reals(_embed_compressed(z_pos, pf)))

    return np.concatenate(parts) if parts else np.array([], dtype=float)


"""
Why 8-fold deduplication

    (+u,+v),(-u,+v),(-u,-v),(+u,-v),(+v,+u),(-v,+u),(-v,-u),(+v,-u)

might not be needed.

_encoding_canonical_key(pf) produces
    (tuple(sorted([(op_u.name, fl_u.counts), (op_v.name, fl_v.counts)])), overlap).

Two Operation objects with the same .name are equal by the Operation dataclass
equality (name, rank, parity, argument_symmetry all match). So different
operations → different names → different side-keys.
Two Flavour objects with the same .counts are equal. So different flavours →
different count-tuples → different side-keys.
Therefore: if two PairFlavours differ in any of (op_u, fl_u, op_v, fl_v, overlap),
their canonical keys also differ.
Since canonical_pair_flavours already returns a duplicate-free list, seen_keys
is injective over that list — it never fires.
The argument rests on the assumption that no two registered operations share a
name with different behaviour (which would be a user error). Under normal use,
the deduplication is provably redundant given the current symatom machinery.
The one hole: if someone registers two different EvaluableOperation objects with
the same name and different eval_fn, those would be considered equal by symatom
(which is intentional — eval_fn is excluded from equality) but would also produce
the same side-key. In that case seen_keys would suppress a PairFlavour
incorrectly. But that's a pre-existing design constraint on how operations are
registered, not a new risk introduced here.
"""
