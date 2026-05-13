from __future__ import annotations
import importlib.util
import numpy as np
from itertools import groupby
from pathlib import Path
from symatom.atoms import ArgumentSymmetry
from symatom.group import SignCorrelationType
from symatom.rep import canonical_pair_flavours
from symatom import repS
from .pairs import eval_pair_orbit, eval_pair_orbit_positive, eval_single_orbit, eval_single_orbit_compressed, _is_self_pair
from .encoders import AtomOrbitEncoderRegistry, OrbitSpec

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


def _sign_correlation_type_from_pf(pf) -> SignCorrelationType:
    """
    Compute the SignCorrelationType for a PairFlavour from its structural
    parameters alone (no specific atom instances required).

    This is equivalent to TheGroup.sign_correlation_type(u, v) for any
    representative (u, v) of this PairFlavour, since the type depends only
    on the per-group label counts, overlap, and operation antisymmetry —
    all of which are encoded in the PairFlavour.

    The u_g == v_g condition (per-group label sets identical) is equivalent
    to full per-group overlap: overlap_g == ku_g == kv_g.

    Algorithm: per-group achievable-set product (same as TheGroup.sign_correlation_type).
    """
    antisym_u = pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
    antisym_v = pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

    achievable: set[tuple[int, int]] = {(1, 1)}

    for ku_g, kv_g, s_g in zip(
        pf.flavour_u.counts, pf.flavour_v.counts, pf.overlap
    ):
        pu_can_flip = antisym_u and ku_g >= 2
        pv_can_flip = antisym_v and kv_g >= 2
        u_g_eq_v_g  = (s_g == ku_g == kv_g)   # full per-group overlap ↔ same label set

        if not pu_can_flip and not pv_can_flip:
            g_ach: set[tuple[int, int]] = {(1, 1)}
        elif pu_can_flip and not pv_can_flip:
            g_ach = {(1, 1), (-1, 1)}
        elif not pu_can_flip and pv_can_flip:
            g_ach = {(1, 1), (1, -1)}
        elif u_g_eq_v_g:
            g_ach = {(1, 1), (-1, -1)}
        else:
            g_ach = {(1, 1), (1, -1), (-1, 1), (-1, -1)}

        achievable = {
            (pu * qu, pv * qv)
            for (pu, pv) in achievable
            for (qu, qv) in g_ach
        }

    u_flips   = (-1,  1) in achievable
    v_flips   = ( 1, -1) in achievable
    both_flip = (-1, -1) in achievable

    if u_flips and v_flips:
        return SignCorrelationType.TYPE_22
    elif u_flips:
        return SignCorrelationType.TYPE_21
    elif v_flips:
        return SignCorrelationType.TYPE_12
    elif both_flip:
        return SignCorrelationType.TYPE_NEG
    else:
        return SignCorrelationType.TYPE_11


def _embed_compressed(z_pos: np.ndarray, pf) -> np.ndarray:
    """
    Embed the sign=(+1,+1) orbit evaluations z_pos into a permutation-invariant
    vector of exactly pf.count(...) complex numbers, exploiting the polynomial
    symmetry structure imposed by the G-orbit's sign correlation type.

    Background
    ----------
    The characteristic polynomial of the full orbit multiset {z_k} has forced
    structure depending on the SignCorrelationType of the pair.  The type is
    determined by _sign_correlation_type_from_pf(pf), NOT by the raw
    ArgumentSymmetry flags of the operations.

    CRITICAL: ArgumentSymmetry alone does NOT determine the sign correlation
    type.  A cross-group ANTISYMMETRIC atom (≤1 label per group) has a
    permanently +1 sign and gives TYPE_11 rather than TYPE_12/21.  Similarly
    AA pairs with full per-group overlap give TYPE_NEG, not TYPE_22.  Always
    use _sign_correlation_type_from_pf to determine the correct branch.

    Why the sign=(+1,+1) subset
    ---------------------------
    For each base label combination, the Atom constructor sorts labels and
    absorbs permutation parity into the stored sign.  Exactly one of the
    (±u, ±v) canonical sign variants has stored (u.sign=+1, v.sign=+1) —
    so z_pos always has exactly n = pf.count() elements.

    Embedder coefficient convention
    --------------------------------
    The Echinocoder zip embedder computes
        polyfromroots(-data)[-2::-1]
    i.e. the polynomial with roots at -data, coefficients in DESCENDING degree
    order with the monic leading term dropped.  For n input values the output
    has n complex entries.  For a degree-2n polynomial the output has 2n entries;
    entry i has degree-parity (2n-1-i) mod 2: i even → odd degree, i odd → even.

    Algorithm per SignCorrelationType
    ----------------------------------
    TYPE_11  (achievable = {(+1,+1)} — no sign changes):
        Orbit: {z_k}  (n elements)
        Embed z_pos directly (n roots → n complex coefficients).

    TYPE_NEG (achievable = {(+1,+1),(−1,−1)} — correlated negation only):
        Orbit: {z_k, −z_k}  (2n elements, closed under z → −z)
        Invariant: {z_k²}   (n complex, since (−z_k)² = z_k²)
        Embed z_pos² directly (n roots → n complex coefficients).
        Note: forming {z_k², conj(z_k²)} (as in TYPE_22) would be wrong here —
        conj(z_k²) is NOT in the orbit, and adding it loses information about
        Im(z_k²) = 2·Re(z_k)·Im(z_k).

    TYPE_12  (achievable = {(+1,+1),(+1,−1)} — v-sign flips freely):
        Orbit: {z_k, conj(z_k)}  (2n elements, conjugate-closed)
        Invariant multiset: {z_k, conj(z_k)}  (conjugate-closed, 2n roots)
        → polynomial has real coefficients → coeffs[i].imag ≈ 0
        Pack n independent reals as n complex: coeffs.real[0::2] + 1j*coeffs.real[1::2]

    TYPE_21  (achievable = {(+1,+1),(−1,+1)} — u-sign flips freely):
        Orbit: {z_k, −conj(z_k)}  (2n elements, anti-conjugate-closed)
        Invariant multiset: {z_k, −conj(z_k)}
        → even-degree coefficients are real, odd-degree are purely imaginary
        Pack: coeffs.imag[0::2] + 1j*coeffs.real[1::2]

    TYPE_22  (achievable = {all four sign combinations}):
        Orbit: {z_k, conj(z_k), −z_k, −conj(z_k)}  (4n elements)
        Invariant: {z_k², conj(z_k²)}  (2n, conjugate-closed)
        → polynomial has real coefficients
        Pack: w_coeffs.real[0::2] + 1j*w_coeffs.real[1::2]

    All branches output n complex numbers (= 2n reals when combined via
    _complex_to_reals), so the total output size is the same as SYM×SYM.

    Validation
    ----------
    Cross-check tests in test_encode.py verify permutation invariance and
    that the output length equals pf.count() for all symmetry classes.
    """
    sct = _sign_correlation_type_from_pf(pf)

    if sct == SignCorrelationType.TYPE_11:
        # No sign changes: embed z_pos directly.
        coeffs, _, _ = _zip_embed(z_pos)
        return coeffs

    elif sct == SignCorrelationType.TYPE_NEG:
        # Correlated negation: invariant is {z_k²}, embed directly.
        w = z_pos ** 2
        coeffs, _, _ = _zip_embed(w)
        return coeffs

    elif sct == SignCorrelationType.TYPE_12:
        # v-sign flips: conjugate-closed multiset → real polynomial.
        z_full = np.concatenate([z_pos, np.conj(z_pos)])
        coeffs, _, _ = _zip_embed(z_full)
        r = coeffs.real   # all imaginary parts ≈ 0
        return r[0::2] + 1j * r[1::2]

    elif sct == SignCorrelationType.TYPE_21:
        # u-sign flips: anti-conjugate-closed → even-real/odd-imag polynomial.
        z_full = np.concatenate([z_pos, -np.conj(z_pos)])
        coeffs, _, _ = _zip_embed(z_full)
        # even indices → odd-degree (purely imaginary); odd indices → even-degree (real)
        return coeffs.imag[0::2] + 1j * coeffs.real[1::2]

    else:
        # TYPE_22: both signs flip freely: squared conjugate-closed → real polynomial.
        assert sct == SignCorrelationType.TYPE_22
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
    fo_list = repS(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    if not pf_list:
        return np.array([], dtype=float)

    seen_keys = set()
    type_sizes = tuple(g.size for g in plan.context.types)
    parts = []
    for pf in pf_list:
        key = _encoding_canonical_key(pf)
        if key in seen_keys:
            print("BEWARE: a deduplication event fired in encode.py despite it being thought that this could/should never happen with the implementation at time of coding.")
            continue
        seen_keys.add(key)
        z_pos = np.array(eval_pair_orbit_positive(pf, plan, event), dtype=complex)
        assert len(z_pos) == pf.count(type_sizes), (
            f"eval_pair_orbit_positive returned {len(z_pos)} values "
            f"but expected pf.count={pf.count(type_sizes)} for {pf!r}"
        )
        parts.append(_complex_to_reals(_embed_compressed(z_pos, pf)))

    return np.concatenate(parts)


def encode(plan, event: dict, registry: AtomOrbitEncoderRegistry | None = None) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector (optimised: 5a + 5b).

    Two-phase encoding
    ------------------
    Phase 1 — single-orbit sort encoding:
        For each distinct FlavouredOperator fo in repS(plan), sort the real
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
    fo_list = repS(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    type_sizes = tuple(g.size for g in plan.context.types)

    parts = []

    # Phase 1: sort-encode the orbit of each distinct FlavouredOperator.
    for fo in fo_list:
        assert fo.count()>0, "No FlavouredOperator generaed by repS should have an empty orbit."
        if registry is not None:
            spec = OrbitSpec.from_flavoured_operator(fo)
            capable = registry.query_all(spec, plan)
            # Diagnostic: show what each capable encoder claims it would produce
            for enc, cap in capable:
                print(f"  [{type(enc).__name__}] proposes  method={cap.method_name}  "
                      f"output_dim={cap.output_dim}  priority={cap.priority} when given {spec} in {plan}")
            # Use the first capable encoder in the list
            first_enc, first_cap = capable[0]
            parts.append(first_enc.encode(spec, event, plan).values)
        else:
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
        comp_drop = (max(non_self_idx, key=lambda i: block[i].count(type_sizes))
                     if non_self_idx else None)
        for i, pf in enumerate(block):
            if _is_self_pair(pf):
                continue   # NULL_SELF: deducible from Phase 1
            if i == comp_drop:
                continue   # NULL_COMP: deducible by complementation
            z_pos = np.array(eval_pair_orbit_positive(pf, plan, event), dtype=complex)
            assert len(z_pos) == pf.count(type_sizes), (
                f"eval_pair_orbit_positive returned {len(z_pos)} values "
                f"but expected pf.count={pf.count(type_sizes)} for {pf!r}"
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
