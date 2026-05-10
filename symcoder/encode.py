from __future__ import annotations
import importlib.util
import numpy as np
from pathlib import Path
from symatom.atoms import ArgumentSymmetry
from symatom.rep import canonical_pair_flavours
from symatom import repL
from .pairs import eval_pair_orbit, eval_pair_orbit_positive

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


def encode(plan, event: dict) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector.

    Algorithm:
      1. Enumerate all canonical PairFlavours for this plan.
      2. Deduplicate under encoding equivalence (seen_keys check).
      3. For each kept PairFlavour, evaluate the sign=(+1,+1) orbit subset:
         z_k = eval(u'_k, E) + i * eval(v'_k, E)  (pf.count() values).
      4. Form the symmetry-class-appropriate invariant multiset and embed via
         _embed_compressed, yielding exactly pf.count() complex coefficients.
      5. Concatenate all parts.

    The output length is sum(pf.count(group_sizes) for pf in kept PairFlavours),
    which equals sum(pf.orbit_size(group_sizes) / symmetry_factor) — a factor
    2 smaller than the naive orbit embedding for pairs involving one ANTISYMMETRIC
    operation, and factor 4 smaller for ANTISYMMETRIC×ANTISYMMETRIC pairs.

    Returns a 1-D numpy array of complex dtype.
    If the plan has no registered operations, returns an empty array.
    """
    fo_list = repL(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    if not pf_list:
        return np.array([], dtype=complex)

    # Encoding-level deduplication: skip PairFlavours whose z-vectors carry
    # the same information as one already embedded (related by u↔v swap or
    # negation of a row for ANTISYMMETRIC operations).
    #
    # NOTE: It is believed — but not yet formally proved — that this check
    # never fires given the current symatom machinery.  The argument is:
    # _encoding_canonical_key is injective on PairFlavours (same-named
    # operations are equal, same-counted Flavours are equal), and
    # canonical_pair_flavours already returns a duplicate-free list — so
    # seen_keys would never contain a collision.  The check is retained
    # because (a) the proof has not been written down rigorously, and
    # (b) if the upstream generation ever changes, this safety net matters.
    #
    # Longer exposition in docstring at end of file.

    group_sizes = tuple(g.size for g in plan.context.groups)
    seen_keys = set()
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
        parts.append(_embed_compressed(z_pos, pf))

    return np.concatenate(parts)


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
