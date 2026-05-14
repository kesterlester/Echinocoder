from __future__ import annotations
import importlib.util
import numpy as np
from itertools import groupby
from pathlib import Path
from symatom.atoms import ArgumentSymmetry
from symatom.group import SignCorrelationType
from symatom.rep import canonical_pair_flavours
from symatom import repS
## from .pairs import eval_pair_orbit, eval_pair_orbit_positive, eval_single_orbit, eval_single_orbit_compressed, _is_self_pair
from .pairs import eval_pair_orbit, eval_pair_orbit_positive, _is_self_pair
from .encoders import AtomOrbitEncoderRegistry, OrbitSpec
from .describe import SegmentInfo, _orbit_example, _assoc_example, _symmetry_class

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
    atoms in the orbit (fo.count_of_atoms_one_per_sign()), one real per atom.

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


def encode_described(
    plan,
    event: dict,
    registry: AtomOrbitEncoderRegistry | None = None,
) -> tuple[np.ndarray, list[SegmentInfo]]:
    """
    Core encoding function.  Runs the full two-phase encoding and simultaneously
    collects a SegmentInfo descriptor for every segment produced.

    Returns (values, segments) where:
      values   — 1-D float64 numpy array (the permutation-invariant embedding)
      segments — list[SegmentInfo] mirroring the structure of values exactly

    Phase 1 metadata (method_name, output_dim) is sourced from the registry's
    assess() responses when a registry is supplied, making describe_encoding()
    automatically consistent with encode() — both use the same selection logic.

    Phase 2 is currently always the algebraic _embed_compressed path (no registry
    connection yet); its SegmentInfo fields are computed algebraically.

    Note: the diagnostic print in Phase 1 is temporary scaffolding and will be
    gated by a bool flag in a future performance pass.
    """
    fo_list    = repS(plan.context, plan.operations)
    pf_list    = canonical_pair_flavours(fo_list, plan.context)
    type_sizes = tuple(g.size for g in plan.context.types)
    types      = plan.context.types

    parts    = []
    segments = []
    cursor   = 0

    # ------------------------------------------------------------------
    # Phase 1: encode each FlavouredOperator orbit.
    # ------------------------------------------------------------------
    for fo in fo_list:
        assert fo.count_of_atoms_one_per_sign() > 0, "No FlavouredOperator generated by repS should have an empty orbit."

        if registry is not None:
            spec   = OrbitSpec.from_flavoured_operator(fo)
            capable = registry.query_all(spec, plan)
            # Diagnostic: show what each capable encoder claims it would produce.
            # TODO: gate this behind a bool flag when performance matters.
            for enc, cap in capable:
                print(f"  [{type(enc).__name__}] proposes  method={cap.method_name}  "
                      f"output_dim={cap.output_dim}  priority={cap.priority}"
                      f"  fo={fo}")
            chosen_enc, chosen_cap = capable[0]
            length = chosen_cap.output_dim   # authoritative from assess()
            # OLD WRONG: result = chosen_enc.encode(spec, event, plan)
            result = chosen_enc.encode(chosen_cap, event, plan)
            assert len(result.values) == length # Checks contract was met!
            parts.append(result.values)
            segments.append(SegmentInfo(
                kind        = "ORBIT",
                start       = cursor,
                length      = length,
                op_u        = fo.operation.name,
                flavour_u   = tuple(fo.flavour.counts),
                method_name = chosen_cap.method_name,
                example     = _orbit_example(fo.operation.name, fo.flavour.counts, types),
            ))
        else:
            raise NotImplementedError("We no longer support evals other than by the encoder registry.") # TODO .. DELETE THIS BRANCH ALTOGETHER when confident we don't need to put it back.
            ## raw    = np.array(eval_single_orbit_compressed(fo, plan, event), dtype=float)
            ## sorted_vals = _sort_encode(raw)
            ## parts.append(sorted_vals)
            ## length = len(sorted_vals)
            ## antisym = fo.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
            ## segments.append(SegmentInfo(
            ##     kind            = "ORBIT",
            ##     start           = cursor,
            ##     length          = length,
            ##     op_u            = fo.operation.name,
            ##     flavour_u       = tuple(fo.flavour.counts),
            ##     sign_compressed = antisym,
            ##     example         = _orbit_example(fo.operation.name, fo.flavour.counts, types),
            ## ))

        cursor += length

    # ------------------------------------------------------------------
    # Phase 2: compressed pair encoding with largest association dropped.
    # Drop order: NULL_SELF first (Phase-1 redundant), then NULL_COMP
    # (complementation drop of the largest remaining per block).
    # ------------------------------------------------------------------
    for _bk, block_iter in groupby(pf_list, key=_overlap_block_key):
        block        = list(block_iter)
        non_self_idx = [i for i, pf in enumerate(block) if not _is_self_pair(pf)]
        comp_drop    = (max(non_self_idx, key=lambda i: block[i].count(type_sizes))
                        if non_self_idx else None)

        for i, pf in enumerate(block):
            sc       = _symmetry_class(pf)
            ex       = _assoc_example(
                pf.op_u.name, pf.flavour_u.counts,
                pf.op_v.name, pf.flavour_v.counts,
                pf.overlap, types,
            )
            notional = 2 * pf.count(type_sizes)
            common   = dict(
                op_u            = pf.op_u.name,
                flavour_u       = tuple(pf.flavour_u.counts),
                op_v            = pf.op_v.name,
                flavour_v       = tuple(pf.flavour_v.counts),
                overlap         = tuple(pf.overlap),
                symmetry_class  = sc,
                notional_length = notional,
                example         = ex,
            )
            if _is_self_pair(pf):
                segments.append(SegmentInfo(kind="NULL_SELF", start=cursor, length=0, **common))
            elif i == comp_drop:
                segments.append(SegmentInfo(kind="NULL_COMP", start=cursor, length=0, **common))
            else:
                z_pos = np.array(eval_pair_orbit_positive(pf, plan, event), dtype=complex)
                assert len(z_pos) == pf.count(type_sizes), (
                    f"eval_pair_orbit_positive returned {len(z_pos)} values "
                    f"but expected pf.count={pf.count(type_sizes)} for {pf!r}"
                )
                parts.append(_complex_to_reals(_embed_compressed(z_pos, pf)))
                segments.append(SegmentInfo(kind="ASSOC", start=cursor, length=notional, **common))
                cursor += notional

    values = np.concatenate(parts) if parts else np.array([], dtype=float)
    return values, segments


def encode(plan, event: dict, registry: AtomOrbitEncoderRegistry | None = None) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector.
    See encode_described() for full documentation of the two-phase algorithm.
    Returns a 1-D float64 numpy array.
    """
    values, _segments = encode_described(plan, event, registry)
    return values


def describe_encoding(
    plan,
    registry: AtomOrbitEncoderRegistry | None = None,
) -> list[SegmentInfo]:
    """
    Return a list of SegmentInfo objects describing the full structure of the
    vector produced by encode(plan, event, registry).

    When a registry is supplied, Phase 1 metadata (method_name, output_dim) is
    derived from the registry's assess() responses — the same selection logic
    used by encode() — so the description is always in sync with the encoding.
    Phase 2 remains algebraic for now.

    When no registry is supplied, falls back to the purely algebraic description
    (compatible with the pre-registry legacy behaviour).

    Usage
    -----
    for seg in describe_encoding(plan, registry):
        print(seg)

    import json
    json.dumps([s.to_dict() for s in describe_encoding(plan, registry)])

    sum(s.length for s in describe_encoding(plan, registry))  # == len(encode(...))
    """
    fo_list    = repS(plan.context, plan.operations)
    pf_list    = canonical_pair_flavours(fo_list, plan.context)
    type_sizes = tuple(g.size for g in plan.context.types)
    types      = plan.context.types

    segments = []
    cursor   = 0

    # Phase 1: one ORBIT segment per FlavouredOperator.
    for fo in fo_list:
        if fo.count_of_atoms_one_per_sign() == 0:
            continue

        if registry is not None:
            spec    = OrbitSpec.from_flavoured_operator(fo)
            capable = registry.query_all(spec, plan)
            if capable:
                chosen_enc, chosen_cap = capable[0]
                length = chosen_cap.output_dim
                segments.append(SegmentInfo(
                    kind        = "ORBIT",
                    start       = cursor,
                    length      = length,
                    op_u        = fo.operation.name,
                    flavour_u   = tuple(fo.flavour.counts),
                    method_name = chosen_cap.method_name,
                    example     = _orbit_example(fo.operation.name, fo.flavour.counts, types),
                ))
                cursor += length
            else:
                # No encoder capable — record with zero length so the segment
                # still appears in the description (length 0, method "FAILED").
                segments.append(SegmentInfo(
                    kind        = "ORBIT",
                    start       = cursor,
                    length      = 0,
                    op_u        = fo.operation.name,
                    flavour_u   = tuple(fo.flavour.counts),
                    method_name = "********* FAILED **********",
                    example     = _orbit_example(fo.operation.name, fo.flavour.counts, types),
                ))
        else:
            # Legacy algebraic path (no registry).
            antisym = fo.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
            n = fo.count_of_atoms_one_per_sign()
            segments.append(SegmentInfo(
                kind            = "ORBIT",
                start           = cursor,
                length          = n,
                notional_length = n,
                op_u            = fo.operation.name,
                flavour_u       = tuple(fo.flavour.counts),
                sign_compressed = antisym,
                example         = _orbit_example(fo.operation.name, fo.flavour.counts, types),
            ))
            cursor += n

    # Phase 2: algebraic (unchanged; no registry connection yet).
    for _bk, block_iter in groupby(pf_list, key=_overlap_block_key):
        block        = list(block_iter)
        non_self_idx = [i for i, pf in enumerate(block) if not _is_self_pair(pf)]
        comp_drop    = (max(non_self_idx, key=lambda i: block[i].count(type_sizes))
                        if non_self_idx else None)

        for i, pf in enumerate(block):
            sc       = _symmetry_class(pf)
            ex       = _assoc_example(
                pf.op_u.name, pf.flavour_u.counts,
                pf.op_v.name, pf.flavour_v.counts,
                pf.overlap, types,
            )
            notional = 2 * pf.count(type_sizes)
            common   = dict(
                op_u            = pf.op_u.name,
                flavour_u       = tuple(pf.flavour_u.counts),
                op_v            = pf.op_v.name,
                flavour_v       = tuple(pf.flavour_v.counts),
                overlap         = tuple(pf.overlap),
                symmetry_class  = sc,
                notional_length = notional,
                example         = ex,
            )
            if _is_self_pair(pf):
                segments.append(SegmentInfo(kind="NULL_SELF", start=cursor, length=0, **common))
            elif i == comp_drop:
                segments.append(SegmentInfo(kind="NULL_COMP", start=cursor, length=0, **common))
            else:
                segments.append(SegmentInfo(kind="ASSOC", start=cursor, length=notional, **common))
                cursor += notional

    return segments


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
