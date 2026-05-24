"""
symcoder.tests.test_alignment_decoder_adversarial
===================================================
Tests using a deliberately adversarial event that exposes a known limitation
of the current binary-pair-only alignment decoder.

Background
----------
The alignment decoder (decode_alignment) resolves eps(a,b,c)'s column values
by consulting binary pair relations: for each already-placed atom ref, it looks
up pair_rel[(ref, eps(a,b,c))] and intersects the candidate v-values.

For *generic* events this works because each pair relation typically narrows
candidates to a single value.  For the adversarial event below, the special
structure of p and q makes EVERY binary pair 2-to-1 for eps(a,b,c):

  * p = (0,0,+1) and q = (0,0,-1) are antipodal unit vectors.
  * Consequence 1: eps(a,p,q) = a·(p×q) = a·0 = 0  for every electron a.
    The eps(a,p,q) row is identically zero → no discriminating information.
  * Consequence 2: dot_2D(x,p) = x[:2]·p[:2] = 0  for every particle x
    (since p[:2] = (0,0)).  All euclidean2 products with p are zero.
  * Consequence 3: every eps(a,b,p) value appears TWICE in the alignment table
    — once from some (σ_even, τ) column and once from (σ_odd', τ') — carrying
    both signs of eps(a,b,c).  No single eps(a,b,p) value determines the sign.

The sign of eps(a,b,c) for each column IS uniquely determined by the JOINT
pair (dot3D(a,p) = u, eps(a,b,p) = w) — this pair identifies the group element
exactly (see test_joint_constraint_resolves_eps_abc_sign).  But the current
decoder only intersects BINARY constraints and falls back to "take smallest
candidate" when both ±0.643 remain, assigning the wrong sign to 6 of 12
columns.

Current status
--------------
test_alignment_decoder_adversarial_event_fails    — EXPECTED TO FAIL.
    Documents the known limitation.  When the decoder is fixed to use joint
    multi-atom constraints, this test should PASS and the xfail marker removed.

test_binary_constraints_are_ambiguous_for_eps_abc — PASSES.
    Proves that no binary pair relation constrains eps(a,b,c)'s sign:
    for every placed atom ref, pair_rel[(ref, eps_abc)] always returns
    both +v and -v for every ref value.

test_joint_constraint_resolves_eps_abc_sign        — PASSES.
    Proves that the COMBINATION (dot3D(a,p), eps(a,b,p)) uniquely determines
    eps(a,b,c) for every column.  This is the information the fixed decoder
    should exploit.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pytest
from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan, repS as _repS_fn
import symcoder.operations.euclidean2
import symcoder.operations.euclidean3
from symcoder.alignment_decoder import decode_alignment, _lookup_v_vals
from symcoder.encoders import (
    OrbitEncoderFactory, SortEncoderFactory, HalfSortEncoderFactory,
    standard_row_pair_factories, OverlapBlockEncoderFactory, Phase2EncoderFactory,
)
from symcoder.eval import evaluate


# ---------------------------------------------------------------------------
# Shared setup
# ---------------------------------------------------------------------------

def _make_plan():
    """Return (plan, ctx) for the 3-electron 2-muon demo plan."""
    mymag = Operation(
        "Mymag", rank=1, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))),
    )
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    plan      = Plan(context=ctx, operations=(
        mymag,
        symcoder.operations.euclidean3.dot,
        symcoder.operations.euclidean2.dot,
        symcoder.operations.euclidean3.eps,
    ))
    return plan, ctx


# Adversarial event: p and q are antipodal unit vectors along z.
# This makes eps(a,p,q) = 0 for every electron, and all 2D dot-products
# with p or q equal zero.
ADVERSARIAL_EVENT = {
    "a": np.array([ 1.2,  0.3,  0.5]),
    "b": np.array([ 0.4,  1.1,  0.2]),
    "c": np.array([-0.5,  1.2, -0.7]),
    "p": np.array([ 0.0,  0.0, +1.0]),
    "q": np.array([ 0.0,  0.0, -1.0]),
}


def _encode_and_decode_phase1_and_phase2(plan, event):
    """Run the full encode/decode pipeline; return (phase1_results, all_pair_decoded)."""
    orbit_fac  = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
    phase2_fac = Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories(),
                                   use_complementarity_drop=False)
    ])
    orbit_enc  = orbit_fac.build(plan)
    phase2_enc = phase2_fac.build(plan)
    ph1_vals   = orbit_enc.encode(event).values
    ph2_vals   = phase2_enc.encode(event).values

    phase1_results = {}
    cursor = 0
    for fo, enc in orbit_enc._selections:
        key    = (fo.operation, tuple(fo.flavour.counts))
        chunk  = ph1_vals[cursor : cursor + enc.output_dim]
        phase1_results[key] = enc.decode(chunk)
        cursor += enc.output_dim

    all_pair_decoded = []
    p2c = 0
    for _, benc in phase2_enc._block_encoders:
        bs = ph2_vals[p2c : p2c + benc.output_dim]
        all_pair_decoded.extend(benc.decode(bs, phase1_results))
        p2c += benc.output_dim

    return phase1_results, all_pair_decoded


def _build_pair_rel(plan, all_pair_decoded):
    """Reproduce the pair_rel dict from decode_alignment internals."""
    the_group = plan.context.the_group
    n_cols    = the_group.order()
    fo_list   = _repS_fn(plan.context, plan.operations)
    repS_atoms    = [fo.canonical_representative() for fo in fo_list]
    repS_atom_set = set(repS_atoms)

    pair_rel = {}
    for decoded in all_pair_decoded:
        if not decoded.atom_pairs or not decoded.pairs:
            continue
        u0, v0 = decoded.atom_pairs[0]
        if u0 == v0 or u0 not in repS_atom_set or v0 not in repS_atom_set:
            continue
        stab = the_group.stabiliser_size_pair(u0, v0)
        fp   = list(decoded.pairs) * stab
        if len(fp) != n_cols:
            continue
        if (u0, v0) not in pair_rel:
            pair_rel[(u0, v0)] = fp
        if (v0, u0) not in pair_rel:
            pair_rel[(v0, u0)] = [(y, x) for x, y in fp]
    return pair_rel, repS_atoms


# ---------------------------------------------------------------------------
# Test 1: document the current decoder failure (expected to fail)
# ---------------------------------------------------------------------------

@pytest.mark.xfail(
    reason=(
        "Known limitation: the binary-pair alignment decoder cannot recover "
        "eps(a,b,c)'s sign for the adversarial event where p=(0,0,1) and "
        "q=(0,0,-1).  Every binary pair relation is 2-to-1 for eps(a,b,c), so "
        "the decoder defaults to assigning -|v| to all 12 columns.  The fix "
        "requires the decoder to use joint multi-atom constraints (e.g. the "
        "combination of dot3D(a,p) and eps(a,b,p) together uniquely determines "
        "eps(a,b,c)'s sign).  See test_joint_constraint_resolves_eps_abc_sign."
    ),
    strict=True,   # must fail; if it passes, the marker should be removed
)
def test_alignment_decoder_adversarial_event_fails():
    """
    The current alignment decoder produces a MISMATCH on the adversarial event.

    Ground truth: 6 columns with eps(a,b,c) = +0.643, 6 with -0.643.
    Decoder output: all 12 columns assigned -0.643 (wrong sign for 6).

    This test is marked xfail(strict=True): it is expected to fail today.
    When the decoder is fixed to handle joint constraints, this test will pass
    and the xfail marker should be removed.
    """
    plan, ctx     = _make_plan()
    event         = ADVERSARIAL_EVENT
    atol          = 1e-8

    phase1_results, all_pair_decoded = _encode_and_decode_phase1_and_phase2(
        plan, event
    )
    decoded = decode_alignment(plan, phase1_results, all_pair_decoded, atol=atol)

    the_group  = plan.context.the_group
    n_cols     = the_group.order()
    fo_list    = _repS_fn(plan.context, plan.operations)
    repS_atoms = [fo.canonical_representative() for fo in fo_list]

    gt_vectors = [
        tuple(evaluate(g.apply(a), event) for a in repS_atoms)
        for g in the_group.all_group_elements()
    ]

    dec_sorted = sorted(decoded.vectors)
    gt_sorted  = sorted(gt_vectors)
    ok = (
        len(dec_sorted) == len(gt_sorted) and
        all(
            max(abs(a - b) for a, b in zip(da, ga)) <= atol
            for da, ga in zip(dec_sorted, gt_sorted)
        )
    )
    assert ok, (
        f"Expected alignment decoder to produce a MISMATCH on the adversarial "
        f"event (this is a known bug).  If it now passes, the xfail marker on "
        f"this test should be removed."
    )


# ---------------------------------------------------------------------------
# Test 2: prove binary constraints are 2-to-1 for eps(a,b,c)
# ---------------------------------------------------------------------------

def test_binary_constraints_are_ambiguous_for_eps_abc():
    """
    For the adversarial event, every binary pair relation involving eps(a,b,c)
    as the 'v' atom is 2-to-1: for every distinct 'u' value, both +|v| and -|v|
    appear in the pair list.

    This formally proves that the current binary-pair decoder cannot determine
    eps(a,b,c)'s column sign from any single pair relation in isolation.
    """
    plan, ctx     = _make_plan()
    event         = ADVERSARIAL_EVENT

    phase1_results, all_pair_decoded = _encode_and_decode_phase1_and_phase2(
        plan, event
    )
    pair_rel, repS_atoms = _build_pair_rel(plan, all_pair_decoded)

    eps_abc = next(
        a for a in repS_atoms
        if "eps" in str(a.operation.name).lower() and "c" in a.labels
    )

    # Collect all pair relations where eps(a,b,c) is the v-atom.
    pairs_to_eps_abc = {
        ref: pair_rel[(ref, eps_abc)]
        for ref in repS_atoms
        if (ref, eps_abc) in pair_rel
    }

    assert pairs_to_eps_abc, "No pair relation found for eps(a,b,c) — encoding gap."

    for ref, pairs in pairs_to_eps_abc.items():
        by_u: dict = {}
        for u, v in pairs:
            by_u.setdefault(round(u, 6), []).append(round(v, 6))

        for u_val, v_vals in by_u.items():
            distinct_v = set(v_vals)
            assert len(distinct_v) >= 2, (
                f"Found a binary constraint that IS unambiguous:\n"
                f"  pair ({ref}, eps(a,b,c)):  u={u_val} → v={distinct_v}\n"
                f"  This would resolve the sign; the decoder should use it."
            )


# ---------------------------------------------------------------------------
# Test 3: prove the JOINT pair (dot3D(a,p), eps(a,b,p)) resolves eps(a,b,c)
# ---------------------------------------------------------------------------

def test_joint_constraint_resolves_eps_abc_sign():
    """
    The combination of dot3D(a,p) and eps(a,b,p) uniquely determines
    eps(a,b,c)'s sign: no two columns share the same (dot3D(a,p), eps(a,b,p))
    value pair.

    This proves the information IS available in the encoded data.  The decoder
    just needs to exploit it — i.e. consult which column the pair (ref_u, ref_w)
    came from, rather than intersecting independent binary lookups.
    """
    plan, ctx     = _make_plan()
    event         = ADVERSARIAL_EVENT

    the_group  = plan.context.the_group
    fo_list    = _repS_fn(plan.context, plan.operations)
    repS_atoms = [fo.canonical_representative() for fo in fo_list]

    eps_abc = next(
        a for a in repS_atoms
        if "eps" in str(a.operation.name).lower() and "c" in a.labels
    )
    eps_abp = next(
        a for a in repS_atoms
        if "eps" in str(a.operation.name).lower()
        and "p" in a.labels and "c" not in a.labels
    )
    dot3_ap = next(
        a for a in repS_atoms
        if a.operation is symcoder.operations.euclidean3.dot
        and "p" in a.labels
    )

    # Build the ground-truth table (no decoding — just direct evaluation).
    gt_table: dict = {}   # (round(dot3_ap), round(eps_abp)) → round(eps_abc)
    for g in the_group.all_group_elements():
        d_val = round(evaluate(g.apply(dot3_ap), event), 6)
        e_val = round(evaluate(g.apply(eps_abp), event), 6)
        c_val = round(evaluate(g.apply(eps_abc), event), 6)
        joint_key = (d_val, e_val)
        assert joint_key not in gt_table or gt_table[joint_key] == c_val, (
            f"Joint key {joint_key} maps to two different eps(a,b,c) values: "
            f"{gt_table[joint_key]} and {c_val}.  Constraint is NOT 1-to-1."
        )
        gt_table[joint_key] = c_val

    # Every joint key should map to exactly ONE eps(a,b,c) value.
    n_keys = len(gt_table)
    n_cols = the_group.order()
    assert n_keys == n_cols, (
        f"Expected {n_cols} distinct (dot3D(a,p), eps(a,b,p)) pairs (one per "
        f"group element), got {n_keys}.  Joint constraint is NOT 1-to-1."
    )

    # Explicitly check that the two eps(a,b,c) values (±|v|) each appear exactly 6 times.
    eps_vals = sorted(gt_table.values())
    counter  = Counter(eps_vals)
    n_distinct = len(counter)
    assert n_distinct == 2, (
        f"Expected exactly 2 distinct eps(a,b,c) values (±|v|), got {n_distinct}: "
        f"{counter}"
    )
    for val, cnt in counter.items():
        assert cnt == n_cols // 2, (
            f"eps(a,b,c) = {val} appears {cnt} times, expected {n_cols // 2}."
        )


# ---------------------------------------------------------------------------
# Test 4: adversarial event passes with a generic (random) event
# ---------------------------------------------------------------------------

def test_alignment_decoder_generic_event_passes():
    """
    Sanity check: the alignment decoder works correctly for a generic random
    event on the same 3+2 plan.  Confirms the bug is specific to the adversarial
    event structure, not a systemic regression.
    """
    plan, ctx = _make_plan()
    rng       = np.random.default_rng(42)
    event     = {lbl: rng.standard_normal(3) for lbl in ctx.all_labels}
    atol      = 1e-8

    phase1_results, all_pair_decoded = _encode_and_decode_phase1_and_phase2(
        plan, event
    )
    decoded = decode_alignment(plan, phase1_results, all_pair_decoded, atol=atol)

    the_group  = plan.context.the_group
    fo_list    = _repS_fn(plan.context, plan.operations)
    repS_atoms = [fo.canonical_representative() for fo in fo_list]

    gt_vectors = [
        tuple(evaluate(g.apply(a), event) for a in repS_atoms)
        for g in the_group.all_group_elements()
    ]

    dec_sorted = sorted(decoded.vectors)
    gt_sorted  = sorted(gt_vectors)
    assert len(dec_sorted) == len(gt_sorted)
    for da, ga in zip(dec_sorted, gt_sorted):
        assert max(abs(a - b) for a, b in zip(da, ga)) <= atol, (
            f"Decoder mismatch on generic random event.\n"
            f"Decoded: {da}\nGround truth: {ga}"
        )
