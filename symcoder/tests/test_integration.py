"""
symcoder.tests.test_integration
================================
End-to-end integration tests exercising the public encode_and_describe() API
and verifying that the EncodingTree segment metadata correctly identifies every
slice of the flat encoded vector.

These tests sit "above" the unit tests in test_decoders.py and test_encode.py:
they call the top-level public API and cross-check the outputs against each
other and against direct atom evaluations, without inspecting internal encoder
state beyond what the tree exposes.

Each test runs in both float and integer event modes via the parametrised
``event_mode`` fixture in conftest.py.
"""
from __future__ import annotations

import numpy as np
import pytest
from symatom import Plan

from symcoder import encode, encode_and_describe, describe_encoding, EvaluableOperation
from symcoder.eval import evaluate
from symcoder.encoders.sort_encoder import SortEncoder, HalfSortEncoder


# ---------------------------------------------------------------------------
# Helper: ground-truth Phase 1 slice for one orbit encoder
# ---------------------------------------------------------------------------

def _phase1_ground_truth(enc, event):
    """Return the sorted values that the orbit encoder would/did store.

    HalfSortEncoder stores sorted absolute values of the positive-sign reps.
    SortEncoder      stores sorted evaluations of the full orbit.
    """
    if isinstance(enc, HalfSortEncoder):
        return sorted(abs(evaluate(a, event)) for a in enc._representatives)
    if isinstance(enc, SortEncoder):
        return sorted(evaluate(a, event) for a in enc._orbit)
    raise TypeError(f"Unknown encoder type: {type(enc)}")


# ---------------------------------------------------------------------------
# Integration test 1: encode_and_describe segment offsets partition the array
# ---------------------------------------------------------------------------

def test_segment_offsets_partition_encoded_array(
    ops, ctx, event, orbit_factory, phase2_factory, atol
):
    """Segment start/stop offsets from the tree exactly partition the encoded array.

    Verifies:
    - total length of encoded equals sum of segment lengths from the tree
    - every segment has stop == start + length
    - non-null segments appear in increasing start order with no gaps
    """
    mag, dot, eps3 = ops
    plan    = Plan(context=ctx, operations=(mag, dot, eps3))
    encoded, tree = encode_and_describe(plan, event, orbit_factory, phase2_factory)

    all_segs = list(tree)

    # Total length consistency.
    assert len(encoded) == sum(s.length for s in all_segs), (
        "encode() length ≠ sum of segment lengths"
    )

    # stop == start + length for every segment.
    for s in all_segs:
        assert s.stop == s.start + s.length, (
            f"Segment {s}: stop={s.stop} ≠ start={s.start} + length={s.length}"
        )

    # Non-null segments are ordered with no gaps between them.
    non_null = [s for s in all_segs if s.length > 0]
    for i in range(len(non_null) - 1):
        assert non_null[i].stop == non_null[i + 1].start, (
            f"Gap between segments {non_null[i]} and {non_null[i+1]}"
        )

    # They must start at 0 and end at len(encoded).
    if non_null:
        assert non_null[0].start == 0
        assert non_null[-1].stop == len(encoded)


# ---------------------------------------------------------------------------
# Integration test 2: encode() == encode_and_describe()[0]
# ---------------------------------------------------------------------------

def test_encode_equals_encode_and_describe_values(
    ops, ctx, event, orbit_factory, phase2_factory
):
    """encode() and encode_and_describe() return identical value arrays.

    encode() and encode_and_describe() share one code path; this test verifies
    that the public wrappers are consistent.
    """
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))

    values_a = encode(plan, event, orbit_factory, phase2_factory)
    values_b, _ = encode_and_describe(plan, event, orbit_factory, phase2_factory)

    np.testing.assert_array_equal(values_a, values_b)


# ---------------------------------------------------------------------------
# Integration test 3: Phase 1 slices match direct atom evaluations
# ---------------------------------------------------------------------------

def test_phase1_slices_match_atom_evaluations(
    ops, ctx, event, orbit_factory, phase2_factory, atol
):
    """Phase 1 slices in the encoded array match direct atom evaluations.

    For each ORBIT segment in tree.phase1, encoded[s.start:s.stop] must equal
    the sorted values produced by evaluating the orbit atoms directly.  This
    cross-checks the tree's start/stop metadata against the actual content
    of the encoded vector, using evaluate() as an independent oracle.
    """
    mag, dot, eps3 = ops
    plan    = Plan(context=ctx, operations=(mag, dot, eps3))
    encoded, tree = encode_and_describe(plan, event, orbit_factory, phase2_factory)

    orbit_enc = orbit_factory.build(plan)

    phase1_segs   = list(tree.phase1)
    enc_selections = orbit_enc._selections

    assert len(phase1_segs) == len(enc_selections), (
        "Number of Phase 1 tree segments ≠ number of orbit encoders"
    )

    for seg, (fo, enc) in zip(phase1_segs, enc_selections):
        assert seg.kind == "ORBIT", f"Expected ORBIT segment, got {seg.kind}"
        assert seg.length == enc.output_dim, (
            f"Segment length {seg.length} ≠ encoder output_dim {enc.output_dim}"
        )

        slice_vals    = list(encoded[seg.start:seg.stop])
        ground_truth  = _phase1_ground_truth(enc, event)

        np.testing.assert_allclose(
            sorted(slice_vals), sorted(ground_truth), atol=atol,
            err_msg=(
                f"Phase 1 slice mismatch for {fo!r}: "
                f"slice={sorted(slice_vals)}, truth={sorted(ground_truth)}"
            ),
        )


# ---------------------------------------------------------------------------
# Integration test 4: describe_encoding() is encode-event-agnostic
# ---------------------------------------------------------------------------

def test_describe_encoding_independent_of_event(
    ops, ctx, event, orbit_factory, phase2_factory
):
    """describe_encoding() produces the same tree regardless of event values.

    The tree (segment kinds, lengths, offsets) is determined purely by the
    plan's algebraic structure; it must not depend on the numerical event.
    Verify this by comparing the tree from describe_encoding() against the
    tree extracted from encode_and_describe() with the actual event.
    """
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))

    tree_from_describe      = describe_encoding(plan, orbit_factory, phase2_factory)
    _, tree_from_encode_desc = encode_and_describe(plan, event, orbit_factory, phase2_factory)

    segs_a = list(tree_from_describe)
    segs_b = list(tree_from_encode_desc)

    assert len(segs_a) == len(segs_b), (
        "describe_encoding() and encode_and_describe() return different numbers of segments"
    )
    for a, b in zip(segs_a, segs_b):
        assert a == b, (
            f"Segment mismatch:\n  describe_encoding: {a}\n  encode_and_describe: {b}"
        )


# ---------------------------------------------------------------------------
# Integration test 5: permutation invariance on the two-group context
# ---------------------------------------------------------------------------

def test_permutation_invariance_two_group_full_plan(
    ops, ctx, orbit_factory, phase2_factory
):
    """encode() is permutation-invariant for the full (mag, dot, eps3) plan.

    Uses the two-group (electrons + muons) context from conftest and swaps
    all electron labels (a↔b↔c cyclic).  The encoded vector must be identical.
    """
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))

    event = {
        "a": np.array([1.0, 0.2, 0.3]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.3, 1.0]),
        "p": np.array([0.8, 0.6, 0.1]),
        "q": np.array([0.2, 0.9, 0.5]),
    }
    # Cyclic permutation a→b→c→a (leaves muons unchanged).
    perm = {"a": "b", "b": "c", "c": "a", "p": "p", "q": "q"}
    event_permuted = {perm[k]: v for k, v in event.items()}

    enc_orig = encode(plan, event,          orbit_factory, phase2_factory)
    enc_perm = encode(plan, event_permuted, orbit_factory, phase2_factory)

    np.testing.assert_array_almost_equal(enc_orig, enc_perm)
