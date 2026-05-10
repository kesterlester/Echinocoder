"""Tests for describe_encoding() and SegmentInfo."""
import json
import numpy as np
import pytest
from symatom import ArgumentSymmetry, VectorGroup, Context, Plan, repL, canonical_pair_flavours
from symcoder import EvaluableOperation, encode, describe_encoding, SegmentInfo


@pytest.fixture
def dot():
    return EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )

@pytest.fixture
def eps3():
    return EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )

@pytest.fixture
def ctx():
    return Context((VectorGroup("electrons", ("a", "b", "c", "d")),))

@pytest.fixture
def plan(dot, eps3, ctx):
    return Plan(context=ctx, operations=(dot, eps3))


# ---------------------------------------------------------------------------
# Return type and basic structure
# ---------------------------------------------------------------------------

def test_returns_list_of_segment_info(plan):
    segs = describe_encoding(plan)
    assert isinstance(segs, list)
    assert all(isinstance(s, SegmentInfo) for s in segs)

def test_non_empty_for_non_trivial_plan(plan):
    assert len(describe_encoding(plan)) > 0

def test_empty_ops_gives_empty_list(ctx):
    plan = Plan(context=ctx, operations=())
    assert describe_encoding(plan) == []

def test_all_empty_groups_gives_empty_list(dot):
    ctx  = Context((VectorGroup("e", ()), VectorGroup("mu", ())))
    plan = Plan(context=ctx, operations=(dot,))
    assert describe_encoding(plan) == []


# ---------------------------------------------------------------------------
# Segment index consistency
# ---------------------------------------------------------------------------

def test_segments_are_contiguous(plan):
    segs = describe_encoding(plan)
    assert segs[0].start == 0
    for a, b in zip(segs, segs[1:]):
        assert a.stop == b.start

def test_total_length_matches_encode_output(plan):
    event = {l: np.random.randn(3) for l in ("a", "b", "c", "d")}
    total = sum(s.length for s in describe_encoding(plan))
    assert total == len(encode(plan, event))

def test_stop_property(plan):
    for s in describe_encoding(plan):
        assert s.stop == s.start + s.length


# ---------------------------------------------------------------------------
# Symmetry class labels
# ---------------------------------------------------------------------------

def test_symmetry_class_ss_for_dot_dot(dot, ctx):
    plan = Plan(context=ctx, operations=(dot,))
    segs = describe_encoding(plan)
    assert all(s.symmetry_class == "SS" for s in segs)

def test_symmetry_class_aa_for_eps3_eps3(eps3, ctx):
    plan = Plan(context=ctx, operations=(eps3,))
    segs = describe_encoding(plan)
    assert all(s.symmetry_class == "AA" for s in segs)

def test_symmetry_class_sa_for_mixed(dot, eps3, ctx):
    plan = Plan(context=ctx, operations=(dot, eps3))
    segs = describe_encoding(plan)
    classes = {s.symmetry_class for s in segs}
    assert "SS" in classes
    assert "AA" in classes
    assert "SA" in classes or "AS" in classes


# ---------------------------------------------------------------------------
# Segment lengths match pf.count()
# ---------------------------------------------------------------------------

def test_segment_lengths_match_pf_count(dot, eps3, ctx):
    plan    = Plan(context=ctx, operations=(dot, eps3))
    segs    = describe_encoding(plan)
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    pf_counts = [pf.count(group_sizes) for pf in canonical_pair_flavours(fo_list, ctx)]
    assert [s.length for s in segs] == pf_counts


# ---------------------------------------------------------------------------
# Human-readable string
# ---------------------------------------------------------------------------

def test_str_contains_op_names(dot, ctx):
    plan = Plan(context=ctx, operations=(dot,))
    for s in describe_encoding(plan):
        assert "dot" in str(s)

def test_str_contains_index_range(plan):
    for s in describe_encoding(plan):
        assert f"[{s.start}:{s.stop}]" in str(s)

def test_str_contains_symmetry_class(plan):
    for s in describe_encoding(plan):
        assert s.symmetry_class in str(s)


# ---------------------------------------------------------------------------
# Machine-readable export
# ---------------------------------------------------------------------------

def test_to_dict_has_required_keys(plan):
    required = {"start", "stop", "length", "op_u", "op_v",
                "flavour_u", "flavour_v", "overlap", "symmetry_class"}
    for s in describe_encoding(plan):
        assert required <= set(s.to_dict().keys())

def test_to_dict_is_json_serialisable(plan):
    dicts = [s.to_dict() for s in describe_encoding(plan)]
    dumped = json.dumps(dicts)
    reloaded = json.loads(dumped)
    assert len(reloaded) == len(dicts)
    assert reloaded[0]["start"] == 0

def test_to_dict_start_stop_consistent(plan):
    for s in describe_encoding(plan):
        d = s.to_dict()
        assert d["stop"] == d["start"] + d["length"]


# ---------------------------------------------------------------------------
# Pure function — no event data needed
# ---------------------------------------------------------------------------

def test_describe_encoding_is_pure(plan):
    """describe_encoding must return the same result every call, no event needed."""
    assert describe_encoding(plan) == describe_encoding(plan)
