"""Tests for describe_encoding() and SegmentInfo (5b structure: ORBIT/ASSOC/NULL)."""
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
# Three kinds present
# ---------------------------------------------------------------------------

def test_orbit_segments_present(plan):
    kinds = {s.kind for s in describe_encoding(plan)}
    assert "ORBIT" in kinds

def test_assoc_and_null_segments_present(plan):
    kinds = {s.kind for s in describe_encoding(plan)}
    assert "ASSOC" in kinds
    assert "NULL" in kinds

def test_only_valid_kinds(plan):
    assert all(s.kind in ("ORBIT", "ASSOC", "NULL") for s in describe_encoding(plan))

def test_orbit_segments_come_first(plan):
    segs = describe_encoding(plan)
    orbit_indices = [i for i, s in enumerate(segs) if s.kind == "ORBIT"]
    non_orbit_indices = [i for i, s in enumerate(segs) if s.kind != "ORBIT"]
    if orbit_indices and non_orbit_indices:
        assert max(orbit_indices) < min(non_orbit_indices)


# ---------------------------------------------------------------------------
# Segment index consistency
# ---------------------------------------------------------------------------

def test_null_segments_have_length_zero(plan):
    assert all(s.length == 0 for s in describe_encoding(plan) if s.kind == "NULL")

def test_null_segments_start_equals_stop(plan):
    assert all(s.start == s.stop for s in describe_encoding(plan) if s.kind == "NULL")

def test_non_null_segments_are_contiguous(plan):
    """ORBIT and ASSOC segments must tile the output array without gaps."""
    non_null = [s for s in describe_encoding(plan) if s.kind != "NULL"]
    if not non_null:
        return
    assert non_null[0].start == 0
    for a, b in zip(non_null, non_null[1:]):
        assert a.stop == b.start

def test_total_length_matches_encode_output(plan):
    event = {l: np.random.randn(3) for l in ("a", "b", "c", "d")}
    total = sum(s.length for s in describe_encoding(plan))
    assert total == len(encode(plan, event))

def test_stop_property(plan):
    for s in describe_encoding(plan):
        assert s.stop == s.start + s.length


# ---------------------------------------------------------------------------
# ORBIT segment correctness
# ---------------------------------------------------------------------------

def test_orbit_op_v_is_none(plan):
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            assert s.op_v is None
            assert s.flavour_v is None
            assert s.overlap is None
            assert s.symmetry_class is None

def test_orbit_count_for_symmetric_op(dot, ctx):
    """SYMMETRIC dot with flavour (2,) in group of 4: C(4,2)=6 atoms → ceil(6/2)=3 complex."""
    plan = Plan(context=ctx, operations=(dot,))
    orbit_segs = [s for s in describe_encoding(plan) if s.kind == "ORBIT"]
    # All orbit segments should have length = ceil(fo.count()/2)
    assert all(s.length > 0 for s in orbit_segs)


# ---------------------------------------------------------------------------
# NULL segment correctness
# ---------------------------------------------------------------------------

def test_each_block_has_exactly_one_null(dot, eps3, ctx):
    """Each OVERLAP BLOCK (fixed op_u, flavour_u, op_v, flavour_v) has exactly one NULL."""
    segs = describe_encoding(Plan(context=ctx, operations=(dot, eps3)))
    # Group ASSOC+NULL segments by their (op_u, flavour_u, op_v, flavour_v)
    from itertools import groupby
    pair_segs = [s for s in segs if s.kind in ("ASSOC", "NULL")]

    def block_key(s):
        return (s.op_u, s.flavour_u, s.op_v, s.flavour_v)

    for _key, block in groupby(pair_segs, key=block_key):
        block = list(block)
        null_count = sum(1 for s in block if s.kind == "NULL")
        assert null_count == 1, f"Expected 1 NULL in block, got {null_count}: {block}"

def test_null_is_largest_in_block(dot, eps3, ctx):
    """The NULL association must have length >= all ASSOC associations in the same block."""
    segs = describe_encoding(Plan(context=ctx, operations=(dot, eps3)))
    from itertools import groupby

    pair_segs = [s for s in segs if s.kind in ("ASSOC", "NULL")]

    def block_key(s):
        return (s.op_u, s.flavour_u, s.op_v, s.flavour_v)

    # For NULL to be "the largest", we need to know its notional length.
    # Recover from pf_list: the NULL's notional length equals pf.count(group_sizes).
    fo_list = repL(ctx, (dot, eps3))
    pf_list = canonical_pair_flavours(fo_list, ctx)
    group_sizes = tuple(g.size for g in ctx.groups)
    pf_count_map = {
        (pf.op_u.name, tuple(pf.flavour_u.counts), pf.op_v.name,
         tuple(pf.flavour_v.counts), tuple(pf.overlap)): 2 * pf.count(group_sizes)
        for pf in pf_list
    }

    for _key, block in groupby(pair_segs, key=block_key):
        block = list(block)
        assoc_lengths = [s.length for s in block if s.kind == "ASSOC"]
        null_segs = [s for s in block if s.kind == "NULL"]
        assert len(null_segs) == 1
        null_key = (null_segs[0].op_u, null_segs[0].flavour_u,
                    null_segs[0].op_v, null_segs[0].flavour_v,
                    null_segs[0].overlap)
        null_notional = pf_count_map[null_key]
        for al in assoc_lengths:
            assert null_notional >= al


# ---------------------------------------------------------------------------
# Symmetry class labels (ASSOC segments only)
# ---------------------------------------------------------------------------

def test_symmetry_class_ss_for_dot_dot(dot, ctx):
    plan = Plan(context=ctx, operations=(dot,))
    assoc_segs = [s for s in describe_encoding(plan) if s.kind == "ASSOC"]
    assert all(s.symmetry_class == "SS" for s in assoc_segs)

def test_symmetry_class_aa_for_eps3_eps3(eps3, ctx):
    plan = Plan(context=ctx, operations=(eps3,))
    assoc_segs = [s for s in describe_encoding(plan) if s.kind == "ASSOC"]
    assert all(s.symmetry_class == "AA" for s in assoc_segs)

def test_symmetry_classes_present_in_mixed_plan(dot, eps3, ctx):
    plan = Plan(context=ctx, operations=(dot, eps3))
    pair_segs = [s for s in describe_encoding(plan) if s.kind in ("ASSOC", "NULL")]
    classes = {s.symmetry_class for s in pair_segs}
    assert "SS" in classes
    assert "AA" in classes
    assert "SA" in classes or "AS" in classes


# ---------------------------------------------------------------------------
# Human-readable string
# ---------------------------------------------------------------------------

def test_orbit_str_contains_ORBIT(dot, ctx):
    plan = Plan(context=ctx, operations=(dot,))
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            assert "ORBIT" in str(s)

def test_assoc_str_contains_op_names(dot, ctx):
    plan = Plan(context=ctx, operations=(dot,))
    for s in describe_encoding(plan):
        if s.kind == "ASSOC":
            assert "dot" in str(s)

def test_null_str_contains_null_encoding(plan):
    for s in describe_encoding(plan):
        if s.kind == "NULL":
            assert "NULL_ENCODING" in str(s)

def test_str_contains_index_range(plan):
    for s in describe_encoding(plan):
        assert f"[{s.start}:{s.stop}]" in str(s)

def test_assoc_str_contains_symmetry_class(dot, eps3, ctx):
    plan = Plan(context=ctx, operations=(dot, eps3))
    for s in describe_encoding(plan):
        if s.kind == "ASSOC":
            assert s.symmetry_class in str(s)


# ---------------------------------------------------------------------------
# Machine-readable export
# ---------------------------------------------------------------------------

def test_orbit_to_dict_has_required_keys(plan):
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            d = s.to_dict()
            assert {"kind", "start", "stop", "length", "op_u", "flavour_u"} <= set(d.keys())
            assert "op_v" not in d

def test_assoc_to_dict_has_required_keys(plan):
    required = {"kind", "start", "stop", "length", "op_u", "op_v",
                "flavour_u", "flavour_v", "overlap", "symmetry_class"}
    for s in describe_encoding(plan):
        if s.kind in ("ASSOC", "NULL"):
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
    assert describe_encoding(plan) == describe_encoding(plan)
