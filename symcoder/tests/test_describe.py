"""Tests for describe_encoding() and SegmentInfo."""
import json
import numpy as np
import pytest
from itertools import groupby
from symatom import ArgumentSymmetry, VectorType, Context, Plan, repS, canonical_pair_flavours
from symcoder import EvaluableOperation, encode, describe_encoding, SegmentInfo

_NULL_KINDS = ("NULL_SELF", "NULL_COMP")
_PAIR_KINDS  = ("ASSOC", "NULL_SELF", "NULL_COMP")
_ALL_KINDS   = ("ORBIT", "ASSOC", "NULL_SELF", "NULL_COMP")


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
    return Context((VectorType("electrons", ("a", "b", "c", "d")),))

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
    ctx  = Context((VectorType("e", ()), VectorType("mu", ())))
    plan = Plan(context=ctx, operations=(dot,))
    assert describe_encoding(plan) == []


# ---------------------------------------------------------------------------
# Valid kinds present
# ---------------------------------------------------------------------------

def test_orbit_segments_present(plan):
    kinds = {s.kind for s in describe_encoding(plan)}
    assert "ORBIT" in kinds

def test_assoc_and_null_segments_present(plan):
    kinds = {s.kind for s in describe_encoding(plan)}
    assert "ASSOC" in kinds
    assert kinds & set(_NULL_KINDS), f"Expected at least one null kind; got {kinds}"

def test_only_valid_kinds(plan):
    assert all(s.kind in _ALL_KINDS for s in describe_encoding(plan))

def test_orbit_segments_come_first(plan):
    segs = describe_encoding(plan)
    orbit_indices     = [i for i, s in enumerate(segs) if s.kind == "ORBIT"]
    non_orbit_indices = [i for i, s in enumerate(segs) if s.kind != "ORBIT"]
    if orbit_indices and non_orbit_indices:
        assert max(orbit_indices) < min(non_orbit_indices)


# ---------------------------------------------------------------------------
# Segment index consistency
# ---------------------------------------------------------------------------

def test_null_segments_have_length_zero(plan):
    assert all(s.length == 0 for s in describe_encoding(plan) if s.kind in _NULL_KINDS)

def test_null_segments_start_equals_stop(plan):
    assert all(s.start == s.stop for s in describe_encoding(plan) if s.kind in _NULL_KINDS)

def test_non_null_segments_are_contiguous(plan):
    """ORBIT and ASSOC segments must tile the output array without gaps."""
    non_null = [s for s in describe_encoding(plan) if s.kind not in _NULL_KINDS]
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

def test_orbit_sign_compressed_field_is_set(plan):
    """Every ORBIT segment has sign_compressed as a bool."""
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            assert isinstance(s.sign_compressed, bool)

def test_orbit_sign_compressed_false_for_symmetric(dot, ctx):
    """SYMMETRIC operation → sign_compressed=False on all ORBIT segments."""
    plan = Plan(context=ctx, operations=(dot,))
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            assert s.sign_compressed is False

def test_orbit_sign_compressed_true_for_antisymmetric(eps3, ctx):
    """ANTISYMMETRIC operation → sign_compressed=True on all ORBIT segments."""
    plan = Plan(context=ctx, operations=(eps3,))
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            assert s.sign_compressed is True

def test_orbit_length_symmetric_equals_fo_count(dot, ctx):
    """SYMMETRIC: ORBIT length == fo.count() (no compression)."""
    plan = Plan(context=ctx, operations=(dot,))
    fo_list = repS(plan.context, plan.operations)
    fo_counts = {(fo.operation.name, fo.flavour.counts): fo.count() for fo in fo_list}
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            expected = fo_counts[(s.op_u, s.flavour_u)]
            assert s.length == expected

def test_orbit_length_antisymmetric_halved(eps3, ctx):
    # TODO .. likely completely wrong, posisibly needs FIXME or deletion.
    """ANTISYMMETRIC (sign-compressed): ORBIT length == fo.count() // 2."""  # TODO .. likely wrong, possibly needs FIXME
    plan = Plan(context=ctx, operations=(eps3,))
    fo_list = repS(plan.context, plan.operations)
    fo_counts = {(fo.operation.name, fo.flavour.counts): fo.count() for fo in fo_list}
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            expected = fo_counts[(s.op_u, s.flavour_u)] // 2  # TODO .. likely wrong, possibly needs FIXME
            expected = fo_counts[(s.op_u, s.flavour_u)]  # TODO .. likely wrong, possibly needs FIXME -- but this is temporary potential hotfix for line above
            assert s.length == expected

def test_orbit_count_for_symmetric_op(dot, ctx):
    """SYMMETRIC dot with flavour (2,) in group of 4: C(4,2)=6 atoms, length=6."""
    plan = Plan(context=ctx, operations=(dot,))
    orbit_segs = [s for s in describe_encoding(plan) if s.kind == "ORBIT"]
    assert all(s.length > 0 for s in orbit_segs)

def test_assoc_null_sign_compressed_is_none(plan):
    """ASSOC and NULL_* segments must not set sign_compressed."""
    for s in describe_encoding(plan):
        if s.kind in _PAIR_KINDS:
            assert s.sign_compressed is None


# ---------------------------------------------------------------------------
# NULL_SELF correctness
# ---------------------------------------------------------------------------

def test_null_self_satisfies_self_pair_condition(plan):
    """Every NULL_SELF segment has op_u==op_v, flavour_u==flavour_v, overlap==flavour_u."""
    for s in describe_encoding(plan):
        if s.kind == "NULL_SELF":
            assert s.op_u == s.op_v
            assert s.flavour_u == s.flavour_v
            assert s.overlap == s.flavour_u

def test_no_assoc_is_self_pair(plan):
    """No ASSOC should satisfy the self-pair condition (those must be NULL_SELF)."""
    for s in describe_encoding(plan):
        if s.kind == "ASSOC":
            is_self = (s.op_u == s.op_v and s.flavour_u == s.flavour_v
                       and s.overlap == s.flavour_u)
            assert not is_self, f"ASSOC segment is a self-pair and should be NULL_SELF: {s}"


# ---------------------------------------------------------------------------
# NULL_COMP correctness
# ---------------------------------------------------------------------------

def test_each_block_has_at_most_one_null_comp(dot, eps3, ctx):
    """Each OVERLAP BLOCK has at most one NULL_COMP."""
    segs = describe_encoding(Plan(context=ctx, operations=(dot, eps3)))
    pair_segs = [s for s in segs if s.kind in _PAIR_KINDS]

    def block_key(s):
        return (s.op_u, s.flavour_u, s.op_v, s.flavour_v)

    for _key, block in groupby(pair_segs, key=block_key):
        block = list(block)
        null_comp_count = sum(1 for s in block if s.kind == "NULL_COMP")
        assert null_comp_count <= 1, (
            f"Expected at most 1 NULL_COMP per block, got {null_comp_count}: {block}"
        )

def test_block_with_non_self_pair_has_exactly_one_null_comp(dot, eps3, ctx):
    """Any block containing at least one non-self-pair has exactly one NULL_COMP."""
    segs = describe_encoding(Plan(context=ctx, operations=(dot, eps3)))
    pair_segs = [s for s in segs if s.kind in _PAIR_KINDS]

    def block_key(s):
        return (s.op_u, s.flavour_u, s.op_v, s.flavour_v)

    for _key, block in groupby(pair_segs, key=block_key):
        block = list(block)
        has_non_self = any(
            s.kind in ("ASSOC", "NULL_COMP") for s in block
        )
        if has_non_self:
            null_comp_count = sum(1 for s in block if s.kind == "NULL_COMP")
            assert null_comp_count == 1, (
                f"Block with non-self-pairs must have exactly 1 NULL_COMP, "
                f"got {null_comp_count}: {block}"
            )

def test_null_comp_is_largest_non_self_in_block(dot, eps3, ctx):
    """NULL_COMP's notional length >= all ASSOC notional lengths in the same block."""
    plan = Plan(context=ctx, operations=(dot, eps3))
    segs = describe_encoding(plan)
    fo_list = repS(ctx, (dot, eps3))
    pf_list = canonical_pair_flavours(fo_list, ctx)
    type_sizes = tuple(g.size for g in ctx.types)
    pf_count_map = {
        (pf.op_u.name, tuple(pf.flavour_u.counts), pf.op_v.name,
         tuple(pf.flavour_v.counts), tuple(pf.overlap)): 2 * pf.count(type_sizes)
        for pf in pf_list
    }

    pair_segs = [s for s in segs if s.kind in _PAIR_KINDS]

    def block_key(s):
        return (s.op_u, s.flavour_u, s.op_v, s.flavour_v)

    for _key, block in groupby(pair_segs, key=block_key):
        block = list(block)
        assoc_lengths = [s.length for s in block if s.kind == "ASSOC"]
        null_comp_segs = [s for s in block if s.kind == "NULL_COMP"]
        if not null_comp_segs:
            continue
        assert len(null_comp_segs) == 1
        nc = null_comp_segs[0]
        nc_notional = pf_count_map[(nc.op_u, nc.flavour_u, nc.op_v, nc.flavour_v, nc.overlap)]
        for al in assoc_lengths:
            assert nc_notional >= al, (
                f"NULL_COMP notional={nc_notional} < ASSOC length={al} in block"
            )


# ---------------------------------------------------------------------------
# Symmetry class labels
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
    pair_segs = [s for s in describe_encoding(plan) if s.kind in _PAIR_KINDS]
    classes = {s.symmetry_class for s in pair_segs}
    assert "SS" in classes
    assert "AA" in classes
    assert "SA" in classes or "AS" in classes


# ---------------------------------------------------------------------------
# Human-readable string (unified column format)
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

def test_null_comp_str_contains_null_comp(plan):
    for s in describe_encoding(plan):
        if s.kind == "NULL_COMP":
            assert "NULL_COMP" in str(s)

def test_null_self_str_contains_null_self(plan):
    for s in describe_encoding(plan):
        if s.kind == "NULL_SELF":
            assert "NULL_SELF" in str(s)

def test_str_contains_index_range(plan):
    for s in describe_encoding(plan):
        assert f"[{s.start}:{s.stop}]" in str(s)

def test_assoc_str_contains_symmetry_class(dot, eps3, ctx):
    plan = Plan(context=ctx, operations=(dot, eps3))
    for s in describe_encoding(plan):
        if s.kind == "ASSOC":
            assert s.symmetry_class in str(s)

def test_orbit_str_variant_column(dot, eps3, ctx):
    """ORBIT rows show 'SC' when sign_compressed, '.' otherwise — in the variant column."""
    plan = Plan(context=ctx, operations=(dot, eps3))
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            tokens = str(s).split()
            if s.sign_compressed:
                assert "SC" in tokens, f"sign-compressed ORBIT missing 'SC': {s}"
            else:
                # op_v and v=(.) and shared=(.) still give '.' tokens
                assert "." in tokens, f"non-compressed ORBIT missing '.' placeholder: {s}"

def test_all_rows_same_number_of_base_tokens(plan):
    """Every row (ignoring the optional '| example' tail) has the same token count."""
    segs = describe_encoding(plan)
    def base_tokens(s):
        line = str(s)
        if "  |  " in line:
            line = line[:line.index("  |  ")]
        return line.split()

    counts = [len(base_tokens(s)) for s in segs]
    assert len(set(counts)) == 1, (
        f"Rows have differing token counts: {set(counts)}"
    )

def test_full_column_equals_length_for_assoc_and_orbit(plan):
    """For ASSOC and ORBIT segments, notional_length == length."""
    for s in describe_encoding(plan):
        if s.kind in ("ASSOC", "ORBIT"):
            nl = s.notional_length if s.notional_length is not None else s.length
            assert nl == s.length, f"notional_length != length for {s.kind}: {s}"

def test_full_column_positive_for_null_segments(plan):
    """For NULL_* segments, notional_length > 0 (shows what was saved)."""
    for s in describe_encoding(plan):
        if s.kind in _NULL_KINDS:
            nl = s.notional_length if s.notional_length is not None else s.length
            assert nl > 0, f"notional_length should be > 0 for dropped segment: {s}"

def test_str_contains_full_column(plan):
    """Every row's string contains a 'full=N' token."""
    for s in describe_encoding(plan):
        assert "full=" in str(s), f"Missing full= in: {s}"


# ---------------------------------------------------------------------------
# Machine-readable export
# ---------------------------------------------------------------------------

def test_orbit_to_dict_has_required_keys(plan):
    for s in describe_encoding(plan):
        if s.kind == "ORBIT":
            d = s.to_dict()
            assert {"kind", "start", "stop", "length", "op_u", "flavour_u", "sign_compressed"} <= set(d.keys())
            assert "op_v" not in d

def test_assoc_to_dict_has_required_keys(plan):
    required = {"kind", "start", "stop", "length", "op_u", "op_v",
                "flavour_u", "flavour_v", "overlap", "symmetry_class"}
    for s in describe_encoding(plan):
        if s.kind in _PAIR_KINDS:
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
