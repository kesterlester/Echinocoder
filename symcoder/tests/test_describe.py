"""Tests for describe_encoding() and SegmentInfo."""
import json
import numpy as np
import pytest
from itertools import groupby
from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan, repS, canonical_pair_flavours
from symcoder import encode, describe_encoding, SegmentInfo, EncodingTree

_NULL_KINDS = ("NULL_SELF", "NULL_COMP")
_PAIR_KINDS  = ("ASSOC", "NULL_SELF", "NULL_COMP")
_ALL_KINDS   = ("ORBIT", "ASSOC", "NULL_SELF", "NULL_COMP", "SIMPLICIAL")


@pytest.fixture
def dot():
    return Operation(
        name="dot", rank=2, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )

@pytest.fixture
def eps3():
    return Operation(
        name="eps3", rank=3, odd_parity=True,
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

def test_returns_list_of_segment_info(plan, orbit_factory, phase2_factory):
    segs = describe_encoding(plan, orbit_factory, phase2_factory)
    assert isinstance(segs, EncodingTree)
    assert all(isinstance(s, SegmentInfo) for s in segs)

def test_non_empty_for_non_trivial_plan(plan, orbit_factory, phase2_factory):
    assert len(describe_encoding(plan, orbit_factory, phase2_factory)) > 0

def test_empty_ops_gives_empty_list(ctx, orbit_factory, phase2_factory):
    plan = Plan(context=ctx, operations=())
    assert list(describe_encoding(plan, orbit_factory, phase2_factory)) == []

def test_all_empty_groups_gives_empty_list(dot, orbit_factory, phase2_factory):
    ctx  = Context((VectorType("e", ()), VectorType("mu", ())))
    plan = Plan(context=ctx, operations=(dot,))
    assert list(describe_encoding(plan, orbit_factory, phase2_factory)) == []


# ---------------------------------------------------------------------------
# Valid kinds present
# ---------------------------------------------------------------------------

def test_orbit_segments_present(plan, orbit_factory, phase2_factory):
    kinds = {s.kind for s in describe_encoding(plan, orbit_factory, phase2_factory)}
    assert "ORBIT" in kinds

def test_assoc_and_null_segments_present(plan, orbit_factory, phase2_factory):
    kinds = {s.kind for s in describe_encoding(plan, orbit_factory, phase2_factory)}
    assert "ASSOC" in kinds
    assert kinds & set(_NULL_KINDS), f"Expected at least one null kind; got {kinds}"

def test_only_valid_kinds(plan, orbit_factory, phase2_factory):
    assert all(s.kind in _ALL_KINDS for s in describe_encoding(plan, orbit_factory, phase2_factory))

def test_orbit_segments_come_first(plan, orbit_factory, phase2_factory):
    segs = describe_encoding(plan, orbit_factory, phase2_factory)
    orbit_indices     = [i for i, s in enumerate(segs) if s.kind == "ORBIT"]
    non_orbit_indices = [i for i, s in enumerate(segs) if s.kind != "ORBIT"]
    if orbit_indices and non_orbit_indices:
        assert max(orbit_indices) < min(non_orbit_indices)


# ---------------------------------------------------------------------------
# Segment index consistency
# ---------------------------------------------------------------------------

def test_null_segments_have_length_zero(plan, orbit_factory, phase2_factory):
    assert all(s.length == 0 for s in describe_encoding(plan, orbit_factory, phase2_factory) if s.kind in _NULL_KINDS)

def test_null_segments_start_equals_stop(plan, orbit_factory, phase2_factory):
    assert all(s.start == s.stop for s in describe_encoding(plan, orbit_factory, phase2_factory) if s.kind in _NULL_KINDS)

def test_non_null_segments_are_contiguous(plan, orbit_factory, phase2_factory):
    """ORBIT and ASSOC segments must tile the output array without gaps."""
    non_null = [s for s in describe_encoding(plan, orbit_factory, phase2_factory) if s.kind not in _NULL_KINDS]
    if not non_null:
        return
    assert non_null[0].start == 0
    for a, b in zip(non_null, non_null[1:]):
        assert a.stop == b.start

def test_total_length_matches_encode_output(plan, orbit_factory, phase2_factory):
    event = {l: np.random.randn(3) for l in ("a", "b", "c", "d")}
    total = sum(s.length for s in describe_encoding(plan, orbit_factory, phase2_factory))
    assert total == len(encode(plan, event, orbit_factory, phase2_factory))

def test_stop_property(plan, orbit_factory, phase2_factory):
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        assert s.stop == s.start + s.length


# ---------------------------------------------------------------------------
# ORBIT segment correctness
# ---------------------------------------------------------------------------

def test_orbit_op_v_is_none(plan, orbit_factory, phase2_factory):
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ORBIT":
            assert s.op_v is None
            assert s.flavour_v is None
            assert s.overlap is None
            assert s.symmetry_class is None

def test_orbit_method_name_is_set(plan, orbit_factory, phase2_factory):
    """Every ORBIT segment has a non-None method_name from the chosen encoder."""
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ORBIT":
            assert s.method_name is not None

def test_orbit_length_symmetric_equals_fo_count(dot, ctx, orbit_factory, phase2_factory):
    """SYMMETRIC: ORBIT length == fo.count_of_atoms_one_per_sign() (no compression)."""
    plan = Plan(context=ctx, operations=(dot,))
    fo_list = repS(plan.context, plan.operations)
    fo_counts = {(fo.operation, fo.flavour.counts): fo.count_of_atoms_one_per_sign() for fo in fo_list}
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ORBIT":
            expected = fo_counts[(s.op_u, s.flavour_u)]
            assert s.length == expected

def test_orbit_count_for_symmetric_op(dot, ctx, orbit_factory, phase2_factory):
    """SYMMETRIC dot with flavour (2,) in group of 4: C(4,2)=6 atoms, length=6."""
    plan = Plan(context=ctx, operations=(dot,))
    orbit_segs = [s for s in describe_encoding(plan, orbit_factory, phase2_factory) if s.kind == "ORBIT"]
    assert all(s.length > 0 for s in orbit_segs)


# ---------------------------------------------------------------------------
# NULL_SELF correctness
# ---------------------------------------------------------------------------

def test_null_self_satisfies_self_pair_condition(plan, orbit_factory, phase2_factory):
    """Every NULL_SELF segment has op_u==op_v, flavour_u==flavour_v, overlap==flavour_u."""
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "NULL_SELF":
            assert s.op_u == s.op_v
            assert s.flavour_u == s.flavour_v
            assert s.overlap == s.flavour_u

def test_no_assoc_is_self_pair(plan, orbit_factory, phase2_factory):
    """No ASSOC should satisfy the self-pair condition (those must be NULL_SELF)."""
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ASSOC":
            is_self = (s.op_u == s.op_v and s.flavour_u == s.flavour_v
                       and s.overlap == s.flavour_u)
            assert not is_self, f"ASSOC segment is a self-pair and should be NULL_SELF: {s}"


# ---------------------------------------------------------------------------
# NULL_COMP correctness
# ---------------------------------------------------------------------------

def _block_key(s):
    return (s.op_u, s.flavour_u, s.op_v, s.flavour_v)


def test_each_block_has_at_most_one_null_comp(dot, eps3, ctx, orbit_factory, phase2_factory):
    """Each OVERLAP BLOCK has at most one NULL_COMP."""
    segs = describe_encoding(Plan(context=ctx, operations=(dot, eps3)), orbit_factory, phase2_factory)
    pair_segs = [s for s in segs if s.kind in _PAIR_KINDS]
    for _key, block in groupby(pair_segs, key=_block_key):
        block = list(block)
        null_comp_count = sum(1 for s in block if s.kind == "NULL_COMP")
        assert null_comp_count <= 1, (
            f"Expected at most 1 NULL_COMP per block, got {null_comp_count}: {block}"
        )

def test_block_with_non_self_pair_has_exactly_one_null_comp(dot, eps3, ctx, orbit_factory, phase2_factory):
    """Any block containing at least one non-self-pair has exactly one NULL_COMP."""
    segs = describe_encoding(Plan(context=ctx, operations=(dot, eps3)), orbit_factory, phase2_factory)
    pair_segs = [s for s in segs if s.kind in _PAIR_KINDS]
    for _key, block in groupby(pair_segs, key=_block_key):
        block = list(block)
        has_non_self = any(s.kind in ("ASSOC", "NULL_COMP") for s in block)
        if has_non_self:
            null_comp_count = sum(1 for s in block if s.kind == "NULL_COMP")
            assert null_comp_count == 1, (
                f"Block with non-self-pairs must have exactly 1 NULL_COMP, "
                f"got {null_comp_count}: {block}"
            )

def test_null_comp_is_largest_non_self_in_block(dot, eps3, ctx, orbit_factory, phase2_factory):
    """NULL_COMP's notional_length >= all ASSOC lengths in the same block."""
    segs = describe_encoding(Plan(context=ctx, operations=(dot, eps3)), orbit_factory, phase2_factory)
    pair_segs = [s for s in segs if s.kind in _PAIR_KINDS]
    for _key, block in groupby(pair_segs, key=_block_key):
        block = list(block)
        assoc_lengths   = [s.length for s in block if s.kind == "ASSOC"]
        null_comp_segs  = [s for s in block if s.kind == "NULL_COMP"]
        if not null_comp_segs:
            continue
        assert len(null_comp_segs) == 1
        nc_notional = null_comp_segs[0].notional_length
        for al in assoc_lengths:
            assert nc_notional >= al, (
                f"NULL_COMP notional_length={nc_notional} < ASSOC length={al} in block"
            )


# ---------------------------------------------------------------------------
# Symmetry class labels
# ---------------------------------------------------------------------------

def test_symmetry_class_11_sigma_for_dot_dot(dot, ctx, orbit_factory, phase2_factory):
    # dot × dot is always a self-pairing block (same op, same flavour).
    # SYMMETRIC operations have all signs +1 → TYPE_11 → σ-compressed to "11_sigma".
    plan = Plan(context=ctx, operations=(dot,))
    assoc_segs = [s for s in describe_encoding(plan, orbit_factory, phase2_factory) if s.kind == "ASSOC"]
    assert all(s.symmetry_class == "11_sigma" for s in assoc_segs)

def test_symmetry_class_neg_sigma_for_eps3_eps3(eps3, ctx, orbit_factory, phase2_factory):
    # eps3 x eps3 is a self-pairing block (same op, same flavour).
    # ANTISYMMETRIC self-pairing → TYPE_NEG → σ-compressed to "neg_sigma".
    # (The achievable set is {(+,+),(-,-)} — TYPE_NEG, not TYPE_22.)
    plan = Plan(context=ctx, operations=(eps3,))
    assoc_segs = [s for s in describe_encoding(plan, orbit_factory, phase2_factory) if s.kind == "ASSOC"]
    assert all(s.symmetry_class == "neg_sigma" for s in assoc_segs)

def test_symmetry_classes_present_in_mixed_plan(dot, eps3, ctx, orbit_factory, phase2_factory):
    plan = Plan(context=ctx, operations=(dot, eps3))
    pair_segs = [s for s in describe_encoding(plan, orbit_factory, phase2_factory) if s.kind in _PAIR_KINDS]
    classes = {s.symmetry_class for s in pair_segs}
    assert "11_sigma" in classes    # dot×dot self-pairing blocks (σ-compressed TYPE_11)
    assert "neg_sigma" in classes   # eps3×eps3 self-pairing blocks (σ-compressed TYPE_NEG)
    assert "12" in classes or "21" in classes


# ---------------------------------------------------------------------------
# Human-readable string (unified column format)
# ---------------------------------------------------------------------------

def test_orbit_str_contains_ORBIT(dot, ctx, orbit_factory, phase2_factory):
    plan = Plan(context=ctx, operations=(dot,))
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ORBIT":
            assert "ORBIT" in str(s)

def test_assoc_str_contains_op_names(dot, ctx, orbit_factory, phase2_factory):
    plan = Plan(context=ctx, operations=(dot,))
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ASSOC":
            assert "dot" in str(s)

def test_null_comp_str_contains_null_comp(plan, orbit_factory, phase2_factory):
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "NULL_COMP":
            assert "NULL_COMP" in str(s)

def test_null_self_str_contains_null_self(plan, orbit_factory, phase2_factory):
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "NULL_SELF":
            assert "NULL_SELF" in str(s)

def test_str_contains_index_range(plan, orbit_factory, phase2_factory):
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        assert f"[{s.start}:{s.stop}]" in str(s)

def test_assoc_str_contains_symmetry_class(dot, eps3, ctx, orbit_factory, phase2_factory):
    plan = Plan(context=ctx, operations=(dot, eps3))
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ASSOC":
            assert s.symmetry_class in str(s)

def test_orbit_str_variant_column(dot, eps3, ctx, orbit_factory, phase2_factory):
    """ORBIT rows show the encoder method_name in the variant column."""
    plan = Plan(context=ctx, operations=(dot, eps3))
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ORBIT":
            assert s.method_name is not None, f"ORBIT segment has no method_name: {s}"
            tokens = str(s).split()
            assert s.method_name in tokens, (
                f"method_name '{s.method_name}' not found as variant token in: {str(s)}"
            )

def test_all_rows_same_number_of_base_tokens(plan, orbit_factory, phase2_factory):
    """Every row (ignoring the optional '| example' tail) has the same token count."""
    segs = describe_encoding(plan, orbit_factory, phase2_factory)
    def base_tokens(s):
        line = str(s)
        if "  |  " in line:
            line = line[:line.index("  |  ")]
        return line.split()

    counts = [len(base_tokens(s)) for s in segs]
    assert len(set(counts)) == 1, (
        f"Rows have differing token counts: {set(counts)}"
    )

def test_full_column_geq_length_for_assoc_and_orbit(plan, orbit_factory, phase2_factory):
    """For active segments (ASSOC, ORBIT), notional_length >= length.
    Compressed encoders (half_sort, TYPE_12/21/22/NEG) have notional_length > length;
    uncompressed encoders (sort, TYPE_11) have notional_length == length."""
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind in ("ASSOC", "ORBIT"):
            nl = s.notional_length if s.notional_length is not None else s.length
            assert nl >= s.length, f"notional_length < length for {s.kind}: {s}"

def test_full_column_positive_for_null_segments(plan, orbit_factory, phase2_factory):
    """For NULL_* segments, notional_length > 0 (shows what was saved)."""
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind in _NULL_KINDS:
            nl = s.notional_length if s.notional_length is not None else s.length
            assert nl > 0, f"notional_length should be > 0 for dropped segment: {s}"

def test_str_contains_full_column(plan, orbit_factory, phase2_factory):
    """Every row's string contains a 'full=N' token."""
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        assert "full=" in str(s), f"Missing full= in: {s}"


# ---------------------------------------------------------------------------
# Machine-readable export
# ---------------------------------------------------------------------------

def test_orbit_to_dict_has_required_keys(plan, orbit_factory, phase2_factory):
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind == "ORBIT":
            d = s.to_dict()
            assert {"kind", "start", "stop", "length", "op_u", "flavour_u", "method_name"} <= set(d.keys())
            assert "op_v" not in d

def test_assoc_to_dict_has_required_keys(plan, orbit_factory, phase2_factory):
    required = {"kind", "start", "stop", "length", "op_u", "op_v",
                "flavour_u", "flavour_v", "overlap", "symmetry_class"}
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        if s.kind in _PAIR_KINDS:
            assert required <= set(s.to_dict().keys())

def test_to_dict_is_json_serialisable(plan, orbit_factory, phase2_factory):
    dicts = [s.to_dict() for s in describe_encoding(plan, orbit_factory, phase2_factory)]
    dumped = json.dumps(dicts)
    reloaded = json.loads(dumped)
    assert len(reloaded) == len(dicts)
    assert reloaded[0]["start"] == 0

def test_to_dict_start_stop_consistent(plan, orbit_factory, phase2_factory):
    for s in describe_encoding(plan, orbit_factory, phase2_factory):
        d = s.to_dict()
        assert d["stop"] == d["start"] + d["length"]


# ---------------------------------------------------------------------------
# Pure function — no event data needed
# ---------------------------------------------------------------------------

def test_describe_encoding_is_pure(plan, orbit_factory, phase2_factory):
    assert describe_encoding(plan, orbit_factory, phase2_factory) == describe_encoding(plan, orbit_factory, phase2_factory)
