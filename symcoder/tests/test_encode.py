"""Tests for encode.py: the top-level encode() function."""
import pytest
import numpy as np
from symatom import (
    ArgumentSymmetry, VectorGroup, Context, Plan, SimpleCanonicaliser,
    repL, canonical_pair_flavours,
)
from symcoder import EvaluableOperation, encode, encode_brute, describe_encoding
from symcoder.pairs import eval_pair_orbit, eval_pair_orbit_positive
from symcoder.encode import _embed_compressed


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

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
def electrons():
    return VectorGroup("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def ctx(electrons):
    return Context(groups=(electrons,))

@pytest.fixture
def plan(dot, eps3, ctx):
    return Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot, eps3))

@pytest.fixture
def event_3d():
    return {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "d": np.array([1.0, 1.0, 1.0]) / 3.0**0.5,
    }


# ---------------------------------------------------------------------------
# Basic shape and type checks
# ---------------------------------------------------------------------------

def test_encode_returns_ndarray(dot, plan, ctx, event_3d):
    result = encode(plan, event_3d)
    assert isinstance(result, np.ndarray)

def test_encode_returns_real(dot, plan, ctx, event_3d):
    result = encode(plan, event_3d)
    assert result.dtype == np.float64

def test_encode_length_matches_describe(dot, eps3, plan, ctx, event_3d):
    """encode() output length matches sum of segment lengths from describe_encoding()."""
    result = encode(plan, event_3d)
    assert len(result) == sum(s.length for s in describe_encoding(plan))

def test_encode_dot_only_length(dot, ctx, event_3d):
    """Length check for dot-only plan matches describe_encoding()."""
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))
    result = encode(plan, event_3d)
    assert len(result) == sum(s.length for s in describe_encoding(plan))

def test_encode_brute_length_matches_sum_of_counts(dot, eps3, plan, ctx, event_3d):
    """encode_brute() gives 2*sum(pf.count()) reals — n complex coeffs × 2 reals each."""
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    expected_len = 2 * sum(
        pf.count(group_sizes)
        for pf in canonical_pair_flavours(fo_list, ctx)
    )
    assert len(encode_brute(plan, event_3d)) == expected_len


# ---------------------------------------------------------------------------
# Permutation invariance
# ---------------------------------------------------------------------------

def test_encode_invariant_under_label_permutation(dot, ctx):
    """
    Swapping two particle labels must not change the encoding.
    encode(E) == encode(E') where E' is E with labels a and b swapped.
    """
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))
    event = {
        "a": np.array([1.0, 2.0, 3.0]),
        "b": np.array([4.0, 5.0, 6.0]),
        "c": np.array([7.0, 0.0, 1.0]),
        "d": np.array([0.0, 3.0, 2.0]),
    }
    event_swapped = dict(event)
    event_swapped["a"], event_swapped["b"] = event["b"], event["a"]

    result       = encode(plan, event)
    result_swapped = encode(plan, event_swapped)
    np.testing.assert_array_almost_equal(result, result_swapped)

def test_encode_invariant_under_cyclic_permutation(dot, ctx):
    """Cyclic permutation a→b→c→d→a leaves the encoding unchanged."""
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))
    event = {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "d": np.array([1.0, 1.0, 0.0]) / 2.0**0.5,
    }
    event_cycled = {
        "a": event["b"],
        "b": event["c"],
        "c": event["d"],
        "d": event["a"],
    }
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_cycled),
    )


# ---------------------------------------------------------------------------
# Determinism and reproducibility
# ---------------------------------------------------------------------------

def test_encode_deterministic(dot, plan, ctx, event_3d):
    """encode() called twice on the same event returns the same vector."""
    r1 = encode(plan, event_3d)
    r2 = encode(plan, event_3d)
    np.testing.assert_array_equal(r1, r2)


# ---------------------------------------------------------------------------
# Two-group context
# ---------------------------------------------------------------------------

def test_encode_two_groups(dot, eps3):
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx  = Context(groups=(electrons, muons))
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot, eps3))
    event = {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "p": np.array([1.0, 1.0, 0.0]) / 2.0**0.5,
        "q": np.array([0.0, 1.0, 1.0]) / 2.0**0.5,
    }
    result = encode(plan, event)
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    assert len(result) > 0

def test_encode_two_groups_permutation_invariant(dot):
    """Swapping muon labels p↔q leaves the encoding unchanged."""
    electrons = VectorGroup("electrons", ("a", "b"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx  = Context(groups=(electrons, muons))
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))
    event = {
        "a": np.array([1.0, 2.0, 0.0]),
        "b": np.array([3.0, 0.0, 1.0]),
        "p": np.array([0.0, 1.0, 4.0]),
        "q": np.array([2.0, 0.0, 1.0]),
    }
    event_swapped = dict(event)
    event_swapped["p"], event_swapped["q"] = event["q"], event["p"]
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_swapped),
    )


# ---------------------------------------------------------------------------
# Optimization 5a: compressed polynomial embedding
# ---------------------------------------------------------------------------

def _make_ops_and_event():
    """Shared setup: four 3-D vectors, all four symmetry-class operation pairs."""
    labels = ("a", "b", "c", "d")
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    event = {
        "a": np.array([1.0, 0.2, 0.3]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.3, 1.0]),
        "d": np.array([0.2, 0.7, 0.1]),
    }
    return dot, eps3, labels, event


def test_eval_pair_orbit_positive_count(dot, eps3, ctx):
    """eval_pair_orbit_positive always returns exactly pf.count() values."""
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        plan = Plan(context=ctx, operations=(dot, eps3))
        event = {l: np.random.randn(3) for l in ctx.all_labels}
        pos = eval_pair_orbit_positive(pf, plan, event)
        assert len(pos) == pf.count(group_sizes), (
            f"Expected {pf.count(group_sizes)}, got {len(pos)} for {pf!r}"
        )


def test_eval_pair_orbit_positive_count_unusual_label_order():
    """Count invariant holds even with unusual alphabetical label ordering."""
    # "toast" < "apple" < "zebra" under default str ordering is False,
    # but "apple" < "toast" < "zebra" — the point is the Atom constructor
    # sorts internally and absorbs parity into sign, so count is always n.
    labels = ("apple", "toast", "zebra")
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    g = VectorGroup("particles", labels)
    ctx = Context((g,))
    plan = Plan(context=ctx, operations=(dot, eps3))
    event = {l: np.random.randn(3) for l in labels}
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        pos = eval_pair_orbit_positive(pf, plan, event)
        assert len(pos) == pf.count(group_sizes), (
            f"Expected {pf.count(group_sizes)}, got {len(pos)} for {pf!r}"
        )


@pytest.mark.parametrize("perm", [
    {"a": "b", "b": "a", "c": "c", "d": "d"},
    {"a": "b", "b": "c", "c": "d", "d": "a"},
    {"a": "d", "b": "c", "c": "b", "d": "a"},
])
def test_compressed_encoding_permutation_invariant_dot(perm, ctx):
    """Compressed encoding of dot-only plan is invariant under label permutation."""
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    plan = Plan(context=ctx, operations=(dot,))
    event = {
        "a": np.array([1.0, 0.2, 0.3]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.3, 1.0]),
        "d": np.array([0.2, 0.7, 0.1]),
    }
    event_permuted = {perm[k]: v for k, v in event.items()}
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_permuted),
    )


@pytest.mark.parametrize("perm", [
    {"a": "b", "b": "a", "c": "c", "d": "d"},
    {"a": "b", "b": "c", "c": "d", "d": "a"},
])
def test_compressed_encoding_permutation_invariant_eps3(perm, ctx):
    """Compressed encoding with ANTISYMMETRIC eps3 is permutation-invariant."""
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    plan = Plan(context=ctx, operations=(eps3,))
    event = {
        "a": np.array([1.0, 0.2, 0.3]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.3, 1.0]),
        "d": np.array([0.2, 0.7, 0.1]),
    }
    event_permuted = {perm[k]: v for k, v in event.items()}
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_permuted),
    )


def test_compressed_length_sym_sym(dot, ctx, event_3d):
    """SYM×SYM: compressed length == orbit_size (no compression)."""
    plan = Plan(context=ctx, operations=(dot,))
    fo_list = repL(ctx, (dot,))
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert pf.count(group_sizes) == pf.orbit_size(group_sizes)


def test_compressed_length_antisym_antisym(ctx, event_3d):
    """ANTISYM×ANTISYM: compressed length == orbit_size / 4."""
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    plan = Plan(context=ctx, operations=(eps3,))
    fo_list = repL(ctx, (eps3,))
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert pf.count(group_sizes) == pf.orbit_size(group_sizes) // 4


def test_compressed_length_sym_antisym(ctx):
    """SYM×ANTISYM or ANTISYM×SYM: compressed length == orbit_size / 2."""
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    plan = Plan(context=ctx, operations=(dot, eps3))
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    mixed_pfs = [
        pf for pf in canonical_pair_flavours(fo_list, ctx)
        if (pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
        != (pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
    ]
    assert len(mixed_pfs) > 0, "Need at least one SYM×ANTISYM pair for this test"
    for pf in mixed_pfs:
        assert pf.count(group_sizes) == pf.orbit_size(group_sizes) // 2
