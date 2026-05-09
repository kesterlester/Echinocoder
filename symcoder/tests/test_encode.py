"""Tests for encode.py: the top-level encode() function."""
import pytest
import numpy as np
from symatom import (
    ArgumentSymmetry, VectorGroup, Context, Plan, SimpleCanonicaliser,
    repL, canonical_pair_flavours,
)
from symcoder import EvaluableOperation, encode
from symcoder.pairs import eval_pair_orbit


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

def test_encode_returns_complex(dot, plan, ctx, event_3d):
    result = encode(plan, event_3d)
    assert np.iscomplexobj(result)

def test_encode_length_matches_sum_of_orbit_sizes(dot, eps3, plan, ctx, event_3d):
    """
    Each PairFlavour contributes exactly orbit_size(group_sizes) coefficients
    (the zip polynomial has that degree; the leading 1 is dropped).
    """
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    expected_len = sum(
        pf.orbit_size(group_sizes)
        for pf in canonical_pair_flavours(fo_list, ctx)
    )
    result = encode(plan, event_3d)
    assert len(result) == expected_len

def test_encode_dot_only_length(dot, ctx, event_3d):
    """Same length check for dot-only plan."""
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))
    fo_list = repL(ctx, (dot,))
    group_sizes = tuple(g.size for g in ctx.groups)
    expected_len = sum(
        pf.orbit_size(group_sizes)
        for pf in canonical_pair_flavours(fo_list, ctx)
    )
    result = encode(plan, event_3d)
    assert len(result) == expected_len


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
    assert np.iscomplexobj(result)
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
