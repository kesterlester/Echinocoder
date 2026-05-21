"""Tests for encode.py: the top-level encode() function."""
import pytest
import numpy as np
from symatom import (
    ArgumentSymmetry, VectorType, Context, Plan,
    repS, canonical_pair_flavours,
)
from symcoder import EvaluableOperation, encode, describe_encoding


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
    return VectorType("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def ctx(electrons):
    return Context(types=(electrons,))

@pytest.fixture
def plan(dot, eps3, ctx):
    return Plan(context=ctx, operations=(dot, eps3))

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

def test_encode_returns_ndarray(plan, event_3d, orbit_factory, phase2_factory):
    result = encode(plan, event_3d, orbit_factory, phase2_factory)
    assert isinstance(result, np.ndarray)

def test_encode_returns_real(plan, event_3d, orbit_factory, phase2_factory):
    result = encode(plan, event_3d, orbit_factory, phase2_factory)
    assert result.dtype == np.float64

def test_encode_length_matches_describe(plan, event_3d, orbit_factory, phase2_factory):
    """encode() output length matches sum of segment lengths from describe_encoding()."""
    result = encode(plan, event_3d, orbit_factory, phase2_factory)
    assert len(result) == sum(s.length for s in describe_encoding(plan, orbit_factory, phase2_factory))

def test_encode_dot_only_length(dot, ctx, event_3d, orbit_factory, phase2_factory):
    """Length check for dot-only plan matches describe_encoding()."""
    plan = Plan(context=ctx, operations=(dot,))
    result = encode(plan, event_3d, orbit_factory, phase2_factory)
    assert len(result) == sum(s.length for s in describe_encoding(plan, orbit_factory, phase2_factory))


# ---------------------------------------------------------------------------
# Permutation invariance
# ---------------------------------------------------------------------------

def test_encode_invariant_under_label_permutation(dot, ctx, orbit_factory, phase2_factory):
    """
    Swapping two particle labels must not change the encoding.
    encode(E) == encode(E') where E' is E with labels a and b swapped.
    """
    plan = Plan(context=ctx, operations=(dot,))
    event = {
        "a": np.array([1.0, 2.0, 3.0]),
        "b": np.array([4.0, 5.0, 6.0]),
        "c": np.array([7.0, 0.0, 1.0]),
        "d": np.array([0.0, 3.0, 2.0]),
    }
    event_swapped = dict(event)
    event_swapped["a"], event_swapped["b"] = event["b"], event["a"]

    np.testing.assert_array_almost_equal(
        encode(plan, event, orbit_factory, phase2_factory),
        encode(plan, event_swapped, orbit_factory, phase2_factory),
    )

def test_encode_invariant_under_cyclic_permutation(dot, ctx, orbit_factory, phase2_factory):
    """Cyclic permutation a→b→c→d→a leaves the encoding unchanged."""
    plan = Plan(context=ctx, operations=(dot,))
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
        encode(plan, event, orbit_factory, phase2_factory),
        encode(plan, event_cycled, orbit_factory, phase2_factory),
    )


# ---------------------------------------------------------------------------
# Determinism and reproducibility
# ---------------------------------------------------------------------------

def test_encode_deterministic(plan, event_3d, orbit_factory, phase2_factory):
    """encode() called twice on the same event returns the same vector."""
    r1 = encode(plan, event_3d, orbit_factory, phase2_factory)
    r2 = encode(plan, event_3d, orbit_factory, phase2_factory)
    np.testing.assert_array_equal(r1, r2)


# ---------------------------------------------------------------------------
# Two-group context
# ---------------------------------------------------------------------------

def test_encode_two_groups(dot, eps3, orbit_factory, phase2_factory):
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(dot, eps3))
    event = {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "p": np.array([1.0, 1.0, 0.0]) / 2.0**0.5,
        "q": np.array([0.0, 1.0, 1.0]) / 2.0**0.5,
    }
    result = encode(plan, event, orbit_factory, phase2_factory)
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    assert len(result) > 0

def test_encode_two_groups_permutation_invariant(dot, orbit_factory, phase2_factory):
    """Swapping muon labels p↔q leaves the encoding unchanged."""
    electrons = VectorType("electrons", ("a", "b"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(dot,))
    event = {
        "a": np.array([1.0, 2.0, 0.0]),
        "b": np.array([3.0, 0.0, 1.0]),
        "p": np.array([0.0, 1.0, 4.0]),
        "q": np.array([2.0, 0.0, 1.0]),
    }
    event_swapped = dict(event)
    event_swapped["p"], event_swapped["q"] = event["q"], event["p"]
    np.testing.assert_array_almost_equal(
        encode(plan, event, orbit_factory, phase2_factory),
        encode(plan, event_swapped, orbit_factory, phase2_factory),
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



@pytest.mark.parametrize("perm", [
    {"a": "b", "b": "a", "c": "c", "d": "d"},
    {"a": "b", "b": "c", "c": "d", "d": "a"},
    {"a": "d", "b": "c", "c": "b", "d": "a"},
])
def test_compressed_encoding_permutation_invariant_dot(perm, ctx, orbit_factory, phase2_factory):
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
        encode(plan, event, orbit_factory, phase2_factory),
        encode(plan, event_permuted, orbit_factory, phase2_factory),
    )


@pytest.mark.parametrize("perm", [
    {"a": "b", "b": "a", "c": "c", "d": "d"},
    {"a": "b", "b": "c", "c": "d", "d": "a"},
])
def test_compressed_encoding_permutation_invariant_eps3(perm, ctx, orbit_factory, phase2_factory):
    # TODO Make this a much better test that looks at cases that can't be compressed too, e.g. eps3(a,p,v)
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
        encode(plan, event, orbit_factory, phase2_factory),
        encode(plan, event_permuted, orbit_factory, phase2_factory),
    )


def test_compressed_length_sym_sym(dot, ctx, event_3d):
    """SYM×SYM: compressed length == orbit_size (no compression)."""
    plan = Plan(context=ctx, operations=(dot,))
    fo_list = repS(ctx, (dot,))
    type_sizes = tuple(g.size for g in ctx.types)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert pf.count(type_sizes) == pf.orbit_size(ctx)



def test_compressed_length_sym_antisym(ctx):
    # TODO FIXME : This is wrong, referring to ANTISYMM etc rather than TYPE_11,12,21,22
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
    fo_list = repS(ctx, (dot, eps3))
    type_sizes = tuple(g.size for g in ctx.types)
    mixed_pfs = [
        pf for pf in canonical_pair_flavours(fo_list, ctx)
        if (pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
        != (pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
    ]
    assert len(mixed_pfs) > 0, "Need at least one SYM×ANTISYM pair for this test"
    for pf in mixed_pfs:
        assert pf.count(type_sizes) == pf.orbit_size(ctx) // 2


# ---------------------------------------------------------------------------
# TYPE_NEG: correlated-negation embedding in a 2-group context
# ---------------------------------------------------------------------------

def _make_eps2():
    return EvaluableOperation(
        name="eps2", rank=2, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(vecs[0][0]*vecs[1][1] - vecs[0][1]*vecs[1][0]),
    )


def test_encode_two_groups_permutation_invariant_eps2(orbit_factory, phase2_factory):
    """encode() with an ANTISYMMETRIC rank-2 op in a 2-group context is permutation-invariant.

    This exercises the TYPE_NEG branch (AA pairs with full electron overlap) and
    ensures the full encode pipeline handles it correctly.
    """
    eps2 = _make_eps2()
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(eps2,))
    event = {
        "a": np.array([1.3, 0.2]),
        "b": np.array([0.7, 1.1]),
        "c": np.array([0.4, 0.9]),
        "p": np.array([0.8, 0.3]),
        "q": np.array([0.1, 1.4]),
    }
    event_swapped = dict(event)
    event_swapped["a"], event_swapped["b"] = event["b"], event["a"]
    np.testing.assert_array_almost_equal(
        encode(plan, event, orbit_factory, phase2_factory),
        encode(plan, event_swapped, orbit_factory, phase2_factory),
        err_msg="encode() not invariant under a↔b swap in 2-group eps2 context"
    )
