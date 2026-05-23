"""
Tests verifying that symcoder.encode() handles contexts with empty VectorTypes.
"""
import numpy as np
import pytest
from symatom import ArgumentSymmetry, VectorType, Context, Plan
from symcoder import EvaluableOperation, encode


@pytest.fixture
def dot():
    return EvaluableOperation(
        name="dot", rank=2, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )

@pytest.fixture
def event():
    return {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
    }


def test_encode_with_empty_group_returns_ndarray(dot, event, orbit_factory, phase2_factory):
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ())
    ctx  = Context((electrons, muons))
    plan = Plan(context=ctx, operations=(dot,))
    result = encode(plan, event, orbit_factory, phase2_factory)
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    assert len(result) > 0

def test_encode_with_empty_group_same_as_without(dot, event, orbit_factory, phase2_factory):
    """
    Encoding with an empty muons group gives the same result as encoding
    with no muons group at all (the empty group contributes nothing).
    """
    electrons = VectorType("electrons", ("a", "b", "c"))
    ctx_with    = Context((electrons, VectorType("muons", ())))
    ctx_without = Context((electrons,))
    plan_with    = Plan(context=ctx_with,    operations=(dot,))
    plan_without = Plan(context=ctx_without, operations=(dot,))
    np.testing.assert_array_equal(
        encode(plan_with,    event, orbit_factory, phase2_factory),
        encode(plan_without, event, orbit_factory, phase2_factory),
    )

def test_encode_with_empty_group_permutation_invariant(dot, event, orbit_factory, phase2_factory):
    """Swapping electron labels leaves the encoding unchanged."""
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ())
    ctx  = Context((electrons, muons))
    plan = Plan(context=ctx, operations=(dot,))
    event_swapped = dict(event)
    event_swapped["a"], event_swapped["b"] = event["b"], event["a"]
    np.testing.assert_array_almost_equal(
        encode(plan, event,         orbit_factory, phase2_factory),
        encode(plan, event_swapped, orbit_factory, phase2_factory),
    )

def test_encode_all_empty_groups_returns_empty(dot, orbit_factory, phase2_factory):
    """With all groups empty there are no atoms, so encode returns an empty array."""
    ctx  = Context((VectorType("electrons", ()), VectorType("muons", ())))
    plan = Plan(context=ctx, operations=(dot,))
    result = encode(plan, {}, orbit_factory, phase2_factory)
    assert isinstance(result, np.ndarray)
    assert len(result) == 0
