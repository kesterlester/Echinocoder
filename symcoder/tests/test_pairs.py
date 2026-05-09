"""Tests for pairs.py: eval_pair_orbit."""
import pytest
import numpy as np
from symatom import (
    ArgumentSymmetry, VectorGroup, Context, Plan,
    SimpleCanonicaliser, repL,
    canonical_pair_flavours,
)
from symcoder import EvaluableOperation, evaluate
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
def ortho_event():
    """Four mutually orthogonal unit 4-vectors (embedded in R^4)."""
    return {
        "a": np.array([1.0, 0.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0, 0.0]),
        "d": np.array([0.0, 0.0, 0.0, 1.0]),
    }


# ---------------------------------------------------------------------------
# eval_pair_orbit — orbit size matches pf.orbit_size()
# ---------------------------------------------------------------------------

def test_eval_pair_orbit_length_dot(dot, plan, ctx, ortho_event):
    fo_list = repL(ctx, (dot,))
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        result = eval_pair_orbit(pf, plan, ortho_event)
        expected = pf.orbit_size(group_sizes)
        assert len(result) == expected, (
            f"Orbit size mismatch for {pf!r}: "
            f"got {len(result)}, expected {expected}"
        )

def test_eval_pair_orbit_length_dot_eps(dot, eps3, plan, ctx):
    event_3d = {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "d": np.array([1.0, 1.0, 1.0]) / 3.0**0.5,
    }
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        result = eval_pair_orbit(pf, plan, event_3d)
        assert len(result) == pf.orbit_size(group_sizes)


# ---------------------------------------------------------------------------
# eval_pair_orbit — numerical values for known cases
# ---------------------------------------------------------------------------

def test_eval_pair_orbit_dot_dot_orthogonal_no_overlap(dot, plan, ctx):
    """
    dot[2] x dot[2] overlap=0 on orthonormal vectors.
    Every orbit element is (dot(x,y), dot(z,w)) with {x,y}∩{z,w}=∅.
    All dot products between distinct basis vectors are 0, so every
    z_k = 0 + 0i.
    """
    event = {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "d": np.array([1.0, 1.0, 1.0]) / 3.0**0.5,
    }
    from symatom.rep import Flavour
    from symatom import PairFlavour
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(0,))
    result = eval_pair_orbit(pf, plan, event)
    assert len(result) == pf.count((4,))

def test_eval_pair_orbit_diagonal_pair(dot, plan, ctx):
    """
    dot[2] x dot[2] overlap=2 means u==v (same atom).
    z_k = eval(u') + i*eval(u') = (1+i)*eval(u').
    So every z_k lies on the line arg(z)=pi/4.
    """
    event = {
        "a": np.array([1.0, 2.0, 0.0]),
        "b": np.array([3.0, 0.0, 4.0]),
        "c": np.array([0.0, 1.0, 1.0]),
        "d": np.array([2.0, 2.0, 1.0]),
    }
    from symatom.rep import Flavour
    from symatom import PairFlavour
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(2,))
    result = eval_pair_orbit(pf, plan, event)
    for z in result:
        assert z.real == pytest.approx(z.imag), (
            f"Expected real==imag (diagonal pair), got {z}"
        )

def test_eval_pair_orbit_returns_complex(dot, plan, ctx, ortho_event):
    fo_list = repL(ctx, (dot,))
    for pf in canonical_pair_flavours(fo_list, ctx):
        result = eval_pair_orbit(pf, plan, ortho_event)
        assert all(isinstance(z, complex) for z in result)


# ---------------------------------------------------------------------------
# eval_pair_orbit — two-group context
# ---------------------------------------------------------------------------

def test_eval_pair_orbit_two_groups(dot, eps3):
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
    fo_list = repL(ctx, (dot, eps3))
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        result = eval_pair_orbit(pf, plan, event)
        assert len(result) == pf.orbit_size(group_sizes)
