"""Tests for EvaluableOperation and evaluate()."""
import math
import pytest
import numpy as np
from symatom import ArgumentSymmetry, Atom, Context, VectorType
from symcoder import EvaluableOperation, evaluate


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
def event_3d():
    """Three orthonormal 3-vectors labelled a, b, c."""
    return {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
    }


# ---------------------------------------------------------------------------
# EvaluableOperation is a valid Operation
# ---------------------------------------------------------------------------

def test_evaluable_operation_is_operation(dot):
    from symatom.atoms import Operation
    assert isinstance(dot, Operation)

def test_evaluable_operation_repr(dot):
    r = repr(dot)
    assert "dot" in r
    assert "rank=2" in r
    assert "SYMMETRIC" in r

def test_evaluable_operations_equal_by_symbolic_fields(dot):
    dot2 = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: 999.0,   # different callable
    )
    assert dot == dot2

def test_evaluable_operation_rank_validation():
    with pytest.raises(ValueError):
        EvaluableOperation(
            name="bad", rank=0, parity=+1,
            argument_symmetry=ArgumentSymmetry.SYMMETRIC,
            eval_fn=lambda vecs: 0.0,
        )


# ---------------------------------------------------------------------------
# evaluate() — dot product
# ---------------------------------------------------------------------------

def test_evaluate_dot_orthogonal(dot, event_3d):
    atom = Atom(dot, ("a", "b"), sign=+1)
    assert evaluate(atom, event_3d) == pytest.approx(0.0)

def test_evaluate_dot_parallel(dot, event_3d):
    # Atom constructor raises on duplicate labels
    with pytest.raises(ValueError):
        Atom(dot, ("a", "a"), sign=+1)

def test_evaluate_dot_self_via_rank1():
    # dot(a,a) is ill-formed; use a rank-1 squared-length op instead
    lenSq = EvaluableOperation(
        name="lenSq", rank=1, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[0])),
    )
    event = {"a": np.array([3.0, 4.0, 0.0])}
    atom = Atom(lenSq, ("a",), sign=+1)
    assert evaluate(atom, event) == pytest.approx(25.0)

def test_evaluate_dot_value(dot, event_3d):
    event = {"a": np.array([1.0, 2.0, 0.0]), "b": np.array([3.0, 4.0, 0.0])}
    atom = Atom(dot, ("a", "b"), sign=+1)
    assert evaluate(atom, event) == pytest.approx(11.0)   # 1*3 + 2*4


# ---------------------------------------------------------------------------
# evaluate() — eps3 (scalar triple product)
# ---------------------------------------------------------------------------

def test_evaluate_eps3_right_handed_frame(eps3, event_3d):
    # eps3(a,b,c) on right-handed orthonormal frame = +1
    atom = Atom(eps3, ("a", "b", "c"), sign=+1)
    assert evaluate(atom, event_3d) == pytest.approx(1.0)

def test_evaluate_eps3_sign_applied(eps3, event_3d):
    # The atom with sign=-1 gives -1 * (triple product)
    atom = Atom(eps3, ("a", "b", "c"), sign=-1)
    assert evaluate(atom, event_3d) == pytest.approx(-1.0)

def test_evaluate_eps3_argument_sort_sign(eps3, event_3d):
    # Atom(eps3, ("b","a","c"), +1) auto-sorts to ("a","b","c") with sign -1.
    # evaluate should return -1 * triple_product = -1.
    atom = Atom(eps3, ("b", "a", "c"), sign=+1)
    assert atom.labels == ("a", "b", "c")
    assert atom.sign == -1
    assert evaluate(atom, event_3d) == pytest.approx(-1.0)


# ---------------------------------------------------------------------------
# evaluate() — error on plain Operation
# ---------------------------------------------------------------------------

def test_evaluate_raises_on_plain_operation(event_3d):
    from symatom.atoms import Operation
    plain_dot = Operation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
    )
    atom = Atom(plain_dot, ("a", "b"), sign=+1)
    with pytest.raises(TypeError, match="eval_fn"):
        evaluate(atom, event_3d)
