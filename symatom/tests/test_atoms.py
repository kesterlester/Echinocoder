"""Tests for Operation, VectorGroup, Atom construction and are_negatives."""
import pytest
from symatom import ArgumentSymmetry, Operation, VectorGroup, Atom, are_negatives

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def dot():
    return Operation("dot", rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

@pytest.fixture
def eps3():
    return Operation("eps3", rank=3, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)

@pytest.fixture
def unstructured_op():
    return Operation("foo", rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.UNSTRUCTURED)


# ---------------------------------------------------------------------------
# Operation
# ---------------------------------------------------------------------------

def test_operation_valid(dot):
    assert dot.name == "dot"
    assert dot.rank == 2
    assert dot.parity == +1
    assert dot.argument_symmetry == ArgumentSymmetry.SYMMETRIC

def test_operation_invalid_rank():
    with pytest.raises(ValueError, match="rank"):
        Operation("bad", rank=0, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

def test_operation_invalid_parity():
    with pytest.raises(ValueError, match="parity"):
        Operation("bad", rank=1, parity=0, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

def test_operation_frozen(dot):
    with pytest.raises(Exception):
        dot.rank = 99


# ---------------------------------------------------------------------------
# VectorGroup
# ---------------------------------------------------------------------------

def test_vector_group_valid():
    g = VectorGroup("electrons", ("a", "b", "c", "d"))
    assert g.size == 4
    assert g.labels == ("a", "b", "c", "d")

def test_vector_group_empty():
    with pytest.raises(ValueError):
        VectorGroup("empty", ())

def test_vector_group_duplicate_labels():
    with pytest.raises(ValueError, match="distinct"):
        VectorGroup("bad", ("a", "a", "b"))

def test_vector_group_requires_tuple():
    with pytest.raises(TypeError):
        VectorGroup("bad", ["a", "b"])

def test_vector_group_frozen():
    g = VectorGroup("g", ("a", "b"))
    with pytest.raises(Exception):
        g.name = "other"


# ---------------------------------------------------------------------------
# Atom — well-formedness
# ---------------------------------------------------------------------------

def test_atom_valid_dot(dot):
    a = Atom(dot, ("a", "b"), sign=+1)
    assert a.sign == +1
    assert a.labels == ("a", "b")

def test_atom_valid_eps3_positive(eps3):
    a = Atom(eps3, ("a", "b", "c"), sign=+1)
    assert a.sign == +1

def test_atom_valid_eps3_negative(eps3):
    a = Atom(eps3, ("a", "b", "c"), sign=-1)
    assert a.sign == -1

def test_atom_wrong_rank(dot):
    with pytest.raises(ValueError, match="rank"):
        Atom(dot, ("a",), sign=+1)

def test_atom_sign_minus_one_on_symmetric_op_raises(dot):
    with pytest.raises(ValueError, match="sign=-1"):
        Atom(dot, ("a", "b"), sign=-1)

def test_atom_sign_minus_one_on_unstructured_op_raises(unstructured_op):
    with pytest.raises(ValueError, match="sign=-1"):
        Atom(unstructured_op, ("a", "b"), sign=-1)

def test_atom_invalid_sign(dot):
    with pytest.raises(ValueError, match="sign"):
        Atom(dot, ("a", "b"), sign=0)

def test_atom_requires_tuple_labels(dot):
    with pytest.raises(TypeError):
        Atom(dot, ["a", "b"], sign=+1)

def test_atom_duplicate_labels_raises(dot, eps3):
    with pytest.raises(ValueError, match="distinct"):
        Atom(dot, ("a", "a"), sign=+1)
    with pytest.raises(ValueError, match="distinct"):
        Atom(eps3, ("a", "a", "b"), sign=+1)

def test_atom_frozen(dot):
    a = Atom(dot, ("a", "b"), sign=+1)
    with pytest.raises(Exception):
        a.sign = -1


# ---------------------------------------------------------------------------
# are_negatives
# ---------------------------------------------------------------------------

def test_are_negatives_true(eps3):
    a = Atom(eps3, ("a", "b", "c"), sign=+1)
    b = Atom(eps3, ("a", "b", "c"), sign=-1)
    assert are_negatives(a, b)
    assert are_negatives(b, a)

def test_are_negatives_same_sign(eps3):
    a = Atom(eps3, ("a", "b", "c"), sign=+1)
    b = Atom(eps3, ("a", "b", "c"), sign=+1)
    assert not are_negatives(a, b)

def test_are_negatives_different_labels(eps3):
    a = Atom(eps3, ("a", "b", "c"), sign=+1)
    b = Atom(eps3, ("a", "b", "d"), sign=-1)
    assert not are_negatives(a, b)

def test_are_negatives_symmetric_op_always_false(dot):
    # Can't even construct a sign=-1 dot, so compare two +1 dots.
    a = Atom(dot, ("a", "b"), sign=+1)
    b = Atom(dot, ("a", "b"), sign=+1)
    assert not are_negatives(a, b)
