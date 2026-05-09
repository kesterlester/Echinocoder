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
# Internal argument canonicalization
# ---------------------------------------------------------------------------

def test_symmetric_op_sorts_labels(dot):
    a = Atom(dot, ("b", "a"), sign=+1)
    assert a.labels == ("a", "b")
    assert a.sign == +1

def test_symmetric_op_equal_regardless_of_input_order(dot):
    assert Atom(dot, ("b", "a"), sign=+1) == Atom(dot, ("a", "b"), sign=+1)

def test_antisymmetric_op_sorts_labels_even_permutation(eps3):
    # (c,a,b) → (a,b,c) requires 2 swaps (even permutation): sign unchanged
    a = Atom(eps3, ("c", "a", "b"), sign=+1)
    assert a.labels == ("a", "b", "c")
    assert a.sign == +1

def test_antisymmetric_op_sorts_labels_odd_permutation(eps3):
    # (b,a,c) → (a,b,c) requires 1 swap (odd permutation): sign flips
    a = Atom(eps3, ("b", "a", "c"), sign=+1)
    assert a.labels == ("a", "b", "c")
    assert a.sign == -1

def test_antisymmetric_op_sign_tracks_parity(eps3):
    # Supplying sign=-1 with an odd input permutation: -1 * -1 = +1
    a = Atom(eps3, ("b", "a", "c"), sign=-1)
    assert a.labels == ("a", "b", "c")
    assert a.sign == +1

def test_antisymmetric_are_negatives_unsorted_input(eps3):
    # (b,a,c)+1 self-canonicalises to (a,b,c)-1; (a,b,c)+1 stays as-is.
    x = Atom(eps3, ("b", "a", "c"), sign=+1)
    y = Atom(eps3, ("a", "b", "c"), sign=+1)
    assert are_negatives(x, y)

def test_unstructured_op_preserves_order(unstructured_op):
    # UNSTRUCTURED atoms are stored in the order given.
    a = Atom(unstructured_op, ("b", "a"), sign=+1)
    assert a.labels == ("b", "a")


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
