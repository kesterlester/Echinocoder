"""Tests for Operation, VectorType, Atom construction and are_negatives."""
import pytest
from symatom import ArgumentSymmetry, Operation, VectorType, Atom, are_negatives

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def dot():
    return Operation("dot", rank=2, odd_parity=False, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

@pytest.fixture
def eps3():
    return Operation("eps3", rank=3, odd_parity=True, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)

@pytest.fixture
def unstructured_op():
    return Operation("foo", rank=2, odd_parity=False, argument_symmetry=ArgumentSymmetry.UNSTRUCTURED)


# ---------------------------------------------------------------------------
# Operation
# ---------------------------------------------------------------------------

def test_operation_valid(dot):
    assert dot.name == "dot"
    assert dot.rank == 2
    assert dot.odd_parity == False          # dot is a scalar (even parity)
    assert dot.parity == +1                 # property: -1 if odd_parity else +1
    assert dot.argument_symmetry == ArgumentSymmetry.SYMMETRIC

def test_operation_invalid_rank():
    with pytest.raises(ValueError, match="rank"):
        Operation("bad", rank=0, odd_parity=False, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

def test_operation_frozen(dot):
    with pytest.raises(Exception):
        dot.rank = 99

def test_operation_tex_none_by_default(dot):
    assert dot.tex is None

def test_operation_tex_stored():
    op = Operation("eps", rank=2, odd_parity=True,
                   argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                   tex=r"\epsilon_{#1#2}")
    assert op.tex == r"\epsilon_{#1#2}"

def test_operation_tex_excluded_from_equality():
    op_a = Operation("eps", rank=2, odd_parity=True,
                     argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                     tex=r"\epsilon_{#1#2}")
    op_b = Operation("eps", rank=2, odd_parity=True,
                     argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                     tex=r"\varepsilon_{#1#2}")
    op_c = Operation("eps", rank=2, odd_parity=True,
                     argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)
    assert op_a == op_b
    assert op_a == op_c

def test_operation_tex_excluded_from_hash():
    op_a = Operation("eps", rank=2, odd_parity=True,
                     argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                     tex=r"\epsilon_{#1#2}")
    op_b = Operation("eps", rank=2, odd_parity=True,
                     argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)
    assert hash(op_a) == hash(op_b)

def test_operation_tex_rank_9_ok():
    op = Operation("big", rank=9, odd_parity=False,
                   argument_symmetry=ArgumentSymmetry.UNSTRUCTURED,
                   tex="#1#2#3#4#5#6#7#8#9")
    assert op.tex is not None

def test_operation_tex_rank_10_raises():
    with pytest.raises(ValueError, match="9 LaTeX parameters"):
        Operation("huge", rank=10, odd_parity=False,
                  argument_symmetry=ArgumentSymmetry.UNSTRUCTURED,
                  tex="#1#2#3#4#5#6#7#8#9#10")

def test_operation_tex_repr():
    op = Operation("eps", rank=2, odd_parity=True,
                   argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                   tex=r"\epsilon_{#1#2}")
    assert r"\epsilon_{#1#2}" in repr(op)

def test_operation_no_tex_repr(dot):
    assert "tex=" not in repr(dot)


# ---------------------------------------------------------------------------
# VectorType
# ---------------------------------------------------------------------------

def test_vector_type_valid():
    g = VectorType("electrons", ("a", "b", "c", "d"))
    assert g.size == 4
    assert g.labels == ("a", "b", "c", "d")

def test_vector_type_empty():
    # Empty groups are allowed: S_0 is the trivial group, 0!=1, C(0,0)=1.
    # They carry meaningful identity information (e.g. muon_count=0 vs jet_count=0).
    g = VectorType("empty", ())
    assert g.size == 0
    assert g.labels == ()

def test_vector_type_duplicate_labels():
    with pytest.raises(ValueError, match="distinct"):
        VectorType("bad", ("a", "a", "b"))

def test_vector_type_requires_tuple():
    with pytest.raises(TypeError):
        VectorType("bad", ["a", "b"])

def test_vector_type_frozen():
    g = VectorType("g", ("a", "b"))
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

def test_atom_sign_minus_one_on_symmetric_op_allowed(dot):
    # sign=-1 is valid for all ArgumentSymmetry values; represents -dot(a,b).
    a = Atom(dot, ("a", "b"), sign=-1)
    assert a.sign == -1
    assert a.labels == ("a", "b")

def test_atom_sign_minus_one_on_symmetric_op_sorts_labels(dot):
    # Label sort still happens; sign is preserved (not perm-parity absorbed).
    a = Atom(dot, ("b", "a"), sign=-1)
    assert a.labels == ("a", "b")
    assert a.sign == -1

def test_atom_sign_minus_one_on_unstructured_op_allowed(unstructured_op):
    # sign=-1 is valid for UNSTRUCTURED too.
    a = Atom(unstructured_op, ("a", "b"), sign=-1)
    assert a.sign == -1

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

def test_are_negatives_symmetric_op_works(dot):
    # SYMMETRIC atoms may carry sign=-1, so are_negatives works for them too.
    a = Atom(dot, ("a", "b"), sign=+1)
    b = Atom(dot, ("a", "b"), sign=-1)
    assert are_negatives(a, b)
    assert are_negatives(b, a)

def test_are_negatives_symmetric_op_same_sign(dot):
    a = Atom(dot, ("a", "b"), sign=+1)
    b = Atom(dot, ("a", "b"), sign=+1)
    assert not are_negatives(a, b)
