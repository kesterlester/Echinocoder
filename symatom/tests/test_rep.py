"""Tests for Flavour, FlavouredOperator, repL, repS."""
import math
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Atom, Context,
)
from symatom.rep import Flavour, FlavouredOperator, repL, repS


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mass():
    return Operation("mass", rank=1, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

@pytest.fixture
def dot():
    return Operation("dot", rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

@pytest.fixture
def eps3():
    return Operation("eps3", rank=3, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)

@pytest.fixture
def electrons():
    return VectorGroup("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def muons():
    return VectorGroup("muons", ("p", "q"))

@pytest.fixture
def ctx1(electrons):
    """Single-group context: 4 electrons."""
    return Context((electrons,))

@pytest.fixture
def jets():
    return VectorGroup("jets", ("u", "v", "w"))

@pytest.fixture
def ctx2(electrons, muons):
    """Two-group context: 4 electrons + 2 muons."""
    return Context((electrons, muons))

@pytest.fixture
def ctx3(electrons, muons, jets):
    """Three-group context: 4 electrons + 2 muons + 3 jets."""
    return Context((electrons, muons, jets))


# ---------------------------------------------------------------------------
# Flavour
# ---------------------------------------------------------------------------

def test_flavour_valid():
    f = Flavour((2, 1))
    assert f.counts == (2, 1)
    assert f.rank == 3
    assert len(f) == 2
    assert list(f) == [2, 1]

def test_flavour_zero_counts():
    f = Flavour((3, 0))
    assert f.rank == 3

def test_flavour_requires_tuple():
    with pytest.raises(TypeError):
        Flavour([2, 1])

def test_flavour_no_negative_counts():
    with pytest.raises(ValueError):
        Flavour((2, -1))


# ---------------------------------------------------------------------------
# FlavouredOperator — construction and validation
# ---------------------------------------------------------------------------

def test_flavoured_operator_valid(dot, ctx1):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1, signed=False)
    assert fo.operation is dot
    assert fo.flavour == Flavour((2,))

def test_flavoured_operator_wrong_group_count(dot, ctx1):
    with pytest.raises(ValueError, match="groups"):
        FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx1, signed=False)

def test_flavoured_operator_wrong_rank(dot, ctx1):
    with pytest.raises(ValueError, match="rank"):
        FlavouredOperator(operation=dot, flavour=Flavour((3,)), context=ctx1, signed=False)

def test_flavoured_operator_count_exceeds_group(dot):
    # A group of size 1 cannot supply 2 labels for a rank-2 dot.
    tiny_ctx = Context((VectorGroup("tiny", ("x",)),))
    with pytest.raises(ValueError, match="size"):
        FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=tiny_ctx, signed=False)


# ---------------------------------------------------------------------------
# FlavouredOperator.count()
# ---------------------------------------------------------------------------

def test_count_mass_single_group(mass, ctx1):
    # C(4,1) = 4; symmetric so repL == repS
    fo = FlavouredOperator(operation=mass, flavour=Flavour((1,)), context=ctx1, signed=False)
    assert fo.count() == 4

def test_count_dot_single_group(dot, ctx1):
    # C(4,2) = 6
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1, signed=False)
    assert fo.count() == 6

def test_count_eps_repS(eps3, ctx1):
    # C(4,3) = 4; repS → no sign doubling
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1, signed=False)
    assert fo.count() == 4

def test_count_eps_repL(eps3, ctx1):
    # C(4,3) = 4; repL → doubled for ANTISYMMETRIC
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1, signed=True)
    assert fo.count() == 8

def test_count_dot_two_groups_mixed_flavour(dot, ctx2):
    # C(4,1) * C(2,1) = 4 * 2 = 8; symmetric so repL == repS
    fo = FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx2, signed=True)
    assert fo.count() == 8

def test_count_eps_two_groups_repL(eps3, ctx2):
    # flavour (2,1): C(4,2) * C(2,1) = 6 * 2 = 12; repL doubles → 24
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((2, 1)), context=ctx2, signed=True)
    assert fo.count() == 24


# ---------------------------------------------------------------------------
# FlavouredOperator.atoms()
# ---------------------------------------------------------------------------

def test_atoms_mass_single_group(mass, ctx1):
    fo = FlavouredOperator(operation=mass, flavour=Flavour((1,)), context=ctx1, signed=False)
    atoms = list(fo.atoms())
    assert len(atoms) == 4
    labels_seen = {a.labels for a in atoms}
    assert labels_seen == {("a",), ("b",), ("c",), ("d",)}
    assert all(a.sign == +1 for a in atoms)

def test_atoms_dot_single_group(dot, ctx1):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1, signed=False)
    atoms = list(fo.atoms())
    assert len(atoms) == 6
    # Every combination of 2 from {a,b,c,d} appears exactly once
    from itertools import combinations
    expected = {tuple(pair) for pair in combinations(("a", "b", "c", "d"), 2)}
    assert {a.labels for a in atoms} == expected

def test_atoms_eps_repS_single_group(eps3, ctx1):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1, signed=False)
    atoms = list(fo.atoms())
    assert len(atoms) == 4
    assert all(a.sign == +1 for a in atoms)

def test_atoms_eps_repL_single_group(eps3, ctx1):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1, signed=True)
    atoms = list(fo.atoms())
    assert len(atoms) == 8
    # Each label combination appears twice: once +1, once -1
    from itertools import combinations
    expected_labels = {tuple(t) for t in combinations(("a", "b", "c", "d"), 3)}
    pos = [a for a in atoms if a.sign == +1]
    neg = [a for a in atoms if a.sign == -1]
    assert {a.labels for a in pos} == expected_labels
    assert {a.labels for a in neg} == expected_labels

def test_atoms_count_matches_count_method(dot, eps3, ctx2):
    for op, flavour, signed in [
        (dot,  Flavour((2, 0)), False),
        (dot,  Flavour((1, 1)), True),
        (eps3, Flavour((2, 1)), False),
        (eps3, Flavour((2, 1)), True),
    ]:
        fo = FlavouredOperator(operation=op, flavour=flavour, context=ctx2, signed=signed)
        assert len(list(fo.atoms())) == fo.count()

def test_atoms_labels_distinct(eps3, ctx2):
    # All generated atoms must satisfy Rule 4 (distinct labels).
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((2, 1)), context=ctx2, signed=True)
    for atom in fo.atoms():
        assert len(set(atom.labels)) == len(atom.labels)


# ---------------------------------------------------------------------------
# FlavouredOperator.contains()
# ---------------------------------------------------------------------------

def test_contains_dot_true(dot, ctx1):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1, signed=False)
    assert fo.contains(Atom(dot, ("a", "b"), sign=+1))

def test_contains_dot_wrong_operation(mass, dot, ctx1):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1, signed=False)
    assert not fo.contains(Atom(mass, ("a",), sign=+1))

def test_contains_eps_repS_positive(eps3, ctx1):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1, signed=False)
    assert fo.contains(Atom(eps3, ("a", "b", "c"), sign=+1))

def test_contains_eps_repS_rejects_negative(eps3, ctx1):
    # repS does not include sign=-1 atoms
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1, signed=False)
    assert not fo.contains(Atom(eps3, ("a", "b", "c"), sign=-1))

def test_contains_eps_repL_accepts_negative(eps3, ctx1):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1, signed=True)
    assert fo.contains(Atom(eps3, ("a", "b", "c"), sign=-1))

def test_contains_wrong_flavour(dot, ctx2):
    # dot(a,b) has flavour (2,0); FlavouredOperator has flavour (1,1)
    fo = FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx2, signed=False)
    assert not fo.contains(Atom(dot, ("a", "b"), sign=+1))

def test_contains_cross_group_dot(dot, ctx2):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx2, signed=False)
    assert fo.contains(Atom(dot, ("a", "p"), sign=+1))
    assert fo.contains(Atom(dot, ("b", "q"), sign=+1))  # label order doesn't matter


# ---------------------------------------------------------------------------
# repL and repS functions
# ---------------------------------------------------------------------------

def test_repL_single_group_single_op(eps3, ctx1):
    # eps3 rank 3, 4 electrons: only flavour (3,) is valid → one FlavouredOperator
    fos = repL(ctx1, [eps3])
    assert len(fos) == 1
    assert fos[0].signed is True
    assert fos[0].flavour == Flavour((3,))

def test_repS_single_group_single_op(eps3, ctx1):
    fos = repS(ctx1, [eps3])
    assert len(fos) == 1
    assert fos[0].signed is False

def test_repL_two_groups_eps3(eps3, ctx2):
    # eps3 rank 3, groups Elecs(4) + Muons(2):
    # valid flavours: (3,0), (2,1), (1,2) — (0,3) excluded since Muons has only 2
    fos = repL(ctx2, [eps3])
    flavours = {fo.flavour for fo in fos}
    assert flavours == {Flavour((3, 0)), Flavour((2, 1)), Flavour((1, 2))}
    assert all(fo.signed is True for fo in fos)

def test_repS_two_groups_eps3(eps3, ctx2):
    fos = repS(ctx2, [eps3])
    flavours = {fo.flavour for fo in fos}
    assert flavours == {Flavour((3, 0)), Flavour((2, 1)), Flavour((1, 2))}
    assert all(fo.signed is False for fo in fos)

def test_repL_two_groups_multiple_ops(mass, dot, eps3, ctx2):
    fos = repL(ctx2, [mass, dot, eps3])
    # mass rank 1: flavours (1,0), (0,1) → 2
    # dot  rank 2: flavours (2,0), (1,1), (0,2) → 3
    # eps3 rank 3: flavours (3,0), (2,1), (1,2) → 3
    assert len(fos) == 8

def test_repL_total_atom_count(eps3, ctx1):
    # repL with one eps3 and 4 electrons: C(4,3)*2 = 8 atoms total
    fos = repL(ctx1, [eps3])
    total = sum(fo.count() for fo in fos)
    assert total == 8

def test_repS_total_atom_count(eps3, ctx1):
    # repS: C(4,3) = 4 atoms total
    fos = repS(ctx1, [eps3])
    total = sum(fo.count() for fo in fos)
    assert total == 4

def test_repL_repS_same_length(mass, dot, eps3, ctx2):
    # repL and repS always produce the same number of FlavouredOperators
    assert len(repL(ctx2, [mass, dot, eps3])) == len(repS(ctx2, [mass, dot, eps3]))


# ---------------------------------------------------------------------------
# Three-group context: Electrons(4) + Muons(2) + Jets(3)
#
# Expected FlavouredOperator counts (one per valid flavour):
#   mass  rank 1:  (1,0,0) (0,1,0) (0,0,1)                          → 3
#   dot   rank 2:  (2,0,0) (1,1,0) (1,0,1) (0,2,0) (0,1,1) (0,0,2) → 6
#   eps3  rank 3:  (3,0,0) (2,1,0) (2,0,1) (1,2,0) (1,1,1)          → 5
#                  (1,0,2) (0,2,1) (0,1,2) (0,0,3)                   → 4
#                                                             total eps3 → 9
#   total: 3 + 6 + 9 = 18
#
# Expected repS atom counts per FlavouredOperator (C(n_i, k_i) product):
#   mass:  C(4,1)=4,  C(2,1)=2,  C(3,1)=3                  → total 9
#   dot:   C(4,2)=6,  C(4,1)C(2,1)=8, C(4,1)C(3,1)=12,
#          C(2,2)=1,  C(2,1)C(3,1)=6, C(3,2)=3             → total 36
#   eps3:  C(4,3)=4,  C(4,2)C(2,1)=12, C(4,2)C(3,1)=18,
#          C(4,1)C(2,2)=4, C(4,1)C(2,1)C(3,1)=24,
#          C(4,1)C(3,2)=12, C(2,2)C(3,1)=3,
#          C(2,1)C(3,2)=6,  C(3,3)=1                        → total 84
#   repS grand total: 9 + 36 + 84 = 129
#   repL grand total: 9 + 36 + 168 = 213  (eps3 atoms doubled)
# ---------------------------------------------------------------------------

def test_three_group_fo_count(mass, dot, eps3, ctx3):
    assert len(repL(ctx3, [mass, dot, eps3])) == 18

def test_three_group_repL_repS_same_fo_count(mass, dot, eps3, ctx3):
    assert len(repL(ctx3, [mass, dot, eps3])) == len(repS(ctx3, [mass, dot, eps3]))

def test_three_group_mass_fo_count(mass, ctx3):
    fos = repS(ctx3, [mass])
    assert len(fos) == 3
    flavours = {fo.flavour for fo in fos}
    assert flavours == {Flavour((1,0,0)), Flavour((0,1,0)), Flavour((0,0,1))}

def test_three_group_dot_fo_count(dot, ctx3):
    fos = repS(ctx3, [dot])
    assert len(fos) == 6

def test_three_group_eps3_fo_count(eps3, ctx3):
    # (0,0,3) excluded because Jets has only 3 labels — exactly valid.
    # (0,3,0) excluded because Muons has only 2 labels.
    fos = repS(ctx3, [eps3])
    assert len(fos) == 9

def test_three_group_repS_atom_totals(mass, dot, eps3, ctx3):
    mass_total = sum(fo.count() for fo in repS(ctx3, [mass]))
    dot_total  = sum(fo.count() for fo in repS(ctx3, [dot]))
    eps3_total = sum(fo.count() for fo in repS(ctx3, [eps3]))
    assert mass_total == 9
    assert dot_total  == 36
    assert eps3_total == 84

def test_three_group_repL_atom_totals(mass, dot, eps3, ctx3):
    # mass and dot are SYMMETRIC — repL count equals repS count.
    # eps3 is ANTISYMMETRIC — repL count is double repS.
    mass_total = sum(fo.count() for fo in repL(ctx3, [mass]))
    dot_total  = sum(fo.count() for fo in repL(ctx3, [dot]))
    eps3_total = sum(fo.count() for fo in repL(ctx3, [eps3]))
    assert mass_total == 9
    assert dot_total  == 36
    assert eps3_total == 168

def test_three_group_grand_total_repS(mass, dot, eps3, ctx3):
    assert sum(fo.count() for fo in repS(ctx3, [mass, dot, eps3])) == 129

def test_three_group_grand_total_repL(mass, dot, eps3, ctx3):
    assert sum(fo.count() for fo in repL(ctx3, [mass, dot, eps3])) == 213

def test_three_group_atoms_matches_count(mass, dot, eps3, ctx3):
    # For every FlavouredOperator, atoms() must yield exactly count() atoms.
    for fo in repL(ctx3, [mass, dot, eps3]):
        assert len(list(fo.atoms())) == fo.count()

def test_three_group_eps3_mixed_flavour_spot_check(eps3, ctx3):
    # Flavour (1,1,1): C(4,1)*C(2,1)*C(3,1) = 4*2*3 = 24 atoms in repS, 48 in repL.
    fo_s = FlavouredOperator(operation=eps3, flavour=Flavour((1,1,1)), context=ctx3, signed=False)
    fo_l = FlavouredOperator(operation=eps3, flavour=Flavour((1,1,1)), context=ctx3, signed=True)
    assert fo_s.count() == 24
    assert fo_l.count() == 48
    assert len(list(fo_s.atoms())) == 24
    assert len(list(fo_l.atoms())) == 48

def test_three_group_atoms_labels_distinct(eps3, ctx3):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((1,1,1)), context=ctx3, signed=True)
    for atom in fo.atoms():
        assert len(set(atom.labels)) == len(atom.labels)
