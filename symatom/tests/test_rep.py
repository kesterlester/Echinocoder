"""Tests for Flavour, FlavouredOperator, repL, repS, PairFlavour."""
import math
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Atom, Context,
    Plan, SimpleCanonicaliser,
)
from symatom.rep import (
    Flavour, FlavouredOperator, repL, repS,
    PairFlavour, pair_flavour_of,
    canonical_pair_flavours, brute_force_canonical_pair_flavours,
)


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
    # (0,3,0) is excluded because Muons has only 2 labels.
    # (0,0,3) is valid: Jets has exactly 3 labels, C(3,3)=1.
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


# ---------------------------------------------------------------------------
# PairFlavour — construction and canonical ordering
# ---------------------------------------------------------------------------

def test_pair_flavour_construction_basic(dot, electrons):
    """PairFlavour builds and stores its fields."""
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(1,))
    assert pf.op_u == dot
    assert pf.flavour_u == fl2
    assert pf.op_v == dot
    assert pf.flavour_v == fl2
    assert pf.overlap == (1,)

def test_pair_flavour_canonical_ordering_swaps_sides(dot, eps3):
    """If (op_v, fl_v) < (op_u, fl_u), the sides are swapped at construction."""
    fl2 = Flavour((2,))
    fl3 = Flavour((3,))
    # "dot" < "eps3" alphabetically, so dot should always end up as op_u.
    pf_natural  = PairFlavour(op_u=dot,  flavour_u=fl2, op_v=eps3, flavour_v=fl3, overlap=(1,))
    pf_reversed = PairFlavour(op_u=eps3, flavour_u=fl3, op_v=dot,  flavour_v=fl2, overlap=(1,))
    assert pf_natural == pf_reversed
    assert pf_natural.op_u == dot
    assert pf_natural.op_v == eps3

def test_pair_flavour_same_op_same_fl_is_symmetric(dot):
    """When both sides are identical, swapping changes nothing."""
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(1,))
    assert pf.op_u == dot and pf.op_v == dot
    assert pf.flavour_u == fl2 and pf.flavour_v == fl2

def test_pair_flavour_frozen(dot):
    """PairFlavour is immutable."""
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(0,))
    with pytest.raises((AttributeError, TypeError)):
        pf.overlap = (1,)  # type: ignore[misc]


# ---------------------------------------------------------------------------
# PairFlavour — validation errors
# ---------------------------------------------------------------------------

def test_pair_flavour_overlap_exceeds_min_flavour(dot):
    """overlap[i] > min(k_u_i, k_v_i) must raise ValueError."""
    fl2 = Flavour((2,))
    fl1 = Flavour((1,))
    with pytest.raises(ValueError, match="overlap"):
        PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl1, overlap=(2,))

def test_pair_flavour_negative_overlap(dot):
    """Negative overlap must raise ValueError."""
    fl2 = Flavour((2,))
    with pytest.raises(ValueError, match="overlap"):
        PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(-1,))

def test_pair_flavour_mismatched_length(dot):
    """overlap length != flavour length must raise ValueError."""
    fl2 = Flavour((2,))
    with pytest.raises(ValueError):
        PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(0, 0))


# ---------------------------------------------------------------------------
# PairFlavour — count()
# ---------------------------------------------------------------------------

def test_pair_flavour_count_dot_dot_single_group(dot):
    """
    With 4 electrons, dot (rank 2) paired with itself at each valid overlap:
      overlap=0: C(4,0)*C(4,2)*C(2,2) = 1*6*1 = 6
      overlap=1: C(4,1)*C(3,1)*C(2,1) = 4*3*2 = 24
      overlap=2: C(4,2)*C(2,0)*C(2,0) = 6*1*1 = 6
    Total = 36 = 6² (6 dot atoms, all ordered pairs).
    """
    fl2 = Flavour((2,))
    group_sizes = (4,)
    counts = {
        s: PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2,
                       overlap=(s,)).count(group_sizes)
        for s in (0, 1, 2)
    }
    assert counts[0] == 6
    assert counts[1] == 24
    assert counts[2] == 6
    assert sum(counts.values()) == 36

def test_pair_flavour_count_two_groups(dot, muons):
    """
    Electrons (n=4, k_u=k_v=2) × Muons (n=2, k_u=k_v=1), overlap (0,0).
    ways_electrons = C(4,0)*C(4,2)*C(2,2) = 6
    ways_muons     = C(2,0)*C(2,1)*C(1,1) = 2
    Total = 12.
    """
    fl = Flavour((2, 1))
    pf = PairFlavour(op_u=dot, flavour_u=fl, op_v=dot, flavour_v=fl, overlap=(0, 0))
    assert pf.count((4, 2)) == 12


# ---------------------------------------------------------------------------
# pair_flavour_of()
# ---------------------------------------------------------------------------

def test_pair_flavour_of_overlap_zero(dot, ctx1):
    """Disjoint atoms: overlap should be (0,)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    pf = pair_flavour_of(u, v, ctx1)
    assert pf.overlap == (0,)
    assert pf.flavour_u == Flavour((2,))
    assert pf.flavour_v == Flavour((2,))

def test_pair_flavour_of_overlap_one(dot, ctx1):
    """Atoms sharing one label: overlap should be (1,)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    pf = pair_flavour_of(u, v, ctx1)
    assert pf.overlap == (1,)

def test_pair_flavour_of_overlap_two_same_atom(dot, ctx1):
    """Same atom paired with itself: overlap equals the flavour count."""
    u = Atom(dot, ("a", "b"), sign=+1)
    pf = pair_flavour_of(u, u, ctx1)
    assert pf.overlap == (2,)

def test_pair_flavour_of_two_groups(dot, eps3, ctx2):
    """Cross-group atom: dot(a, p) has flavour (1,1), self-overlap (1,1)."""
    u = Atom(dot, ("a", "p"), sign=+1)
    pf = pair_flavour_of(u, u, ctx2)
    assert pf.flavour_u == Flavour((1, 1))
    assert pf.overlap == (1, 1)

def test_pair_flavour_of_is_symmetric(dot, ctx1):
    """pair_flavour_of(u, v) == pair_flavour_of(v, u) (canonical ordering)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    assert pair_flavour_of(u, v, ctx1) == pair_flavour_of(v, u, ctx1)

def test_pair_flavour_of_mixed_ops(dot, eps3, ctx1):
    """pair_flavour_of works for atoms with different operations."""
    u = Atom(dot,  ("a", "b"),      sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    pf = pair_flavour_of(u, v, ctx1)
    assert pf.op_u == dot   # dot < eps3 alphabetically
    assert pf.op_v == eps3
    assert pf.overlap == (2,)  # 'a' and 'b' are shared


# ---------------------------------------------------------------------------
# canonical_pair_flavours vs brute_force: cross-validation
# ---------------------------------------------------------------------------

def test_cpf_matches_brute_force_single_group_dot_only(dot, ctx1):
    """With one group and one op, both methods agree."""
    fo_list = repL(ctx1, [dot])
    fast   = set(canonical_pair_flavours(fo_list, ctx1))
    brute  = brute_force_canonical_pair_flavours(fo_list, ctx1)
    assert fast == brute

def test_cpf_matches_brute_force_single_group_two_ops(dot, eps3, ctx1):
    """With one group, dot + eps3: both methods agree."""
    fo_list = repL(ctx1, [dot, eps3])
    fast   = set(canonical_pair_flavours(fo_list, ctx1))
    brute  = brute_force_canonical_pair_flavours(fo_list, ctx1)
    assert fast == brute

def test_cpf_matches_brute_force_two_groups(dot, eps3, ctx2):
    """With two groups (electrons + muons): both methods agree."""
    fo_list = repL(ctx2, [dot, eps3])
    fast   = set(canonical_pair_flavours(fo_list, ctx2))
    brute  = brute_force_canonical_pair_flavours(fo_list, ctx2)
    assert fast == brute

def test_cpf_matches_brute_force_three_groups(mass, dot, eps3, ctx3):
    """With three groups (electrons + muons + jets): both methods agree."""
    fo_list = repL(ctx3, [mass, dot, eps3])
    fast   = set(canonical_pair_flavours(fo_list, ctx3))
    brute  = brute_force_canonical_pair_flavours(fo_list, ctx3)
    assert fast == brute

def test_cpf_empty_fo_list(ctx1):
    """Empty input gives empty output."""
    assert canonical_pair_flavours([], ctx1) == []
    assert brute_force_canonical_pair_flavours([], ctx1) == set()


# ---------------------------------------------------------------------------
# canonical_pair_flavours — specific orbit-structure checks
# ---------------------------------------------------------------------------

def test_cpf_overlap_zero_and_one_are_distinct(dot, ctx1):
    """(dot, dot, overlap=0) and (dot, dot, overlap=1) are different PairFlavours."""
    fl2 = Flavour((2,))
    pf0 = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(0,))
    pf1 = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(1,))
    assert pf0 != pf1
    fo_list = repL(ctx1, [dot])
    result = set(canonical_pair_flavours(fo_list, ctx1))
    assert pf0 in result
    assert pf1 in result

def test_cpf_count_single_group_dot(dot, ctx1):
    """
    4 electrons, dot only.  Valid overlaps for (dot, dot): 0, 1, 2.
    So canonical_pair_flavours should return exactly 3 PairFlavours.
    """
    fo_list = repL(ctx1, [dot])
    result = canonical_pair_flavours(fo_list, ctx1)
    assert len(result) == 3

def test_cpf_count_single_group_dot_eps3(dot, eps3, ctx1):
    """
    4 electrons, dot (rank 2) + eps3 (rank 3).
    dot-dot pairs:   overlap ∈ {0,1,2}          → 3 PairFlavours
    dot-eps3 pairs:  overlap ∈ {max(0,2+3-4)=1, min(2,3)=2} = {1,2} → 2
    eps3-eps3 pairs: overlap ∈ {max(0,3+3-4)=2, min(3,3)=3} = {2,3} → 2
    Total: 7 distinct PairFlavours.
    """
    fo_list = repL(ctx1, [dot, eps3])
    result = canonical_pair_flavours(fo_list, ctx1)
    assert len(result) == 7

def test_cpf_count_satisfies_total_pair_count(dot, ctx1):
    """
    Sum of count() over all PairFlavours must equal (total dot atoms)².
    For 4 electrons with dot: 6 atoms, 36 ordered pairs.
    """
    fo_list = repL(ctx1, [dot])
    group_sizes = (4,)
    total = sum(pf.count(group_sizes) for pf in canonical_pair_flavours(fo_list, ctx1))
    assert total == 36   # 6 dot atoms → 6² ordered pairs

def test_cpf_is_sorted_deterministically(dot, eps3, ctx1):
    """canonical_pair_flavours returns the same sorted list on repeated calls."""
    fo_list = repL(ctx1, [dot, eps3])
    r1 = canonical_pair_flavours(fo_list, ctx1)
    r2 = canonical_pair_flavours(fo_list, ctx1)
    assert r1 == r2


# ---------------------------------------------------------------------------
# SimpleCanonicaliser: atom-tuple order does not affect canonical form
# ---------------------------------------------------------------------------

def test_canon_pair_order_independent(dot, eps3, electrons):
    """
    canon((dot(a,b), eps3(a,b,c))) == canon((eps3(a,b,c), dot(a,b))).
    The SimpleCanonicaliser sorts atoms within the tuple, so input order
    must not affect the canonical form.
    """
    ctx  = Context((electrons,))
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot, eps3))
    u = Atom(dot,  ("a", "b"),      sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    assert plan.canonicalise((u, v)) == plan.canonicalise((v, u))

def test_canon_dot_dot_pair_order_independent(dot, electrons):
    """canon((dot(a,b), dot(c,d))) == canon((dot(c,d), dot(a,b)))."""
    ctx  = Context((electrons,))
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    assert plan.canonicalise((u, v)) == plan.canonicalise((v, u))


# ---------------------------------------------------------------------------
# PairFlavour.orbit_size and orbit_elements
# ---------------------------------------------------------------------------

@pytest.fixture
def ctx4(electrons):
    return Context((electrons,))   # electrons has labels a,b,c,d → size 4

def test_orbit_size_dot_dot_matches_count(dot, ctx4):
    """For SYMMETRIC ops, orbit_size == count (no antisymmetric doubling)."""
    fl2 = Flavour((2,))
    group_sizes = (4,)
    for s in (0, 1, 2):
        pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(s,))
        assert pf.orbit_size(group_sizes) == pf.count(group_sizes)

def test_orbit_size_eps_dot_doubles_for_eps(eps3, dot, ctx4):
    """For one ANTISYMMETRIC op, orbit_size == 2 * count."""
    fl3 = Flavour((3,))
    fl2 = Flavour((2,))
    group_sizes = (4,)
    pf = PairFlavour(op_u=eps3, flavour_u=fl3, op_v=dot, flavour_v=fl2, overlap=(1,))
    assert pf.orbit_size(group_sizes) == 2 * pf.count(group_sizes)

def test_orbit_size_eps_eps_quadruples(eps3, ctx4):
    """For two ANTISYMMETRIC ops, orbit_size == 4 * count."""
    fl3 = Flavour((3,))
    group_sizes = (4,)
    pf = PairFlavour(op_u=eps3, flavour_u=fl3, op_v=eps3, flavour_v=fl3, overlap=(2,))
    assert pf.orbit_size(group_sizes) == 4 * pf.count(group_sizes)

def test_orbit_elements_length_matches_orbit_size(dot, ctx4):
    """orbit_elements returns exactly orbit_size pairs."""
    fl2 = Flavour((2,))
    group_sizes = (4,)
    for s in (0, 1, 2):
        pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(s,))
        elems = pf.orbit_elements(ctx4)
        assert len(elems) == pf.orbit_size(group_sizes)

def test_orbit_elements_length_with_eps(eps3, dot):
    """orbit_elements length matches orbit_size for mixed ANTISYMMETRIC pair."""
    electrons = VectorGroup("electrons", ("a", "b", "c", "d"))
    ctx = Context((electrons,))
    fl3 = Flavour((3,))
    fl2 = Flavour((2,))
    group_sizes = (4,)
    pf = PairFlavour(op_u=eps3, flavour_u=fl3, op_v=dot, flavour_v=fl2, overlap=(1,))
    elems = pf.orbit_elements(ctx)
    assert len(elems) == pf.orbit_size(group_sizes)

def test_orbit_elements_all_have_correct_pair_flavour(dot, ctx4):
    """Every element returned by orbit_elements has the correct PairFlavour."""
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(1,))
    for u, v in pf.orbit_elements(ctx4):
        assert pair_flavour_of(u, v, ctx4) == pf

def test_orbit_elements_sum_over_all_pf_equals_total_pairs(dot, ctx4):
    """
    Sum of orbit_size over all canonical PairFlavours == (total repL atoms)².
    This mirrors the count() consistency check but for orbit_size (which equals
    count() when all ops are SYMMETRIC, as is the case here).
    """
    fo_list = repL(ctx4, [dot])
    group_sizes = tuple(g.size for g in ctx4.groups)
    total_orbit = sum(pf.orbit_size(group_sizes) for pf in canonical_pair_flavours(fo_list, ctx4))
    total_atoms = sum(fo.count() for fo in fo_list)
    assert total_orbit == total_atoms ** 2
