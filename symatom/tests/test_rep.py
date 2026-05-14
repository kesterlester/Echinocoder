"""Tests for Flavour, FlavouredOperator, repS, PairFlavour."""
import math
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorType, Atom, Context,
)
from symatom.rep import (
    Flavour, FlavouredOperator, repS,
    PairFlavour, pair_flavour_of,
    canonical_pair_flavours, brute_force_canonical_pair_flavours,
    _valid_flavours_rear_first, _valid_flavours_front_first, _valid_flavours,
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
    return VectorType("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def muons():
    return VectorType("muons", ("p", "q"))

@pytest.fixture
def ctx1(electrons):
    """Single-group context: 4 electrons."""
    return Context((electrons,))

@pytest.fixture
def jets():
    return VectorType("jets", ("u", "v", "w"))

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
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1)
    assert fo.operation is dot
    assert fo.flavour == Flavour((2,))

def test_flavoured_operator_wrong_group_count(dot, ctx1):
    with pytest.raises(ValueError, match="vector types"):
        FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx1)

def test_flavoured_operator_wrong_rank(dot, ctx1):
    with pytest.raises(ValueError, match="rank"):
        FlavouredOperator(operation=dot, flavour=Flavour((3,)), context=ctx1)

def test_flavoured_operator_count_exceeds_group(dot):
    # A group of size 1 cannot supply 2 labels for a rank-2 dot.
    tiny_ctx = Context((VectorType("tiny", ("x",)),))
    with pytest.raises(ValueError, match="size"):
        FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=tiny_ctx)


# ---------------------------------------------------------------------------
# FlavouredOperator.count_of_atoms_one_per_sign()
# ---------------------------------------------------------------------------

def test_count_of_atoms_one_per_sign_mass_single_group(mass, ctx1):
    fo = FlavouredOperator(operation=mass, flavour=Flavour((1,)), context=ctx1)
    assert fo.count_of_atoms_one_per_sign() == 4

def test_count_of_atoms_one_per_sign_dot_single_group(dot, ctx1):
    # C(4,2) = 6
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1)
    assert fo.count_of_atoms_one_per_sign() == 6

def test_count_of_atoms_one_per_sign_eps(eps3, ctx1):
    # C(4,3) = 4;
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1)
    assert fo.count_of_atoms_one_per_sign() == 4

def test_count_of_atoms_one_per_sign_dot_two_groups_mixed_flavour(dot, ctx2):
    # C(4,1) * C(2,1) = 4 * 2 = 8; 
    fo = FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx2)
    assert fo.count_of_atoms_one_per_sign() == 8

def test_count_of_atoms_one_per_sign_eps_two_groupsL(eps3, ctx2):
    # flavour (2,1): C(4,2) * C(2,1) = 6 * 2 = 12; 
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((2, 1)), context=ctx2)
    assert fo.count_of_atoms_one_per_sign() == 12


# ---------------------------------------------------------------------------
# FlavouredOperator.atoms_one_per_sign()
# ---------------------------------------------------------------------------

def test_atoms_one_per_sign_mass_single_group(mass, ctx1):
    fo = FlavouredOperator(operation=mass, flavour=Flavour((1,)), context=ctx1)
    atoms = list(fo.atoms_one_per_sign())
    assert len(atoms) == 4
    labels_seen = {a.labels for a in atoms}
    assert labels_seen == {("a",), ("b",), ("c",), ("d",)}
    assert all(a.sign == +1 for a in atoms)

def test_atoms_one_per_sign_dot_single_group(dot, ctx1):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1)
    atoms = list(fo.atoms_one_per_sign())
    assert len(atoms) == 6
    # Every combination of 2 from {a,b,c,d} appears exactly once
    from itertools import combinations
    expected = {tuple(pair) for pair in combinations(("a", "b", "c", "d"), 2)}
    assert {a.labels for a in atoms} == expected

def test_atoms_one_per_sign_eps_repS_single_group(eps3, ctx1):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1)
    atoms = list(fo.atoms_one_per_sign())
    assert len(atoms) == 4
    assert all(a.sign == +1 for a in atoms)

def test_atoms_one_per_sign_count_matches_count_method(dot, eps3, ctx2):
    for op, flavour in [
        (dot,  Flavour((2, 0))),
        (dot,  Flavour((1, 1))),
        (eps3, Flavour((2, 1))),
        (eps3, Flavour((2, 1))),
    ]:
        fo = FlavouredOperator(operation=op, flavour=flavour, context=ctx2)
        assert len(list(fo.atoms_one_per_sign())) == fo.count_of_atoms_one_per_sign()

def test_atoms_one_per_sign_labels_distinct(eps3, ctx2):
    # All generated atoms must satisfy Rule 4 (distinct labels).
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((2, 1)), context=ctx2)
    for atom in fo.atoms_one_per_sign():
        assert len(set(atom.labels)) == len(atom.labels)


# ---------------------------------------------------------------------------
# FlavouredOperator.matches_ignoring_sign()
# ---------------------------------------------------------------------------

def test_matches_ignoring_sign_dot_true(dot, ctx1):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1)
    assert fo.matches_ignoring_sign(Atom(dot, ("a", "b"), sign=+1))

def test_matches_ignoring_sign_dot_wrong_operation(mass, dot, ctx1):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((2,)), context=ctx1)
    assert not fo.matches_ignoring_sign(Atom(mass, ("a",), sign=+1))

def test_matches_ignoring_sign_eps_positive(eps3, ctx1):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1)
    assert fo.matches_ignoring_sign(Atom(eps3, ("a", "b", "c"), sign=+1))

def test_matches_ignoring_sign_eps_negative(eps3, ctx1):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((3,)), context=ctx1)
    assert fo.matches_ignoring_sign(Atom(eps3, ("a", "b", "c"), sign=-1))

def test_matches_ignoring_sign_wrong_flavour(dot, ctx2):
    # dot(a,b) has flavour (2,0); FlavouredOperator has flavour (1,1)
    fo = FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx2)
    assert not fo.matches_ignoring_sign(Atom(dot, ("a", "b"), sign=+1))

def test_matches_ignoring_sign_cross_group_dot(dot, ctx2):
    fo = FlavouredOperator(operation=dot, flavour=Flavour((1, 1)), context=ctx2)
    assert fo.matches_ignoring_sign(Atom(dot, ("a", "p"), sign=+1))
    assert fo.matches_ignoring_sign(Atom(dot, ("b", "q"), sign=+1))  # label order doesn't matter


# ---------------------------------------------------------------------------
# repS function
# ---------------------------------------------------------------------------

def test_repS_single_group_single_op(eps3, ctx1):
    # eps3 rank 3, 4 electrons: only flavour (3,) is valid → one FlavouredOperator
    fos = repS(ctx1, [eps3])
    assert len(fos) == 1
    assert fos[0].flavour == Flavour((3,))

def test_repS_single_group_single_op(eps3, ctx1):
    fos = repS(ctx1, [eps3])
    assert len(fos) == 1

def test_repS_two_groups_eps3(eps3, ctx2):
    fos = repS(ctx2, [eps3])
    flavours = {fo.flavour for fo in fos}
    assert flavours == {Flavour((3, 0)), Flavour((2, 1)), Flavour((1, 2))}

def test_repS_two_groups_multiple_ops(mass, dot, eps3, ctx2):
    fos = repS(ctx2, [mass, dot, eps3])
    # mass rank 1: flavours (1,0), (0,1) → 2
    # dot  rank 2: flavours (2,0), (1,1), (0,2) → 3
    # eps3 rank 3: flavours (3,0), (2,1), (1,2) → 3
    assert len(fos) == 8

def test_repS_total_atom_count(eps3, ctx1):
    # repS with one eps3 and 4 electrons: C(4,3) = 4 atoms total
    fos = repS(ctx1, [eps3])
    total = sum(fo.count_of_atoms_one_per_sign() for fo in fos)
    assert total == 4

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
    assert len(repS(ctx3, [mass, dot, eps3])) == 18

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
    mass_total = sum(fo.count_of_atoms_one_per_sign() for fo in repS(ctx3, [mass]))
    dot_total  = sum(fo.count_of_atoms_one_per_sign() for fo in repS(ctx3, [dot]))
    eps3_total = sum(fo.count_of_atoms_one_per_sign() for fo in repS(ctx3, [eps3]))
    assert mass_total == 9
    assert dot_total  == 36
    assert eps3_total == 84

def test_three_group_grand_total_repS(mass, dot, eps3, ctx3):
    assert sum(fo.count_of_atoms_one_per_sign() for fo in repS(ctx3, [mass, dot, eps3])) == 129

def test_three_group_atoms_matches_count(mass, dot, eps3, ctx3):
    # For every FlavouredOperator, atoms() must yield exactly count() atoms.
    for fo in repS(ctx3, [mass, dot, eps3]):
        assert len(list(fo.atoms_one_per_sign())) == fo.count_of_atoms_one_per_sign()

def test_three_group_eps3_mixed_flavour_spot_check(eps3, ctx3):
    # Flavour (1,1,1): C(4,1)*C(2,1)*C(3,1) = 4*2*3 = 24 atoms in repS, 48 in repL.
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((1,1,1)), context=ctx3)
    assert fo.count_of_atoms_one_per_sign() == 24

def test_three_group_atoms_labels_distinct(eps3, ctx3):
    fo = FlavouredOperator(operation=eps3, flavour=Flavour((1,1,1)), context=ctx3)
    for atom in fo.atoms_one_per_sign():
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
    type_sizes = (4,)
    counts = {
        s: PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2,
                       overlap=(s,)).count(type_sizes)
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
    fo_list = repS(ctx1, [dot])
    fast   = set(canonical_pair_flavours(fo_list, ctx1))
    brute  = brute_force_canonical_pair_flavours(fo_list, ctx1)
    assert fast == brute

def test_cpf_matches_brute_force_single_group_two_ops(dot, eps3, ctx1):
    """With one group, dot + eps3: both methods agree."""
    fo_list = repS(ctx1, [dot, eps3])
    fast   = set(canonical_pair_flavours(fo_list, ctx1))
    brute  = brute_force_canonical_pair_flavours(fo_list, ctx1)
    assert fast == brute

def test_cpf_matches_brute_force_two_groups(dot, eps3, ctx2):
    """With two groups (electrons + muons): both methods agree."""
    fo_list = repS(ctx2, [dot, eps3])
    fast   = set(canonical_pair_flavours(fo_list, ctx2))
    brute  = brute_force_canonical_pair_flavours(fo_list, ctx2)
    assert fast == brute

def test_cpf_matches_brute_force_three_groups(mass, dot, eps3, ctx3):
    """With three groups (electrons + muons + jets): both methods agree."""
    fo_list = repS(ctx3, [mass, dot, eps3])
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
    fo_list = repS(ctx1, [dot])
    result = set(canonical_pair_flavours(fo_list, ctx1))
    assert pf0 in result
    assert pf1 in result

def test_cpf_count_single_group_dot(dot, ctx1):
    """
    4 electrons, dot only.  Valid overlaps for (dot, dot): 0, 1, 2.
    So canonical_pair_flavours should return exactly 3 PairFlavours.
    """
    fo_list = repS(ctx1, [dot])
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
    fo_list = repS(ctx1, [dot, eps3])
    result = canonical_pair_flavours(fo_list, ctx1)
    assert len(result) == 7

def test_cpf_count_satisfies_total_pair_count(dot, ctx1):
    """
    Sum of count() over all PairFlavours must equal (total dot atoms)².
    For 4 electrons with dot: 6 atoms, 36 ordered pairs.
    """
    fo_list = repS(ctx1, [dot])
    type_sizes = (4,)
    total = sum(pf.count(type_sizes) for pf in canonical_pair_flavours(fo_list, ctx1))
    assert total == 36   # 6 dot atoms → 6² ordered pairs

def test_cpf_is_sorted_deterministically(dot, eps3, ctx1):
    """canonical_pair_flavours returns the same sorted list on repeated calls."""
    fo_list = repS(ctx1, [dot, eps3])
    r1 = canonical_pair_flavours(fo_list, ctx1)
    r2 = canonical_pair_flavours(fo_list, ctx1)
    assert r1 == r2


# ---------------------------------------------------------------------------
# PairFlavour.orbit_size and orbit_elements
# ---------------------------------------------------------------------------

@pytest.fixture
def ctx4(electrons):
    return Context((electrons,))   # electrons has labels a,b,c,d → size 4

def test_orbit_size_dot_dot_matches_count(dot, ctx4):
    """For SYMMETRIC ops, orbit_size == count (no antisymmetric doubling)."""
    fl2 = Flavour((2,))
    type_sizes = (4,)
    for s in (0, 1, 2):
        pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(s,))
        assert pf.orbit_size(type_sizes) == pf.count(type_sizes)

def test_orbit_size_eps_dot_doubles_for_eps(eps3, dot, ctx4):
    """For one ANTISYMMETRIC op, orbit_size == 2 * count."""
    fl3 = Flavour((3,))
    fl2 = Flavour((2,))
    type_sizes = (4,)
    pf = PairFlavour(op_u=eps3, flavour_u=fl3, op_v=dot, flavour_v=fl2, overlap=(1,))
    assert pf.orbit_size(type_sizes) == 2 * pf.count(type_sizes)

def test_orbit_size_eps_eps_nonzero_overlap(eps3, ctx4):
    """For two ANTISYMMETRIC ops with non-zero overlap, orbit_size == 2 * count.

    Non-zero overlap couples the sign flips on both sides: any permutation that
    flips the sign of u (by creating an odd permutation on u's labels) must also
    flip the sign of v (the shared labels contribute the same odd-even parity to
    both).  The orbit therefore contains only (+u,+v) and (-u,-v) variants, not
    all four sign combinations.  orbit_size = 2 * count, not 4 * count.

    4 * count arises only for zero-overlap AA pairs with rank >= 2 on both sides;
    that requires n >= ku+kv which is impossible here with ku=kv=3 and n=4.
    """
    fl3 = Flavour((3,))
    type_sizes = (4,)
    pf = PairFlavour(op_u=eps3, flavour_u=fl3, op_v=eps3, flavour_v=fl3, overlap=(2,))
    assert pf.orbit_size(type_sizes) == 2 * pf.count(type_sizes)

def test_orbit_elements_length_matches_orbit_size(dot, ctx4):
    """orbit_elements returns exactly orbit_size pairs."""
    fl2 = Flavour((2,))
    type_sizes = (4,)
    for s in (0, 1, 2):
        pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(s,))
        elems = pf.orbit_elements(ctx4)
        assert len(elems) == pf.orbit_size(type_sizes)

def test_orbit_elements_length_with_eps(eps3, dot):
    """orbit_elements length matches orbit_size for mixed ANTISYMMETRIC pair."""
    electrons = VectorType("electrons", ("a", "b", "c", "d"))
    ctx = Context((electrons,))
    fl3 = Flavour((3,))
    fl2 = Flavour((2,))
    type_sizes = (4,)
    pf = PairFlavour(op_u=eps3, flavour_u=fl3, op_v=dot, flavour_v=fl2, overlap=(1,))
    elems = pf.orbit_elements(ctx)
    assert len(elems) == pf.orbit_size(type_sizes)

def test_orbit_elements_all_have_correct_pair_flavour(dot, ctx4):
    """Every element returned by orbit_elements has the correct PairFlavour."""
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(1,))
    for u, v in pf.orbit_elements(ctx4):
        assert pair_flavour_of(u, v, ctx4) == pf

def test_orbit_elements_sum_over_all_pf_equals_total_pairs(dot, ctx4):
    """
    Sum of orbit_size over all canonical PairFlavours == (total repS atoms)².
    This mirrors the count() consistency check but for orbit_size (which equals
    count() when all ops are SYMMETRIC, as is the case here).
    """
    fo_list = repS(ctx4, [dot])
    type_sizes = tuple(g.size for g in ctx4.types)
    total_orbit = sum(pf.orbit_size(type_sizes) for pf in canonical_pair_flavours(fo_list, ctx4))
    total_atoms = sum(fo.count_of_atoms_one_per_sign() for fo in fo_list)
    assert total_orbit == total_atoms ** 2


# ---------------------------------------------------------------------------
# canonical_pair_flavours: self-pairing and block-consecutiveness guarantees
# ---------------------------------------------------------------------------

def test_cpf_includes_self_pairing(dot, eps3, ctx2):
    """
    canonical_pair_flavours includes self-pairings (fo_u == fo_v).
    This is intentional: pairing a FO with itself captures correlations between
    distinct atoms of the same type, and full-overlap self-pairs degenerate to
    the single-atom orbit (handled by Phase 1 encoding under 5b).
    """
    fo_list = repS(ctx2, [dot, eps3])
    pf_list = canonical_pair_flavours(fo_list, ctx2)
    # At least one PairFlavour must pair an operation with itself
    self_pairs = [pf for pf in pf_list if pf.op_u == pf.op_v and pf.flavour_u == pf.flavour_v]
    assert len(self_pairs) > 0, "Expected self-pairings (fo_u == fo_v) to be present"

def test_cpf_every_fo_self_paired(dot, eps3, ctx2):
    """Every FlavouredOperator in repS is paired with itself at some overlap."""
    fo_list = repS(ctx2, [dot, eps3])
    pf_list = canonical_pair_flavours(fo_list, ctx2)
    for fo in fo_list:
        found = any(
            pf.op_u == fo.operation and pf.flavour_u == fo.flavour and
            pf.op_v == fo.operation and pf.flavour_v == fo.flavour
            for pf in pf_list
        )
        assert found, f"No self-pairing found for {fo!r}"

def test_cpf_overlap_blocks_are_contiguous(dot, eps3, ctx2):
    """
    PairFlavours with the same (op_u, flavour_u, op_v, flavour_v) — an OVERLAP
    BLOCK — must appear as a contiguous run in the canonical_pair_flavours output.
    encode() and describe_encoding() rely on this via itertools.groupby.
    This test will catch any future change to the sort order that breaks it.
    """
    fo_list = repS(ctx2, [dot, eps3])
    pf_list = canonical_pair_flavours(fo_list, ctx2)

    def block_key(pf):
        return (pf.op_u.name, pf.flavour_u.counts, pf.op_v.name, pf.flavour_v.counts)

    seen_blocks = {}
    for idx, pf in enumerate(pf_list):
        key = block_key(pf)
        if key not in seen_blocks:
            seen_blocks[key] = idx
        else:
            # If we've seen this block key before, the previous occurrence must
            # be the immediately preceding index (no gap allowed).
            prev_idx = seen_blocks[key]
            assert prev_idx == idx - 1, (
                f"OVERLAP BLOCK {key} is not contiguous: "
                f"last seen at index {prev_idx}, reappears at index {idx}"
            )
            seen_blocks[key] = idx

def test_cpf_overlap_blocks_contiguous_three_groups(dot, eps3, ctx3):
    """Block-consecutiveness holds in a three-group context too."""
    fo_list = repS(ctx3, [dot, eps3])
    pf_list = canonical_pair_flavours(fo_list, ctx3)

    def block_key(pf):
        return (pf.op_u.name, pf.flavour_u.counts, pf.op_v.name, pf.flavour_v.counts)

    seen_blocks = {}
    for idx, pf in enumerate(pf_list):
        key = block_key(pf)
        if key in seen_blocks:
            prev_idx = seen_blocks[key]
            assert prev_idx == idx - 1, (
                f"OVERLAP BLOCK {key} is not contiguous in three-group context: "
                f"last seen at index {prev_idx}, reappears at index {idx}"
            )
        seen_blocks[key] = idx


# ---------------------------------------------------------------------------
# _valid_flavours_rear_first and _valid_flavours_front_first
# ---------------------------------------------------------------------------

# Shared helper: build VectorType objects from (name, labels) pairs.
def _types(*specs):
    return tuple(VectorType(name, labels) for name, labels in specs)


def test_valid_flavours_unconstrained_rear_first():
    """
    With three groups each of size ≥ rank, rear_first iterates k=0..max at each
    level, so later types fill first.
    Types (2,2,2), rank 2 → expected ascending-lex order.
    """
    types = _types(("e", ("a","b")), ("m", ("c","d")), ("j", ("u","v")))
    result = list(_valid_flavours_rear_first(types, 2))
    expected = [
        Flavour((0,0,2)), Flavour((0,1,1)), Flavour((0,2,0)),
        Flavour((1,0,1)), Flavour((1,1,0)), Flavour((2,0,0)),
    ]
    assert result == expected


def test_valid_flavours_unconstrained_front_first():
    """
    With three groups each of size ≥ rank, front_first iterates k=max..0 at each
    level, so earlier types fill first.
    Types (2,2,2), rank 2 → expected descending-lex order.
    """
    types = _types(("e", ("a","b")), ("m", ("c","d")), ("j", ("u","v")))
    result = list(_valid_flavours_front_first(types, 2))
    expected = [
        Flavour((2,0,0)), Flavour((1,1,0)), Flavour((1,0,1)),
        Flavour((0,2,0)), Flavour((0,1,1)), Flavour((0,0,2)),
    ]
    assert result == expected


def test_valid_flavours_front_is_exact_reverse_of_rear_unconstrained():
    """front_first is always the exact reverse of rear_first (unconstrained sizes)."""
    types = _types(("e", ("a","b")), ("m", ("c","d")), ("j", ("u","v")))
    rear  = list(_valid_flavours_rear_first(types, 2))
    front = list(_valid_flavours_front_first(types, 2))
    assert front == list(reversed(rear))


def test_valid_flavours_constrained_sizes():
    """
    With types of sizes (1, 2, 2) and rank 2, the flavour (2,*,*) is excluded.
    rear_first:  (0,0,2), (0,1,1), (0,2,0), (1,0,1), (1,1,0)
    front_first: (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)
    """
    types = _types(("e", ("a",)), ("m", ("c","d")), ("j", ("u","v")))
    rear  = list(_valid_flavours_rear_first(types, 2))
    front = list(_valid_flavours_front_first(types, 2))
    expected_rear = [
        Flavour((0,0,2)), Flavour((0,1,1)), Flavour((0,2,0)),
        Flavour((1,0,1)), Flavour((1,1,0)),
    ]
    assert rear  == expected_rear
    assert front == list(reversed(expected_rear))
    # (2,*,*) must never appear — size of first group is 1
    assert all(f.counts[0] <= 1 for f in rear)
    assert all(f.counts[0] <= 1 for f in front)


def test_valid_flavours_front_is_exact_reverse_of_rear_constrained():
    """front_first is the exact reverse of rear_first even when sizes are constrained."""
    types = _types(("e", ("a",)), ("m", ("c","d")), ("j", ("u","v")))
    rear  = list(_valid_flavours_rear_first(types, 2))
    front = list(_valid_flavours_front_first(types, 2))
    assert front == list(reversed(rear))


def test_valid_flavours_same_set_both_orders():
    """Both variants enumerate exactly the same set of Flavours."""
    types = _types(("e", ("a","b","c","d")), ("m", ("p","q")), ("j", ("u","v","w")))
    for rank in (1, 2, 3):
        rear  = set(_valid_flavours_rear_first(types, rank))
        front = set(_valid_flavours_front_first(types, rank))
        assert rear == front, f"Sets differ at rank {rank}"


def test_valid_flavours_rank_zero():
    """Rank 0 yields exactly one Flavour: all zeros."""
    types = _types(("e", ("a","b")), ("m", ("c","d")))
    assert list(_valid_flavours_rear_first(types, 0))  == [Flavour((0, 0))]
    assert list(_valid_flavours_front_first(types, 0)) == [Flavour((0, 0))]


def test_valid_flavours_alias_is_front_first():
    """The _valid_flavours alias must point to the front-first variant."""
    types = _types(("e", ("a","b")), ("m", ("c","d")), ("j", ("u","v")))
    assert (list(_valid_flavours(types, 2))
            == list(_valid_flavours_front_first(types, 2)))
