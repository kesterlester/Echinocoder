import math
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorType, Atom, Context,
)
from symatom.rep import (
    Flavour, FlavouredOperator, repS,
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
    return VectorType("electrons", ("a", "b", "c",))

@pytest.fixture
def muons():
    return VectorType("muons", ("p", "q"))

@pytest.fixture
def ctx1(electrons):
    """Single-group context: 3 electrons."""
    return Context((electrons,))

@pytest.fixture
def jets():
    return VectorType("jets", ("u", "v", "w"))

@pytest.fixture
def ctx2(electrons, muons):
    """Two-group context: 2 electrons + 2 muons."""
    return Context((electrons, muons))

@pytest.fixture
def ctx3(electrons, muons, jets):
    """Three-group context: 3 electrons + 2 muons + 3 jets."""
    return Context((electrons, muons, jets))


# ---------------------------------------------------------------------------
# PairFlavour — construction and canonical ordering
# ---------------------------------------------------------------------------

# THIS TEST IS CORRECT. DO NOT "FIX" IT IF IT FAILS. FIX THE CODE!!
def test_lester_1(ctx1, eps3, electrons):
    fl2 = Flavour((2,))
    pf = PairFlavour(op_u=eps3, flavour_u=Flavour((3,)), 
                     op_v=eps3, flavour_v=Flavour((3,)), overlap=(3,))

    type_sizes = tuple(g.size for g in ctx1.types)

    print(f"\npf {pf} has atoms:")
    atoms = pf.orbit_elements(ctx1)
    for atom in atoms:
        print("     ",atom)
    os = pf.orbit_size(type_sizes)
    print(f"Orbit size is {os}.")
    print(f"It should have 2 atoms, namely\n    +{atoms[0]} and\n    -{atoms[0]}")
    print("""It should have 2 atoms, namely
      (+eps3(a, b, c), +eps3(a, b, c)) and
      (-eps3(a, b, c), -eps3(a, b, c)).""")

    assert os == 2
    answers =( "(+eps3(a, b, c), +eps3(a, b, c))", "(-eps3(a, b, c), -eps3(a, b, c))" )
    assert str(atoms[0]) in answers
    assert str(atoms[1]) in answers
    assert str(atoms[0]) != str(atoms[1])
    assert atoms[0] != atoms[1]

