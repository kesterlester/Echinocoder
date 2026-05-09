"""Tests for orbit(), stabiliser_size(), orbit_and_stabiliser_size()."""
import math
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Atom,
    Context, Plan, SimpleCanonicaliser,
    orbit, stabiliser_size, orbit_and_stabiliser_size,
)


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
def electrons():
    return VectorGroup("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def plan(dot, eps3, electrons):
    ctx = Context((electrons,))
    return Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot, eps3))


# ---------------------------------------------------------------------------
# Orbit membership and structure
# ---------------------------------------------------------------------------

def test_orbit_contains_canonical_of_input(plan, dot):
    x = (Atom(dot, ("b", "a"), sign=+1),)
    orb = orbit(x, plan)
    assert plan.canonicalise(x) in orb

def test_orbit_all_elements_are_canonical(plan, eps3):
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    for elem in orbit(x, plan):
        assert plan.canonicalise(elem) == elem

def test_orbit_no_duplicates(plan, dot):
    x = (Atom(dot, ("a", "b"), sign=+1),)
    orb = orbit(x, plan)
    assert len(orb) == len(set(orb))

def test_orbit_same_for_all_representatives(plan, dot):
    # dot(a,b) and dot(b,a) are in the same orbit; orbit() should agree.
    x = (Atom(dot, ("a", "b"), sign=+1),)
    y = (Atom(dot, ("b", "a"), sign=+1),)
    assert set(orbit(x, plan)) == set(orbit(y, plan))


# ---------------------------------------------------------------------------
# Orbit sizes and the orbit-stabiliser theorem
# ---------------------------------------------------------------------------

def test_orbit_stabiliser_theorem(plan, dot, electrons):
    """For every test atom tuple: |orbit| * |stabiliser| == |G|."""
    G = math.factorial(electrons.size)
    atoms_to_test = [
        (Atom(dot, ("a", "b"), sign=+1),),
        (Atom(dot, ("a", "a"), sign=+1),),   # self-dot: smaller orbit
    ]
    for x in atoms_to_test:
        orb, stab = orbit_and_stabiliser_size(x, plan)
        assert len(orb) * stab == G

def test_orbit_stabiliser_theorem_eps(plan, eps3, electrons):
    G = math.factorial(electrons.size)
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    orb, stab = orbit_and_stabiliser_size(x, plan)
    assert len(orb) * stab == G

def test_orbit_size_dot_ab(plan, dot, electrons):
    # dot(a,b) with 4 electrons: orbit contains dot(x,y) for all
    # unordered pairs {x,y} with x <= y, plus self-dots with x==y.
    # Unordered pairs from 4 labels: C(4,2) = 6, plus 4 self-dots = 10.
    # But dot(a,b) ~ dot(a,c) ~ dot(b,d) etc all collapse to the same
    # canonical form — so orbit size should be 1 (all cross-dots are equivalent).
    x = (Atom(dot, ("a", "b"), sign=+1),)
    assert len(orbit(x, plan)) == 1

def test_orbit_size_self_dot(plan, dot):
    # dot(a,a) is a self-dot; it's in a different orbit from dot(a,b).
    # All self-dots are equivalent under the group, so orbit size == 1.
    x = (Atom(dot, ("a", "a"), sign=+1),)
    assert len(orbit(x, plan)) == 1

def test_self_dot_and_cross_dot_different_orbits(plan, dot):
    x = (Atom(dot, ("a", "a"), sign=+1),)
    y = (Atom(dot, ("a", "b"), sign=+1),)
    assert set(orbit(x, plan)) != set(orbit(y, plan))

def test_orbit_size_eps_single_atom(plan, eps3, electrons):
    # eps3(a,b,c) with 4 electrons: orbit size == 1.
    # S4 acts transitively on 3-element subsets of {a,b,c,d}: the permutation
    # (c<->d) maps eps3(a,b,c) to eps3(a,b,d), so all four subsets are in the
    # same orbit.  After canonicalisation they all collapse to a single form.
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    assert len(orbit(x, plan)) == 1

def test_orbit_size_eps_stabiliser(plan, eps3, electrons):
    # orbit size = 1, so stabiliser = |G| / 1 = 24.
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    assert stabiliser_size(x, plan) == 24


# ---------------------------------------------------------------------------
# Two-group context
# ---------------------------------------------------------------------------

def test_two_group_orbit(dot):
    electrons = VectorGroup("electrons", ("a", "b"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx  = Context((electrons, muons))
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))

    # dot(a, p) — one electron label, one muon label.
    # Full group = S2 x S2, size 4.
    # Orbit: dot(a,p), dot(a,q), dot(b,p), dot(b,q) all collapse to same
    # canonical form, so orbit size == 1.
    x = (Atom(dot, ("a", "p"), sign=+1),)
    assert len(orbit(x, plan)) == 1

def test_two_group_orbit_stabiliser_theorem(dot):
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx  = Context((electrons, muons))
    plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(dot,))
    G    = math.factorial(3) * math.factorial(2)  # 6 * 2 = 12

    x = (Atom(dot, ("a", "p"), sign=+1),)
    orb, stab = orbit_and_stabiliser_size(x, plan)
    assert len(orb) * stab == G
