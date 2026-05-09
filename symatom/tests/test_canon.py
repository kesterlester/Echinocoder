"""
Tests for the canonicalisation contract (C1–C5) against SimpleCanonicaliser.

The tests are written against the CONTRACT (SPEC.md Section 4.2), not against
any specific canonical form.  A different canonicaliser can be tested by
substituting it into the plan fixture.

NOTE on C4: the sign-consistency property C4 as written in SPEC.md rev 2 is
under review.  For ANTISYMMETRIC operations, atoms (+1, eps, labels) and
(-1, eps, labels) lie in the same orbit (related by an odd label permutation),
so C3 forces canon(x) == canon(y) for both.  This makes the "product of signs"
the same for both canonical forms regardless of the group element connecting
them, which is inconsistent with the stated C4 formula.  C4 will be revised
in a future spec update.  It is not tested here.
"""
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Atom,
    Context, Plan, SimpleCanonicaliser,
)
from symatom.orbits import orbit


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
# C1 — Idempotent:  canon(canon(x)) == canon(x)
# ---------------------------------------------------------------------------

def test_c1_single_dot(plan, dot):
    x = (Atom(dot, ("b", "a"), sign=+1),)
    assert plan.canonicalise(plan.canonicalise(x)) == plan.canonicalise(x)

def test_c1_single_eps_positive(plan, eps3):
    x = (Atom(eps3, ("c", "a", "b"), sign=+1),)
    assert plan.canonicalise(plan.canonicalise(x)) == plan.canonicalise(x)

def test_c1_single_eps_negative(plan, eps3):
    x = (Atom(eps3, ("b", "c", "a"), sign=-1),)
    assert plan.canonicalise(plan.canonicalise(x)) == plan.canonicalise(x)

def test_c1_pair(plan, dot, eps3):
    x = (
        Atom(eps3, ("c", "b", "a"), sign=+1),
        Atom(dot,  ("d", "a"),      sign=+1),
    )
    assert plan.canonicalise(plan.canonicalise(x)) == plan.canonicalise(x)


# ---------------------------------------------------------------------------
# C2 — Representative:  canon(x) lies in the same orbit as x
# (tested as: canon(x) appears in orbit(x, plan))
# ---------------------------------------------------------------------------

def test_c2_single_dot(plan, dot):
    x = (Atom(dot, ("b", "a"), sign=+1),)
    assert plan.canonicalise(x) in orbit(x, plan)

def test_c2_single_eps(plan, eps3):
    x = (Atom(eps3, ("c", "a", "b"), sign=+1),)
    assert plan.canonicalise(x) in orbit(x, plan)

def test_c2_pair(plan, dot, eps3):
    x = (
        Atom(eps3, ("b", "a", "c"), sign=-1),
        Atom(dot,  ("c", "b"),      sign=+1),
    )
    assert plan.canonicalise(x) in orbit(x, plan)


# ---------------------------------------------------------------------------
# C3 — Consistent:  canon(x) == canon(y)  iff  x and y are in the same orbit
# ---------------------------------------------------------------------------

def test_c3_same_orbit_dot_swap(plan, dot):
    # dot(a,b) and dot(b,a) differ only by argument order; dot is symmetric,
    # so arg-sort makes them identical.
    x = (Atom(dot, ("a", "b"), sign=+1),)
    y = (Atom(dot, ("b", "a"), sign=+1),)
    assert plan.canonicalise(x) == plan.canonicalise(y)

def test_c3_same_orbit_label_relabeling(plan, dot):
    # dot(a,b) and dot(a,c) are related by the permutation b<->c.
    x = (Atom(dot, ("a", "b"), sign=+1),)
    y = (Atom(dot, ("a", "c"), sign=+1),)
    assert plan.canonicalise(x) == plan.canonicalise(y)

def test_c3_same_orbit_eps_sign_flip(plan, eps3):
    # eps(a,b,c)+1 and eps(b,a,c)+1 are related by label perm a<->b,
    # which arg-sorts eps(b,a,c) to eps(a,b,c) with sign -1.
    # Both should canonicalise to the same form.
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    y = (Atom(eps3, ("b", "a", "c"), sign=+1),)
    assert plan.canonicalise(x) == plan.canonicalise(y)

def test_c3_same_orbit_eps_positive_and_negative(plan, eps3):
    # eps(a,b,c)+1 and eps(a,b,c)-1 are in the same orbit (via a<->b).
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    y = (Atom(eps3, ("a", "b", "c"), sign=-1),)
    assert plan.canonicalise(x) == plan.canonicalise(y)

def test_c3_different_orbit_different_operation(plan, dot, eps3):
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    y = (Atom(dot,  ("a", "b"),      sign=+1),)
    assert plan.canonicalise(x) != plan.canonicalise(y)

def test_c3_different_orbit_different_label_count(plan, dot, eps3):
    # A 1-atom tuple and a 2-atom tuple are never in the same orbit.
    x = (Atom(dot, ("a", "b"), sign=+1),)
    y = (
        Atom(dot, ("a", "b"), sign=+1),
        Atom(dot, ("a", "c"), sign=+1),
    )
    assert plan.canonicalise(x) != plan.canonicalise(y)


# ---------------------------------------------------------------------------
# C5 — Deterministic:  repeated calls return identical results
# ---------------------------------------------------------------------------

def test_c5_repeated_calls(plan, dot, eps3):
    x = (
        Atom(eps3, ("c", "a", "b"), sign=+1),
        Atom(dot,  ("b", "d"),      sign=+1),
    )
    results = [plan.canonicalise(x) for _ in range(5)]
    assert all(r == results[0] for r in results)
