"""
Tests for the canonicalisation contract (C1–C4).

The contract tests are parametrised over both SimpleCanonicaliser (brute force,
permanent reference) and DirectCanonicaliser (fast, active-label algorithm).
Cross-comparison tests assert they produce identical output on every case.
"""
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Atom,
    Context, Plan, SimpleCanonicaliser, DirectCanonicaliser,
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


# ---------------------------------------------------------------------------
# Parametrised over both canonicalisers
# ---------------------------------------------------------------------------

CANONICALISERS = [
    pytest.param(SimpleCanonicaliser(), id="simple"),
    pytest.param(DirectCanonicaliser(), id="direct"),
]

def _plan(canonicaliser, dot, eps3, electrons):
    ctx = Context((electrons,))
    return Plan(context=ctx, canonicaliser=canonicaliser, operations=(dot, eps3))

@pytest.fixture
def plan(dot, eps3, electrons):
    """Default plan fixture uses SimpleCanonicaliser for backward compat."""
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
# C4 — Deterministic:  repeated calls return identical results
# ---------------------------------------------------------------------------

def test_c4_repeated_calls(plan, dot, eps3):
    x = (
        Atom(eps3, ("c", "a", "b"), sign=+1),
        Atom(dot,  ("b", "d"),      sign=+1),
    )
    results = [plan.canonicalise(x) for _ in range(5)]
    assert all(r == results[0] for r in results)


# ---------------------------------------------------------------------------
# Joint canonicalisation of multi-atom tuples
#
# canon((U, V)) is NOT the same as (canon(U), canon(V)).
# The same group element must be applied to all atoms simultaneously.
# Two pairs that differ only in their label-sharing structure are in
# different orbits even when every individual atom canonicalises identically.
# ---------------------------------------------------------------------------

def test_joint_canon_different_orbits_by_label_sharing(plan, dot):
    # Every cross-dot individually canonicalises to the same form (orbit size 1).
    # But the PAIRS differ: (dot(a,b), dot(c,d)) uses four distinct labels
    # while (dot(a,b), dot(a,c)) shares label 'a'.  No relabelling can turn a
    # 4-label pair into a 3-label pair, so they are in different orbits and
    # must have different canonical forms.
    four_label_pair  = (Atom(dot, ("a", "b"), sign=+1), Atom(dot, ("c", "d"), sign=+1))
    three_label_pair = (Atom(dot, ("a", "b"), sign=+1), Atom(dot, ("a", "c"), sign=+1))
    assert plan.canonicalise(four_label_pair) != plan.canonicalise(three_label_pair)

def test_joint_canon_not_independent(plan, dot):
    # Demonstrates the failure mode of naive independent canonicalisation:
    # canon(U) == canon(V) for all cross-dots, so (canon(U), canon(V)) would
    # be the same tuple for both pairs above — collapsing two distinct orbits.
    four_label_pair  = (Atom(dot, ("a", "b"), sign=+1), Atom(dot, ("c", "d"), sign=+1))
    three_label_pair = (Atom(dot, ("a", "b"), sign=+1), Atom(dot, ("a", "c"), sign=+1))
    naive_four  = tuple(plan.canonicalise((a,))[0] for a in four_label_pair)
    naive_three = tuple(plan.canonicalise((a,))[0] for a in three_label_pair)
    # Naive approach wrongly gives the same result for both:
    assert naive_four == naive_three
    # But joint canonicalisation correctly distinguishes them:
    assert plan.canonicalise(four_label_pair) != plan.canonicalise(three_label_pair)


# ---------------------------------------------------------------------------
# Contract tests parametrised over both canonicalisers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("canon", CANONICALISERS)
def test_c1_idempotent_parametrised(canon, dot, eps3, electrons):
    p = _plan(canon, dot, eps3, electrons)
    for x in [
        (Atom(dot,  ("b", "a"),      sign=+1),),
        (Atom(eps3, ("c", "a", "b"), sign=+1),),
        (Atom(eps3, ("b", "c", "a"), sign=-1),),
        (Atom(eps3, ("c", "b", "a"), sign=+1), Atom(dot, ("d", "a"), sign=+1)),
    ]:
        assert p.canonicalise(p.canonicalise(x)) == p.canonicalise(x)

@pytest.mark.parametrize("canon", CANONICALISERS)
def test_c3_same_orbit_parametrised(canon, dot, eps3, electrons):
    p = _plan(canon, dot, eps3, electrons)
    # Related by label permutation b<->c
    x = (Atom(dot, ("a", "b"), sign=+1),)
    y = (Atom(dot, ("a", "c"), sign=+1),)
    assert p.canonicalise(x) == p.canonicalise(y)

@pytest.mark.parametrize("canon", CANONICALISERS)
def test_c3_different_orbit_parametrised(canon, dot, eps3, electrons):
    p = _plan(canon, dot, eps3, electrons)
    x = (Atom(eps3, ("a", "b", "c"), sign=+1),)
    y = (Atom(dot,  ("a", "b"),      sign=+1),)
    assert p.canonicalise(x) != p.canonicalise(y)

@pytest.mark.parametrize("canon", CANONICALISERS)
def test_joint_canon_parametrised(canon, dot, eps3, electrons):
    p = _plan(canon, dot, eps3, electrons)
    four_label_pair  = (Atom(dot, ("a", "b"), sign=+1), Atom(dot, ("c", "d"), sign=+1))
    three_label_pair = (Atom(dot, ("a", "b"), sign=+1), Atom(dot, ("a", "c"), sign=+1))
    assert p.canonicalise(four_label_pair) != p.canonicalise(three_label_pair)


# ---------------------------------------------------------------------------
# Cross-comparison: DirectCanonicaliser must match SimpleCanonicaliser exactly
# ---------------------------------------------------------------------------

def _all_test_cases(dot, eps3):
    return [
        (Atom(dot,  ("b", "a"),      sign=+1),),
        (Atom(dot,  ("d", "c"),      sign=+1),),
        (Atom(eps3, ("c", "a", "b"), sign=+1),),
        (Atom(eps3, ("b", "c", "a"), sign=-1),),
        (Atom(eps3, ("a", "b", "c"), sign=+1),),
        (Atom(eps3, ("a", "b", "c"), sign=-1),),
        (Atom(eps3, ("c", "b", "a"), sign=+1), Atom(dot, ("d", "a"),      sign=+1)),
        (Atom(dot,  ("a", "b"),      sign=+1), Atom(dot, ("c", "d"),      sign=+1)),
        (Atom(dot,  ("a", "b"),      sign=+1), Atom(dot, ("a", "c"),      sign=+1)),
        (Atom(dot,  ("b", "d"),      sign=+1), Atom(eps3, ("a", "c", "d"), sign=+1)),
    ]

def test_direct_matches_simple_single_group(dot, eps3, electrons):
    """DirectCanonicaliser agrees with SimpleCanonicaliser on all single-group cases."""
    simple = SimpleCanonicaliser()
    direct = DirectCanonicaliser()
    ctx = Context((electrons,))
    for x in _all_test_cases(dot, eps3):
        assert direct.canonicalise(x, ctx) == simple.canonicalise(x, ctx), (
            f"Mismatch on {x!r}:\n"
            f"  simple={simple.canonicalise(x, ctx)}\n"
            f"  direct={direct.canonicalise(x, ctx)}"
        )

def test_direct_matches_simple_two_groups(dot, eps3):
    """DirectCanonicaliser agrees with SimpleCanonicaliser in a two-group context."""
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    simple = SimpleCanonicaliser()
    direct = DirectCanonicaliser()
    cases = [
        (Atom(dot,  ("a", "p"), sign=+1),),
        (Atom(dot,  ("b", "q"), sign=+1),),
        (Atom(eps3, ("a", "b", "p"), sign=+1),),
        (Atom(dot,  ("a", "p"), sign=+1), Atom(dot, ("b", "q"), sign=+1)),
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot, ("p", "q"), sign=+1)),
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(dot, ("p", "q"), sign=+1)),
    ]
    for x in cases:
        assert direct.canonicalise(x, ctx) == simple.canonicalise(x, ctx), (
            f"Mismatch on {x!r}:\n"
            f"  simple={simple.canonicalise(x, ctx)}\n"
            f"  direct={direct.canonicalise(x, ctx)}"
        )

# def test_direct_matches_simple_two_groups_of_which_one_is_empty(dot, eps3):
#     """DirectCanonicaliser agrees with SimpleCanonicaliser in a two-group context."""
#     electrons = VectorGroup("electrons", ("a", "b", "c", "d"))
#     muons     = VectorGroup("muons",     tuple())
#     ctx       = Context((electrons, muons))
#     simple = SimpleCanonicaliser()
#     direct = DirectCanonicaliser()
#     cases = [
#         (Atom(dot,  ("a", "b"), sign=+1),),
#         (Atom(dot,  ("d", "b"), sign=+1),),
#         (Atom(dot,  ("c", "a"), sign=+1),),
#         (Atom(eps3, ("a", "b", "c"), sign=+1),),
#         (Atom(eps3, ("a", "d", "c"), sign=+1),),
#         (Atom(eps3, ("a", "d", "c"), sign=-1),),
#         (Atom(eps3, ("b", "d", "c"), sign=-1),),
#     ]
#     for x in cases:
#         assert direct.canonicalise(x, ctx) == simple.canonicalise(x, ctx), (
#             f"Mismatch on {x!r}:\n"
#             f"  simple={simple.canonicalise(x, ctx)}\n"
#             f"  direct={direct.canonicalise(x, ctx)}"
#         )

def test_direct_matches_simple_three_groups_with_funny_names(dot, eps3):
    """DirectCanonicaliser agrees with SimpleCanonicaliser in a two-group context."""
    electrons = VectorGroup("electrons", ("u", "q", "c"))
    muons     = VectorGroup("muons",     ("p", "b"))
    jets      = VectorGroup("jets",      ("a", "v","w","z"))
    ctx       = Context((electrons, muons))
    simple = SimpleCanonicaliser()
    direct = DirectCanonicaliser()
    cases = [
        (Atom(dot,  ("z", "p"), sign=+1),),
        (Atom(dot,  ("a", "p"), sign=+1),),
        (Atom(dot,  ("b", "q"), sign=+1),),
        (Atom(eps3, ("a", "b", "p"), sign=+1),),
        (Atom(eps3, ("z", "b", "p"), sign=+1),),
        (Atom(eps3, ("q", "p", "w"), sign=-1),),
        (Atom(eps3, ("a", "z", "v"), sign=-1),),
        (Atom(dot,  ("a", "p"), sign=+1), Atom(dot, ("b", "q"), sign=+1)),
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot, ("p", "q"), sign=+1)),
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(dot, ("p", "q"), sign=+1)),
    ]
    for x in cases:
        assert direct.canonicalise(x, ctx) == simple.canonicalise(x, ctx), (
            f"Mismatch on {x!r}:\n"
            f"  simple={simple.canonicalise(x, ctx)}\n"
            f"  direct={direct.canonicalise(x, ctx)}"
        )
