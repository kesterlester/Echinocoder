"""
Tests for symatom.group: TheGroup.

TheGroup wraps a Context and provides orbit enumeration, stabiliser size,
and orbit-membership testing.

All methods are currently brute-force O(∏n_g!).  These tests verify
correctness (not performance) and serve as the permanent reference against
which faster future implementations will be validated.

Note: sign_correlation_type() and SignCorrelationType have moved to
symcoder/encode.py and are no longer part of symatom.
"""
import math
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorType, Atom, Context,
    TheGroup,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def dot():
    return Operation("dot", rank=2, odd_parity=False, argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                   eval_fn=lambda v: 0.0)

@pytest.fixture
def eps3():
    return Operation("eps3", rank=3, odd_parity=True, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                   eval_fn=lambda v: 0.0)

@pytest.fixture
def eps2():
    """Rank-2 antisymmetric operation — useful for small-n overlap tests."""
    return Operation("eps2", rank=2, odd_parity=True, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                   eval_fn=lambda v: 0.0)

@pytest.fixture
def electrons4():
    return VectorType("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def electrons3():
    return VectorType("electrons", ("a", "b", "c"))

@pytest.fixture
def ctx4(electrons4):
    return Context((electrons4,))

@pytest.fixture
def ctx3(electrons3):
    return Context((electrons3,))

@pytest.fixture
def ctx_two_groups():
    e = VectorType("electrons", ("a", "b", "c"))
    mu = VectorType("muons", ("p", "q"))
    return Context((e, mu))


# ---------------------------------------------------------------------------
# order()
# ---------------------------------------------------------------------------

def test_order_single_group_n4(ctx4):
    g = ctx4.the_group
    assert g.order() == math.factorial(4)  # 24

def test_order_single_group_n3(ctx3):
    g = ctx3.the_group
    assert g.order() == math.factorial(3)  # 6

def test_order_two_groups(ctx_two_groups):
    g = ctx_two_groups.the_group
    # 3! * 2! = 12
    assert g.order() == math.factorial(3) * math.factorial(2)

def test_order_empty_context():
    ctx = Context((VectorType("e", ()),))
    g = ctx.the_group
    assert g.order() == 1  # 0! = 1


# ---------------------------------------------------------------------------
# orbit() — basic properties
# ---------------------------------------------------------------------------

def test_orbit_contains_representative(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = ctx4.the_group
    orb = g.orbit_pair(u, v)
    assert (u, v) in orb

def test_orbit_pair_first_element_is_representative(dot, ctx4):
    """orbit_pair() first element is the input pair (from identity permutation)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    g = ctx4.the_group
    orb = g.orbit_pair(u, v)
    assert orb[0] == (u, v)

def test_orbit_pair_no_duplicates(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = ctx4.the_group
    orb = g.orbit_pair(u, v)
    assert len(orb) == len(set(orb))

def test_orbit_pair_length_matches_orbit_size_pair(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = ctx4.the_group
    assert len(g.orbit_pair(u, v)) == g.orbit_size_pair(u, v)

def test_orbit_pair_all_pairs_are_atom_tuples(dot, eps3, ctx4):
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    for pair in g.orbit_pair(u, v):
        assert isinstance(pair, tuple) and len(pair) == 2
        assert isinstance(pair[0], Atom) and isinstance(pair[1], Atom)


# ---------------------------------------------------------------------------
# stabiliser_size() and orbit_size()
# ---------------------------------------------------------------------------

def test_orbit_size_pair_times_stab_equals_order(dot, eps3, ctx4):
    """orbit_size_pair * stabiliser_size_pair == order (orbit-stabiliser theorem)."""
    g = ctx4.the_group
    pairs = [
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot,  ("c", "d"), sign=+1)),
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot,  ("a", "c"), sign=+1)),
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot,  ("a", "b"), sign=+1)),
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(eps3, ("a", "b", "c"), sign=+1)),
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(dot, ("a", "b"), sign=+1)),
    ]
    for u, v in pairs:
        assert g.orbit_size_pair(u, v) * g.stabiliser_size_pair(u, v) == g.order(), (
            f"Orbit-stabiliser theorem violated for {u!r}, {v!r}"
        )

def test_orbit_size_pair_ss_zero_overlap(dot, ctx4):
    """SS no-overlap: C(4,2)*C(2,2) = 6."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = ctx4.the_group
    assert g.orbit_size_pair(u, v) == 6

def test_orbit_size_pair_ss_partial_overlap(dot, ctx4):
    """SS overlap=1: C(4,1)*C(3,1)*C(2,1) = 24? No: count=C(4,1)*C(3,1)*C(2,1)=...

    For SS s=1, a=1, b=1, r=1, n=4:
      stab = 1!*1!*1!*1! = 1
      orbit_size = 4!/1 = 24
    But count = C(4,1)*C(3,1)*C(2,1) = 4*3*2 = 24 — SS count == orbit_size. Correct.
    """
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    g = ctx4.the_group
    assert g.orbit_size_pair(u, v) == 24

def test_orbit_size_pair_ss_full_overlap(dot, ctx4):
    """SS full overlap (u==v): stab = 2!*0!*0!*2! = 4, orbit_size = 24/4 = 6."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    assert g.orbit_size_pair(u, v) == 6

def test_orbit_size_pair_aa_full_overlap(eps3, ctx3):
    """AA full overlap: stab = E(3) = 3, orbit_size = 6/3 = 2."""
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = ctx3.the_group
    assert g.orbit_size_pair(u, v) == 2

def test_orbit_size_pair_aa_zero_overlap(eps2, ctx4):
    """AA zero overlap, rank=2, n=4: stab=1, orbit_size=24."""
    u = Atom(eps2, ("a", "b"), sign=+1)
    v = Atom(eps2, ("c", "d"), sign=+1)
    g = ctx4.the_group
    assert g.orbit_size_pair(u, v) == 24

def test_orbit_size_pair_matches_len_orbit_pair(eps3, ctx3):
    """orbit_size_pair computed algebraically must equal len(orbit_pair()) for all test cases."""
    g = ctx3.the_group
    test_pairs = [
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(eps3, ("a", "b", "c"), sign=+1)),
    ]
    for u, v in test_pairs:
        assert g.orbit_size_pair(u, v) == len(g.orbit_pair(u, v))


# ---------------------------------------------------------------------------
# orbit() agrees with pf.orbit_elements() for canonical pairs
# ---------------------------------------------------------------------------

def test_orbit_pair_agrees_with_pf_orbit_elements(dot, eps3):
    """
    For the canonical representative of every PairFlavour, TheGroup.orbit_pair and
    pf.orbit_elements must return the same set of pairs.
    """
    from symatom import repS, canonical_pair_flavours
    electrons = VectorType("electrons", ("a", "b", "c", "d"))
    ctx = Context((electrons,))
    g = ctx.the_group
    for pf in canonical_pair_flavours(repS(ctx, [dot, eps3]), ctx):
        # canonical representative: first labels of each group up to flavour count
        ku = pf.flavour_u.counts[0]
        kv = pf.flavour_v.counts[0]
        s  = pf.overlap[0]
        u_labels = electrons.labels[:ku]
        v_labels = electrons.labels[:s] + electrons.labels[ku:ku + kv - s]
        u = Atom(pf.op_u, u_labels, sign=+1)
        v = Atom(pf.op_v, v_labels, sign=+1)
        assert set(g.orbit_pair(u, v)) == set(pf.orbit_elements(ctx)), (
            f"orbit_pair disagrees with pf.orbit_elements for {pf!r}"
        )


# ---------------------------------------------------------------------------
# in_orbit()
# ---------------------------------------------------------------------------

def test_in_orbit_pair_self(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = ctx4.the_group
    assert g.in_orbit_pair((u, v), (u, v))

def test_in_orbit_pair_orbit_element(dot, ctx4):
    """Every element of orbit_pair(u,v) is in_orbit_pair((u,v), (u,v))."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    g = ctx4.the_group
    for pair in g.orbit_pair(u, v):
        assert g.in_orbit_pair(pair, (u, v))

def test_not_in_orbit_pair_different_flavour(dot, ctx4):
    """A pair with different label counts is not in the orbit."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    # Different PairFlavour: overlap=1 instead of 0
    u2 = Atom(dot, ("a", "b"), sign=+1)
    v2 = Atom(dot, ("a", "c"), sign=+1)
    g = ctx4.the_group
    assert not g.in_orbit_pair((u2, v2), (u, v))


# ---------------------------------------------------------------------------
# Two-group context
# ---------------------------------------------------------------------------

def test_orbit_size_pair_two_groups(dot, ctx_two_groups):
    """orbit_size_pair * stabiliser_size_pair == order in a two-group context."""
    g = ctx_two_groups.the_group
    u = Atom(dot, ("a", "p"), sign=+1)
    v = Atom(dot, ("b", "q"), sign=+1)
    assert g.orbit_size_pair(u, v) * g.stabiliser_size_pair(u, v) == g.order()

def test_orbit_pair_length_two_groups(dot, ctx_two_groups):
    """len(orbit_pair) == orbit_size_pair in a two-group context."""
    g = ctx_two_groups.the_group
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("p", "q"), sign=+1)
    assert len(g.orbit_pair(u, v)) == g.orbit_size_pair(u, v)


# ---------------------------------------------------------------------------
# Step 3 cross-validation: primary methods must agree with *_brute variants
# ---------------------------------------------------------------------------

@pytest.fixture
def pairs_for_xval(dot, eps3, eps2, ctx4):
    """Representative atom-pairs for cross-validation sweeps."""
    return [
        # SS: dot × dot
        (Atom(dot,  ("a", "b"),      +1), Atom(dot,  ("c", "d"),      +1)),
        (Atom(dot,  ("a", "b"),      +1), Atom(dot,  ("a", "c"),      +1)),
        (Atom(dot,  ("a", "b"),      +1), Atom(dot,  ("a", "b"),      +1)),
        # SA: dot × eps3
        (Atom(dot,  ("a", "b"),      +1), Atom(eps3, ("a", "b", "c"), +1)),
        (Atom(dot,  ("a", "b"),      +1), Atom(eps3, ("b", "c", "d"), +1)),
        # AS: eps3 × dot
        (Atom(eps3, ("a", "b", "c"), +1), Atom(dot,  ("a", "b"),      +1)),
        # AA full overlap
        (Atom(eps3, ("a", "b", "c"), +1), Atom(eps3, ("a", "b", "c"), +1)),
        (Atom(eps2, ("a", "b"),      +1), Atom(eps2, ("a", "b"),      +1)),
        # AA zero overlap
        (Atom(eps2, ("a", "b"),      +1), Atom(eps2, ("c", "d"),      +1)),
        # AA partial overlap
        (Atom(eps2, ("a", "b"),      +1), Atom(eps2, ("a", "c"),      +1)),
        (Atom(eps3, ("a", "b", "c"), +1), Atom(eps3, ("a", "b", "d"), +1)),
    ]


def test_orbit_pair_agrees_with_brute(pairs_for_xval, ctx4):
    """orbit_pair() and orbit_brute_pair() return the same set of pairs."""
    g = ctx4.the_group
    for u, v in pairs_for_xval:
        assert set(g.orbit_pair(u, v)) == set(g.orbit_brute_pair(u, v))


def test_in_orbit_pair_agrees_with_brute(dot, eps3, ctx4):
    """in_orbit_pair() and in_orbit_brute_pair() agree on all orbit elements."""
    g = ctx4.the_group
    u = Atom(dot,  ("a", "b"),      +1)
    v = Atom(eps3, ("a", "b", "c"), +1)
    orb = g.orbit_brute_pair(u, v)
    rep = (u, v)
    for pair in orb:
        assert g.in_orbit_pair(pair, rep) == g.in_orbit_brute_pair(pair, rep)




# ---------------------------------------------------------------------------
# Single-atom orbit, stabiliser_size, orbit_size, in_orbit
# ---------------------------------------------------------------------------

def test_orbit_single_contains_representative(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    assert u in g.orbit(u)

def test_orbit_single_first_element_is_representative(dot, ctx4):
    """orbit() first element is u itself (from identity permutation)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    assert g.orbit(u)[0] == u

def test_orbit_single_no_duplicates(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    orb = g.orbit(u)
    assert len(orb) == len(set(orb))

def test_orbit_single_length_matches_orbit_size(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    assert len(g.orbit(u)) == g.orbit_size(u)

def test_orbit_size_times_stab_equals_order_single(dot, eps2, eps3, ctx4):
    """orbit_size(u) * stabiliser_size(u) == order() for single atoms."""
    g = ctx4.the_group
    atoms = [
        Atom(dot,  ("a", "b"),      sign=+1),
        Atom(dot,  ("a", "c"),      sign=+1),
        Atom(eps2, ("a", "b"),      sign=+1),
        Atom(eps2, ("a", "c"),      sign=+1),
        Atom(eps3, ("a", "b", "c"), sign=+1),
    ]
    for u in atoms:
        assert g.orbit_size(u) * g.stabiliser_size(u) == g.order(), (
            f"Orbit-stabiliser theorem violated for single atom {u!r}"
        )

def test_orbit_single_agrees_with_brute(dot, eps2, ctx4):
    """orbit() and orbit_brute() return the same set for single atoms."""
    g = ctx4.the_group
    for u in [
        Atom(dot,  ("a", "b"),      sign=+1),
        Atom(eps2, ("a", "b"),      sign=+1),
        Atom(eps2, ("b", "c"),      sign=+1),
    ]:
        assert set(g.orbit(u)) == set(g.orbit_brute(u))

def test_in_orbit_single_self(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    assert g.in_orbit(u, u)

def test_in_orbit_single_orbit_element(dot, ctx4):
    """Every element of orbit(u) satisfies in_orbit(elem, u)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    g = ctx4.the_group
    for a in g.orbit(u):
        assert g.in_orbit(a, u)

def test_in_orbit_single_agrees_with_brute(dot, eps2, ctx4):
    """in_orbit() and in_orbit_brute() agree on all single-atom orbit elements."""
    g = ctx4.the_group
    u = Atom(eps2, ("a", "b"), +1)
    for a in g.orbit_brute(u):
        assert g.in_orbit(a, u) == g.in_orbit_brute(a, u)

def test_stabiliser_size_single_symmetric(dot, ctx4):
    """SYMMETRIC atom with k labels: stab = k! * (n-k)!"""
    g = ctx4.the_group  # n=4
    u = Atom(dot, ("a", "b"), sign=+1)  # k=2
    assert g.stabiliser_size(u) == math.factorial(2) * math.factorial(2)

def test_stabiliser_size_single_antisymmetric(eps2, ctx4):
    """ANTISYMMETRIC atom with k labels: stab = E(k) * (n-k)!"""
    g = ctx4.the_group  # n=4
    u = Atom(eps2, ("a", "b"), sign=+1)  # k=2
    assert g.stabiliser_size(u) == 1 * math.factorial(2)  # E(2)=1, (4-2)!=2

def test_orbit_size_single_two_groups(dot, ctx_two_groups):
    """orbit_size(u) * stabiliser_size(u) == order() in a two-group context."""
    g = ctx_two_groups.the_group
    u = Atom(dot, ("a", "p"), sign=+1)
    assert g.orbit_size(u) * g.stabiliser_size(u) == g.order()
