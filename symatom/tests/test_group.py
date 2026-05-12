"""
Tests for symatom.group: TheGroup and SignCorrelationType.

TheGroup wraps a Context and provides orbit enumeration, stabiliser size,
orbit-membership testing, sign-correlation typing, and companion-orbit
enumeration for atom-pairs.

All methods are currently brute-force O(∏n_g!).  These tests verify
correctness (not performance) and serve as the permanent reference against
which faster future implementations will be validated.
"""
import math
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Atom, Context,
    TheGroup, SignCorrelationType,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def dot():
    return Operation("dot", rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)

@pytest.fixture
def eps3():
    return Operation("eps3", rank=3, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)

@pytest.fixture
def eps2():
    """Rank-2 antisymmetric operation — useful for small-n overlap tests."""
    return Operation("eps2", rank=2, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)

@pytest.fixture
def electrons4():
    return VectorGroup("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def electrons3():
    return VectorGroup("electrons", ("a", "b", "c"))

@pytest.fixture
def ctx4(electrons4):
    return Context((electrons4,))

@pytest.fixture
def ctx3(electrons3):
    return Context((electrons3,))

@pytest.fixture
def ctx_two_groups():
    e = VectorGroup("electrons", ("a", "b", "c"))
    mu = VectorGroup("muons", ("p", "q"))
    return Context((e, mu))


# ---------------------------------------------------------------------------
# order()
# ---------------------------------------------------------------------------

def test_order_single_group_n4(ctx4):
    g = TheGroup(ctx4)
    assert g.order() == math.factorial(4)  # 24

def test_order_single_group_n3(ctx3):
    g = TheGroup(ctx3)
    assert g.order() == math.factorial(3)  # 6

def test_order_two_groups(ctx_two_groups):
    g = TheGroup(ctx_two_groups)
    # 3! * 2! = 12
    assert g.order() == math.factorial(3) * math.factorial(2)

def test_order_empty_context():
    ctx = Context((VectorGroup("e", ()),))
    g = TheGroup(ctx)
    assert g.order() == 1  # 0! = 1


# ---------------------------------------------------------------------------
# orbit() — basic properties
# ---------------------------------------------------------------------------

def test_orbit_contains_representative(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    orb = g.orbit(u, v)
    assert (u, v) in orb

def test_orbit_first_element_is_representative(dot, ctx4):
    """orbit() first element is the input pair (from identity permutation)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    g = TheGroup(ctx4)
    orb = g.orbit(u, v)
    assert orb[0] == (u, v)

def test_orbit_no_duplicates(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    orb = g.orbit(u, v)
    assert len(orb) == len(set(orb))

def test_orbit_length_matches_orbit_size(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    assert len(g.orbit(u, v)) == g.orbit_size(u, v)

def test_orbit_all_pairs_are_atom_tuples(dot, eps3, ctx4):
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(dot, ("a", "b"), sign=+1)
    g = TheGroup(ctx4)
    for pair in g.orbit(u, v):
        assert isinstance(pair, tuple) and len(pair) == 2
        assert isinstance(pair[0], Atom) and isinstance(pair[1], Atom)


# ---------------------------------------------------------------------------
# stabiliser_size() and orbit_size()
# ---------------------------------------------------------------------------

def test_orbit_size_times_stab_equals_order(dot, eps3, ctx4):
    """orbit_size * stabiliser_size == order (orbit-stabiliser theorem)."""
    g = TheGroup(ctx4)
    pairs = [
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot,  ("c", "d"), sign=+1)),
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot,  ("a", "c"), sign=+1)),
        (Atom(dot,  ("a", "b"), sign=+1), Atom(dot,  ("a", "b"), sign=+1)),
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(eps3, ("a", "b", "c"), sign=+1)),
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(dot, ("a", "b"), sign=+1)),
    ]
    for u, v in pairs:
        assert g.orbit_size(u, v) * g.stabiliser_size(u, v) == g.order(), (
            f"Orbit-stabiliser theorem violated for {u!r}, {v!r}"
        )

def test_orbit_size_ss_zero_overlap(dot, ctx4):
    """SS no-overlap: C(4,2)*C(2,2) = 6."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    assert g.orbit_size(u, v) == 6

def test_orbit_size_ss_partial_overlap(dot, ctx4):
    """SS overlap=1: C(4,1)*C(3,1)*C(2,1) = 24? No: count=C(4,1)*C(3,1)*C(2,1)=...

    For SS s=1, a=1, b=1, r=1, n=4:
      stab = 1!*1!*1!*1! = 1
      orbit_size = 4!/1 = 24
    But count = C(4,1)*C(3,1)*C(2,1) = 4*3*2 = 24 — SS count == orbit_size. Correct.
    """
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    g = TheGroup(ctx4)
    assert g.orbit_size(u, v) == 24

def test_orbit_size_ss_full_overlap(dot, ctx4):
    """SS full overlap (u==v): stab = 2!*0!*0!*2! = 4, orbit_size = 24/4 = 6."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "b"), sign=+1)
    g = TheGroup(ctx4)
    assert g.orbit_size(u, v) == 6

def test_orbit_size_aa_full_overlap(eps3, ctx3):
    """AA full overlap: stab = E(3) = 3, orbit_size = 6/3 = 2."""
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx3)
    assert g.orbit_size(u, v) == 2

def test_orbit_size_aa_zero_overlap(eps2, ctx4):
    """AA zero overlap, rank=2, n=4: stab=1, orbit_size=24."""
    u = Atom(eps2, ("a", "b"), sign=+1)
    v = Atom(eps2, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    assert g.orbit_size(u, v) == 24

def test_orbit_size_matches_len_orbit(eps3, ctx3):
    """orbit_size computed algebraically must equal len(orbit()) for all test cases."""
    g = TheGroup(ctx3)
    test_pairs = [
        (Atom(eps3, ("a", "b", "c"), sign=+1), Atom(eps3, ("a", "b", "c"), sign=+1)),
    ]
    for u, v in test_pairs:
        assert g.orbit_size(u, v) == len(g.orbit(u, v))


# ---------------------------------------------------------------------------
# orbit() agrees with pf.orbit_elements() for canonical pairs
# ---------------------------------------------------------------------------

def test_orbit_agrees_with_pf_orbit_elements(dot, eps3):
    """
    For the canonical representative of every PairFlavour, TheGroup.orbit and
    pf.orbit_elements must return the same set of pairs.
    """
    from symatom import repL, canonical_pair_flavours
    electrons = VectorGroup("electrons", ("a", "b", "c", "d"))
    ctx = Context((electrons,))
    g = TheGroup(ctx)
    for pf in canonical_pair_flavours(repL(ctx, [dot, eps3]), ctx):
        # canonical representative: first labels of each group up to flavour count
        ku = pf.flavour_u.counts[0]
        kv = pf.flavour_v.counts[0]
        s  = pf.overlap[0]
        u_labels = electrons.labels[:ku]
        v_labels = electrons.labels[:s] + electrons.labels[ku:ku + kv - s]
        u = Atom(pf.op_u, u_labels, sign=+1)
        v = Atom(pf.op_v, v_labels, sign=+1)
        assert set(g.orbit(u, v)) == set(pf.orbit_elements(ctx)), (
            f"orbit disagrees with pf.orbit_elements for {pf!r}"
        )


# ---------------------------------------------------------------------------
# in_orbit()
# ---------------------------------------------------------------------------

def test_in_orbit_self(dot, ctx4):
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    assert g.in_orbit((u, v), (u, v))

def test_in_orbit_orbit_element(dot, ctx4):
    """Every element of orbit(u,v) is in_orbit((u,v), (u,v))."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    g = TheGroup(ctx4)
    for pair in g.orbit(u, v):
        assert g.in_orbit(pair, (u, v))

def test_not_in_orbit_different_flavour(dot, ctx4):
    """A pair with different label counts is not in the orbit."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    # Different PairFlavour: overlap=1 instead of 0
    u2 = Atom(dot, ("a", "b"), sign=+1)
    v2 = Atom(dot, ("a", "c"), sign=+1)
    g = TheGroup(ctx4)
    assert not g.in_orbit((u2, v2), (u, v))


# ---------------------------------------------------------------------------
# sign_correlation_type()
# ---------------------------------------------------------------------------

def test_sign_correlation_type_ss(dot, ctx4):
    """SS pair: TYPE_11 (signs always +1, nothing flips)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    assert g.sign_correlation_type(u, v) == SignCorrelationType.TYPE_11

def test_sign_correlation_type_ss_overlap(dot, ctx4):
    """SS pair with overlap: still TYPE_11."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("a", "c"), sign=+1)
    g = TheGroup(ctx4)
    assert g.sign_correlation_type(u, v) == SignCorrelationType.TYPE_11

def test_sign_correlation_type_aa_full_overlap(eps3, ctx3):
    """AA full overlap: signs couple (only (++) and (--) in orbit) → TYPE_11."""
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx3)
    assert g.sign_correlation_type(u, v) == SignCorrelationType.TYPE_11

def test_sign_correlation_type_aa_zero_overlap(eps2, ctx4):
    """AA zero overlap: signs are independent → TYPE_22."""
    u = Atom(eps2, ("a", "b"), sign=+1)
    v = Atom(eps2, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    assert g.sign_correlation_type(u, v) == SignCorrelationType.TYPE_22

def test_sign_correlation_type_as(eps3, dot, ctx4):
    """AS pair (u antisymmetric, v symmetric): u flips freely → TYPE_21."""
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(dot,  ("a", "b"),      sign=+1)
    g = TheGroup(ctx4)
    assert g.sign_correlation_type(u, v) == SignCorrelationType.TYPE_21

def test_sign_correlation_type_sa(dot, eps3, ctx4):
    """SA pair (u symmetric, v antisymmetric): v flips freely → TYPE_12."""
    u = Atom(dot,  ("a", "b"),      sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx4)
    assert g.sign_correlation_type(u, v) == SignCorrelationType.TYPE_12


# ---------------------------------------------------------------------------
# companion_orbits()
# ---------------------------------------------------------------------------

def test_companion_orbits_ss_one_orbit(dot, ctx4):
    """SS pair has only one companion orbit (signs are always +1)."""
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    companions = g.companion_orbits(u, v)
    assert len(companions) == 1

def test_companion_orbits_primary_is_first(dot, eps3, ctx4):
    """Primary orbit is always the first element."""
    u = Atom(dot,  ("a", "b"),      sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx4)
    companions = g.companion_orbits(u, v)
    primary = g.orbit(u, v)
    assert set(companions[0]) == set(primary)

def test_companion_orbits_aa_full_overlap_two_orbits(eps3, ctx3):
    """AA full overlap: only (++) orbit and (+−)/(−+) orbit are distinct → 2 companions."""
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx3)
    companions = g.companion_orbits(u, v)
    assert len(companions) == 2

def test_companion_orbits_aa_full_overlap_disjoint(eps3, ctx3):
    """The two companion orbits for AA full overlap must be disjoint."""
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx3)
    c0, c1 = g.companion_orbits(u, v)
    assert set(c0).isdisjoint(set(c1))

def test_companion_orbits_aa_full_overlap_union_size(eps3, ctx3):
    """The union of the two AA full-overlap companion orbits covers all sign variants."""
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx3)
    companions = g.companion_orbits(u, v)
    union = set()
    for orb in companions:
        union.update(orb)
    # 2 companions × 2 pairs each = 4 total sign-variant orbit elements
    assert len(union) == 4

def test_companion_orbits_aa_zero_overlap_one_orbit(eps2, ctx4):
    """AA zero overlap: all 4 sign variants are in the same orbit → 1 companion."""
    u = Atom(eps2, ("a", "b"), sign=+1)
    v = Atom(eps2, ("c", "d"), sign=+1)
    g = TheGroup(ctx4)
    companions = g.companion_orbits(u, v)
    assert len(companions) == 1

def test_companion_orbits_no_duplicates(dot, eps3, ctx4):
    """No orbit appears twice in companion_orbits()."""
    u = Atom(dot,  ("a", "b"),      sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx4)
    companions = g.companion_orbits(u, v)
    orbit_keys = [frozenset(orb) for orb in companions]
    assert len(orbit_keys) == len(set(orbit_keys))

def test_companion_orbits_sa_one_orbit_shared_labels(dot, eps3, ctx4):
    """
    SA pair where u and v share labels (overlap>0): still 1 companion orbit.

    u = dot(a,b), v = eps3(a,b,c) with n=4.  The permutation (a↔b) fixes
    dot(a,b) (symmetric) and flips eps3(a,b,c) → -eps3(a,b,c), so
    (dot(a,b), -eps3(a,b,c)) lies in the primary orbit.  Therefore
    orbit(u, -v) == orbit(u, v) → only 1 companion orbit.
    """
    u = Atom(dot,  ("a", "b"),      sign=+1)
    v = Atom(eps3, ("a", "b", "c"), sign=+1)
    g = TheGroup(ctx4)
    companions = g.companion_orbits(u, v)
    assert len(companions) == 1

def test_companion_orbits_as_one_orbit_shared_labels(eps3, dot, ctx4):
    """
    AS pair where u and v share labels: still 1 companion orbit.

    u = eps3(a,b,c), v = dot(a,b) with n=4.  The permutation (a↔b)
    flips eps3(a,b,c) → -eps3(a,b,c) while fixing dot(a,b), so
    (-eps3(a,b,c), dot(a,b)) is in the primary orbit → 1 companion.
    """
    u = Atom(eps3, ("a", "b", "c"), sign=+1)
    v = Atom(dot,  ("a", "b"),      sign=+1)
    g = TheGroup(ctx4)
    companions = g.companion_orbits(u, v)
    assert len(companions) == 1

def test_companion_orbits_sa_two_orbits_rank1(eps2):
    """
    SA pair where no permutation can flip v while keeping the *same* u atom:
    → 2 distinct companion orbits.

    u = rank-1 symmetric op on label 'a' (rank 1, so u only ever uses a single
    label); v = rank-2 antisymmetric op on ('a','b').  With n=2:
      orbit(u, v)  = {(+dot1(a), +eps2(a,b)), (+dot1(b), -eps2(a,b))}
      orbit(u, -v) = {(+dot1(a), -eps2(a,b)), (+dot1(b), +eps2(a,b))}
    The two orbits are disjoint, giving 2 companion orbits.
    """
    dot1 = Operation("dot1", rank=1, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)
    ctx2 = Context((VectorGroup("e", ("a", "b")),))
    g    = TheGroup(ctx2)
    u = Atom(dot1, ("a",),      sign=+1)
    v = Atom(eps2, ("a", "b"), sign=+1)
    companions = g.companion_orbits(u, v)
    assert len(companions) == 2


# ---------------------------------------------------------------------------
# Two-group context
# ---------------------------------------------------------------------------

def test_orbit_size_two_groups(dot, ctx_two_groups):
    """orbit_size * stabiliser_size == order in a two-group context."""
    g = TheGroup(ctx_two_groups)
    u = Atom(dot, ("a", "p"), sign=+1)
    v = Atom(dot, ("b", "q"), sign=+1)
    assert g.orbit_size(u, v) * g.stabiliser_size(u, v) == g.order()

def test_orbit_length_two_groups(dot, ctx_two_groups):
    """len(orbit) == orbit_size in a two-group context."""
    g = TheGroup(ctx_two_groups)
    u = Atom(dot, ("a", "b"), sign=+1)
    v = Atom(dot, ("p", "q"), sign=+1)
    assert len(g.orbit(u, v)) == g.orbit_size(u, v)


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


def test_sign_correlation_type_agrees_with_brute(pairs_for_xval, ctx4):
    """Algebraic sign_correlation_type() == sign_correlation_type_brute() for all pairs."""
    g = TheGroup(ctx4)
    for u, v in pairs_for_xval:
        algebraic = g.sign_correlation_type(u, v)
        brute     = g.sign_correlation_type_brute(u, v)
        assert algebraic == brute, (
            f"sign_correlation_type mismatch for ({u!r}, {v!r}): "
            f"algebraic={algebraic}, brute={brute}"
        )


def test_orbit_agrees_with_brute(pairs_for_xval, ctx4):
    """orbit() and orbit_brute() return the same set of pairs."""
    g = TheGroup(ctx4)
    for u, v in pairs_for_xval:
        assert set(g.orbit(u, v)) == set(g.orbit_brute(u, v))


def test_in_orbit_agrees_with_brute(dot, eps3, ctx4):
    """in_orbit() and in_orbit_brute() agree on all orbit elements."""
    g = TheGroup(ctx4)
    u = Atom(dot,  ("a", "b"),      +1)
    v = Atom(eps3, ("a", "b", "c"), +1)
    orb = g.orbit_brute(u, v)
    rep = (u, v)
    for pair in orb:
        assert g.in_orbit(pair, rep) == g.in_orbit_brute(pair, rep)


def test_companion_orbits_agrees_with_brute(pairs_for_xval, ctx4):
    """companion_orbits() and companion_orbits_brute() return the same orbit sets."""
    g = TheGroup(ctx4)
    for u, v in pairs_for_xval:
        primary     = [frozenset(orb) for orb in g.companion_orbits(u, v)]
        brute       = [frozenset(orb) for orb in g.companion_orbits_brute(u, v)]
        assert primary == brute, (
            f"companion_orbits mismatch for ({u!r}, {v!r})"
        )


def test_sign_correlation_type_brute_aa_partial_overlap(eps2, ctx4):
    """sign_correlation_type_brute gives TYPE_22 for AA partial overlap."""
    g = TheGroup(ctx4)
    u = Atom(eps2, ("a", "b"), +1)
    v = Atom(eps2, ("a", "c"), +1)
    assert g.sign_correlation_type_brute(u, v) == SignCorrelationType.TYPE_22


def test_sign_correlation_type_algebraic_vs_brute_sa_rank1(eps2):
    """Algebraic and brute sign_correlation_type agree for rank-1 SA pairs."""
    dot1 = Operation("dot1", rank=1, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)
    ctx2 = Context((VectorGroup("e", ("a", "b")),))
    g = TheGroup(ctx2)
    u = Atom(dot1, ("a",),      +1)
    v = Atom(eps2, ("a", "b"), +1)
    assert g.sign_correlation_type(u, v) == g.sign_correlation_type_brute(u, v)


def test_sign_correlation_type_sweep_two_groups(dot, eps3, eps2):
    """Algebraic and brute agree across all PairFlavours in a two-group context."""
    from symatom import repL, canonical_pair_flavours
    from symatom.rep import PairFlavour
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    g         = TheGroup(ctx)
    fo_list   = repL(ctx, [dot, eps3, eps2])
    for pf in canonical_pair_flavours(fo_list, ctx):
        # Build the canonical representative pair for this PairFlavour
        u_labels, v_labels = [], []
        for grp, ku, kv, s in zip(ctx.groups, pf.flavour_u.counts,
                                   pf.flavour_v.counts, pf.overlap):
            u_labels.extend(grp.labels[:ku])
            v_labels.extend(grp.labels[:s])
            v_labels.extend(grp.labels[ku:ku + kv - s])
        u = Atom(pf.op_u, tuple(u_labels), +1)
        v = Atom(pf.op_v, tuple(v_labels), +1)
        assert g.sign_correlation_type(u, v) == g.sign_correlation_type_brute(u, v), (
            f"Mismatch for PairFlavour {pf!r}"
        )
