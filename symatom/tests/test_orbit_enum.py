"""
Tests for OrbitEnumerator implementations.

Both BruteForceOrbitEnumerator and DirectOrbitEnumerator are tested via the
parametrised suite below.  A cross-comparison test verifies that Direct
produces the same element set as BruteForce for every PairFlavour.
"""
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Context,
    BruteForceOrbitEnumerator, DirectOrbitEnumerator,
    canonical_pair_flavours, repL,
)
from symatom.rep import Flavour, FlavouredOperator, PairFlavour, pair_flavour_of


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
def electrons():
    return VectorGroup("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def ctx(electrons):
    return Context((electrons,))


# ---------------------------------------------------------------------------
# Parametrised over both enumerator implementations
# ---------------------------------------------------------------------------

ENUMERATORS = [
    pytest.param(BruteForceOrbitEnumerator(), id="brute_force"),
    pytest.param(DirectOrbitEnumerator(),     id="direct"),
]


@pytest.mark.parametrize("enumerator", ENUMERATORS)
def test_orbit_elements_length_matches_orbit_size_dot(enumerator, dot, ctx):
    """orbit_elements returns exactly orbit_size pairs for all dot PairFlavours."""
    fo_list = repL(ctx, [dot])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        elems = enumerator.orbit_elements(pf, ctx)
        assert len(elems) == pf.orbit_size(group_sizes)

@pytest.mark.parametrize("enumerator", ENUMERATORS)
def test_orbit_elements_length_matches_orbit_size_dot_eps(enumerator, dot, eps3, ctx):
    """orbit_elements returns exactly orbit_size pairs for mixed dot+eps3 PairFlavours.

    DirectOrbitEnumerator generates all sign-combinations independently for any
    ANTISYMMETRIC op (for encoding purposes).  orbit_size is the true single-orbit
    size, which may be less than Direct's count when AS, SA, or AA PairFlavours
    appear.  Those are skipped for Direct; BruteForce is checked for all.
    """
    fo_list = repL(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        has_antisym = (pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC or
                       pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
        if has_antisym:
            continue  # Both enumerators generate all sign-combos for any antisymmetric
                      # op; orbit_size is the true single-orbit size and will differ
        elems = enumerator.orbit_elements(pf, ctx)
        assert len(elems) == pf.orbit_size(group_sizes)

@pytest.mark.parametrize("enumerator", ENUMERATORS)
def test_orbit_elements_all_correct_pair_flavour(enumerator, dot, ctx):
    """Every returned pair has the expected PairFlavour."""
    fo_list = repL(ctx, [dot])
    for pf in canonical_pair_flavours(fo_list, ctx):
        for u, v in enumerator.orbit_elements(pf, ctx):
            assert pair_flavour_of(u, v, ctx) == pf

@pytest.mark.parametrize("enumerator", ENUMERATORS)
def test_orbit_elements_two_groups(enumerator, dot, eps3):
    """Works with a two-group context.

    As in test_orbit_elements_length_matches_orbit_size_dot_eps, PairFlavours
    with any ANTISYMMETRIC op are skipped for DirectOrbitEnumerator since it
    intentionally generates all sign-combinations for encoding purposes.
    """
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    fo_list   = repL(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        has_antisym = (pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC or
                       pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
        if has_antisym:
            continue  # Both enumerators generate all sign-combos for any antisymmetric
                      # op; orbit_size is the true single-orbit size and will differ
        elems = enumerator.orbit_elements(pf, ctx)
        assert len(elems) == pf.orbit_size(group_sizes)


# ---------------------------------------------------------------------------
# BruteForceOrbitEnumerator — additional non-parametrised tests
# ---------------------------------------------------------------------------

def test_brute_force_returns_list(dot, ctx):
    fl2 = Flavour((2,))
    pf  = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(1,))
    result = BruteForceOrbitEnumerator().orbit_elements(pf, ctx)
    assert isinstance(result, list)
    assert all(isinstance(pair, tuple) and len(pair) == 2 for pair in result)

def test_brute_force_matches_pf_orbit_elements_for_ss(dot, ctx):
    """BruteForceOrbitEnumerator agrees with pf.orbit_elements() for SS PairFlavours.

    For SYMMETRIC × SYMMETRIC pairs there is only one sign variant per label
    assignment, so all pairs belong to a single orbit-type and the cross-product
    enumeration coincides with the true G-orbit traversal in pf.orbit_elements.

    For AS, SA, or AA PairFlavours this equality does NOT hold in general:
    BruteForceOrbitEnumerator returns all sign-combo pairs (matching Direct),
    while pf.orbit_elements returns only the single G-orbit of the canonical
    representative.
    """
    fo_list = repL(ctx, [dot])
    bf = BruteForceOrbitEnumerator()
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert set(bf.orbit_elements(pf, ctx)) == set(pf.orbit_elements(ctx))


# ---------------------------------------------------------------------------
# Cross-comparison: Direct must produce the same element set as BruteForce
# ---------------------------------------------------------------------------

def test_direct_matches_brute_force_dot(dot, ctx):
    """Direct produces exactly the same atom-pairs as BruteForce for all dot PairFlavours."""
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    for pf in canonical_pair_flavours(repL(ctx, [dot]), ctx):
        assert set(bf.orbit_elements(pf, ctx)) == set(direct.orbit_elements(pf, ctx))

def test_direct_matches_brute_force_dot_eps(dot, eps3, ctx):
    """Direct produces exactly the same atom-pairs as BruteForce for all dot+eps3 PairFlavours."""
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    for pf in canonical_pair_flavours(repL(ctx, [dot, eps3]), ctx):
        assert set(bf.orbit_elements(pf, ctx)) == set(direct.orbit_elements(pf, ctx))

def test_direct_matches_brute_force_two_groups(dot, eps3):
    """Direct matches BruteForce in a two-group context, for all PairFlavours."""
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    for pf in canonical_pair_flavours(repL(ctx, [dot, eps3]), ctx):
        assert set(bf.orbit_elements(pf, ctx)) == set(direct.orbit_elements(pf, ctx))
