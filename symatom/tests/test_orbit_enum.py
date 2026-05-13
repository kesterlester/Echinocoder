"""
Tests for OrbitEnumerator implementations.

BruteForceOrbitEnumerator delegates to pf.orbit_elements() and is the
permanent reference implementation — correct by construction for all
PairFlavours.

DirectOrbitEnumerator (Step 2) uses TheGroup(context).orbit() to return the
true G-orbit.  Both enumerators must return the same set for every PairFlavour;
the parametrised cross-comparison tests at the bottom verify this.
"""
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorGroup, Context,
    BruteForceOrbitEnumerator, DirectOrbitEnumerator,
    canonical_pair_flavours, repS,
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
# BruteForceOrbitEnumerator — correct for all PairFlavours, no caveats
# ---------------------------------------------------------------------------

def test_brute_force_length_matches_orbit_size_dot(dot, ctx):
    """BruteForce returns exactly orbit_size pairs for all dot PairFlavours."""
    fo_list = repS(ctx, [dot])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert len(BruteForceOrbitEnumerator().orbit_elements(pf, ctx)) == pf.orbit_size(group_sizes)

def test_brute_force_length_matches_orbit_size_dot_eps(dot, eps3, ctx):
    """BruteForce returns exactly orbit_size pairs for all dot+eps3 PairFlavours."""
    fo_list = repS(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert len(BruteForceOrbitEnumerator().orbit_elements(pf, ctx)) == pf.orbit_size(group_sizes)

def test_brute_force_length_matches_orbit_size_two_groups(dot, eps3):
    """BruteForce returns exactly orbit_size pairs in a two-group context."""
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    fo_list   = repS(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert len(BruteForceOrbitEnumerator().orbit_elements(pf, ctx)) == pf.orbit_size(group_sizes)

def test_brute_force_all_correct_pair_flavour(dot, ctx):
    """Every pair returned by BruteForce has the expected PairFlavour."""
    fo_list = repS(ctx, [dot])
    for pf in canonical_pair_flavours(fo_list, ctx):
        for u, v in BruteForceOrbitEnumerator().orbit_elements(pf, ctx):
            assert pair_flavour_of(u, v, ctx) == pf

def test_brute_force_equals_pf_orbit_elements(dot, eps3, ctx):
    """BruteForce is exactly pf.orbit_elements() for all PairFlavours.

    BruteForceOrbitEnumerator delegates to pf.orbit_elements, so this is true
    by construction.  The test is kept as a contract: if the delegation ever
    changes, this must still hold.
    """
    bf = BruteForceOrbitEnumerator()
    for pf in canonical_pair_flavours(repS(ctx, [dot, eps3]), ctx):
        assert bf.orbit_elements(pf, ctx) == pf.orbit_elements(ctx)

def test_brute_force_returns_list(dot, ctx):
    fl2 = Flavour((2,))
    pf  = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(1,))
    result = BruteForceOrbitEnumerator().orbit_elements(pf, ctx)
    assert isinstance(result, list)
    assert all(isinstance(pair, tuple) and len(pair) == 2 for pair in result)


# ---------------------------------------------------------------------------
# DirectOrbitEnumerator — correct for SS; xfail for AS/SA/AA (OrbitUnion bug)
# ---------------------------------------------------------------------------

def test_direct_length_matches_orbit_size_dot(dot, ctx):
    """Direct returns exactly orbit_size pairs for SS (dot-only) PairFlavours."""
    fo_list = repS(ctx, [dot])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert len(DirectOrbitEnumerator().orbit_elements(pf, ctx)) == pf.orbit_size(group_sizes)

def test_direct_length_matches_orbit_size_dot_eps(dot, eps3, ctx):
    """Direct returns exactly orbit_size pairs for all dot+eps3 PairFlavours."""
    fo_list = repS(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert len(DirectOrbitEnumerator().orbit_elements(pf, ctx)) == pf.orbit_size(group_sizes)

def test_direct_length_matches_orbit_size_two_groups(dot, eps3):
    """Direct returns exactly orbit_size pairs in a two-group context."""
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    fo_list   = repS(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert len(DirectOrbitEnumerator().orbit_elements(pf, ctx)) == pf.orbit_size(group_sizes)


# ---------------------------------------------------------------------------
# Cross-comparison: Direct must agree with BruteForce for every PairFlavour
# ---------------------------------------------------------------------------

def test_direct_matches_brute_force_dot(dot, ctx):
    """Direct agrees with BruteForce for all SS (dot-only) PairFlavours."""
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    for pf in canonical_pair_flavours(repS(ctx, [dot]), ctx):
        assert set(bf.orbit_elements(pf, ctx)) == set(direct.orbit_elements(pf, ctx))

def test_direct_matches_brute_force_dot_eps(dot, eps3, ctx):
    """Direct agrees with BruteForce for all dot+eps3 PairFlavours."""
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    for pf in canonical_pair_flavours(repS(ctx, [dot, eps3]), ctx):
        assert set(bf.orbit_elements(pf, ctx)) == set(direct.orbit_elements(pf, ctx))

def test_direct_matches_brute_force_two_groups(dot, eps3):
    """Direct agrees with BruteForce in a two-group context."""
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    for pf in canonical_pair_flavours(repS(ctx, [dot, eps3]), ctx):
        assert set(bf.orbit_elements(pf, ctx)) == set(direct.orbit_elements(pf, ctx))
