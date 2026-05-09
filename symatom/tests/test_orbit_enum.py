"""
Tests for OrbitEnumerator implementations.

BruteForceOrbitEnumerator is tested fully.  DirectOrbitEnumerator tests are
parametrised alongside it but skipped until the implementation is complete.
When DirectOrbitEnumerator is implemented, remove the skip markers and the
parametrised tests will serve as the correctness cross-check against brute
force.
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
    pytest.param(DirectOrbitEnumerator(),     id="direct",
                 marks=pytest.mark.skip(reason="DirectOrbitEnumerator not yet implemented")),
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
    """orbit_elements returns exactly orbit_size pairs for mixed dot+eps3 PairFlavours."""
    fo_list = repL(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
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
    """Works with a two-group context."""
    electrons = VectorGroup("electrons", ("a", "b", "c"))
    muons     = VectorGroup("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    fo_list   = repL(ctx, [dot, eps3])
    group_sizes = tuple(g.size for g in ctx.groups)
    for pf in canonical_pair_flavours(fo_list, ctx):
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

def test_brute_force_matches_pf_orbit_elements(dot, ctx):
    """BruteForceOrbitEnumerator must agree with pf.orbit_elements() exactly."""
    fo_list = repL(ctx, [dot])
    bf = BruteForceOrbitEnumerator()
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert bf.orbit_elements(pf, ctx) == pf.orbit_elements(ctx)


# ---------------------------------------------------------------------------
# DirectOrbitEnumerator — raises until implemented
# ---------------------------------------------------------------------------

def test_direct_raises_not_implemented(dot, ctx):
    fl2 = Flavour((2,))
    pf  = PairFlavour(op_u=dot, flavour_u=fl2, op_v=dot, flavour_v=fl2, overlap=(0,))
    with pytest.raises(NotImplementedError):
        DirectOrbitEnumerator().orbit_elements(pf, ctx)
