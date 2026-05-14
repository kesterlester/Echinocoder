"""
Tests verifying that the symatom stack handles VectorTypes with zero labels.

An empty group represents a particle species with zero members.  It is
mathematically harmless (S_0 is the trivial group, 0!=1, C(0,0)=1) and
carries meaningful semantic identity: a context with (electrons=3, muons=0)
is distinct from one with (electrons=0, muons=3) even if the non-zero groups
are numerically identical.

These tests exercise the full machinery — Context, repS,
canonical_pair_flavours, OrbitEnumerators, orbit() — when one or more groups
are empty, ensuring nothing silently breaks.
"""
import pytest
from symatom import (
    ArgumentSymmetry, Operation, VectorType, Atom,
    Context,
    BruteForceOrbitEnumerator, DirectOrbitEnumerator,
    repS, canonical_pair_flavours,
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
def electrons():
    return VectorType("electrons", ("a", "b", "c"))

@pytest.fixture
def empty_muons():
    return VectorType("muons", ())

@pytest.fixture
def ctx_one_empty(electrons, empty_muons):
    return Context((electrons, empty_muons))

@pytest.fixture
def ctx_all_empty():
    return Context((VectorType("electrons", ()), VectorType("muons", ())))


# ---------------------------------------------------------------------------
# Context construction
# ---------------------------------------------------------------------------

def test_context_with_one_empty_group(ctx_one_empty):
    assert len(ctx_one_empty.types) == 2
    assert ctx_one_empty.types[1].size == 0

def test_context_repr_includes_empty_group(ctx_one_empty):
    r = repr(ctx_one_empty)
    assert "muons[0]" in r

def test_context_all_labels_omits_empty_group(ctx_one_empty):
    assert set(ctx_one_empty.all_labels) == {"a", "b", "c"}

def test_context_all_empty_groups():
    ctx = Context((VectorType("e", ()), VectorType("mu", ())))
    assert ctx.all_labels == ()


# ---------------------------------------------------------------------------
# repS — FlavouredOperator generation
# ---------------------------------------------------------------------------

def test_repl_with_empty_group_produces_fos(dot, eps3, ctx_one_empty):
    """repS should still produce FOs; empty group contributes 0-count entries."""
    fo_list = repS(ctx_one_empty, [dot, eps3])
    assert len(fo_list) > 0

def test_repl_with_empty_group_fo_flavour_counts(dot, ctx_one_empty):
    """Every FO for a rank-2 op in a (3, 0) context has second count = 0."""
    fo_list = repS(ctx_one_empty, [dot])
    for fo in fo_list:
        assert fo.flavour.counts[1] == 0   # muon count always 0

def test_repl_all_empty_gives_no_fos(dot, ctx_all_empty):
    """No atoms can be formed when all groups are empty."""
    fo_list = repS(ctx_all_empty, [dot])
    assert fo_list == []


# ---------------------------------------------------------------------------
# canonical_pair_flavours
# ---------------------------------------------------------------------------

def test_cpf_with_empty_group(dot, eps3, ctx_one_empty):
    fo_list = repS(ctx_one_empty, [dot, eps3])
    pfs = canonical_pair_flavours(fo_list, ctx_one_empty)
    assert len(pfs) > 0

def test_cpf_empty_group_overlap_is_zero(dot, ctx_one_empty):
    """The overlap for an empty group is always 0."""
    fo_list = repS(ctx_one_empty, [dot])
    for pf in canonical_pair_flavours(fo_list, ctx_one_empty):
        assert pf.overlap[1] == 0   # muon group is always 0

def test_cpf_all_empty_gives_no_pfs(dot, ctx_all_empty):
    fo_list = repS(ctx_all_empty, [dot])
    pfs = canonical_pair_flavours(fo_list, ctx_all_empty)
    assert pfs == []


# ---------------------------------------------------------------------------
# OrbitEnumerators with an empty group
# ---------------------------------------------------------------------------

def test_orbit_enumerator_size_matches_with_empty_group_brute_force(dot, eps3, ctx_one_empty):
    """BruteForce returns exactly orbit_size pairs for every PairFlavour when a group is empty."""
    fo_list = repS(ctx_one_empty, [dot, eps3])
    type_sizes = tuple(g.size for g in ctx_one_empty.types)
    for pf in canonical_pair_flavours(fo_list, ctx_one_empty):
        assert len(BruteForceOrbitEnumerator().orbit_elements(pf, ctx_one_empty)) == pf.orbit_size(type_sizes)

def test_orbit_enumerator_size_matches_with_empty_group_direct(dot, eps3, ctx_one_empty):
    """Direct returns exactly orbit_size pairs for every PairFlavour when a group is empty."""
    fo_list = repS(ctx_one_empty, [dot, eps3])
    type_sizes = tuple(g.size for g in ctx_one_empty.types)
    for pf in canonical_pair_flavours(fo_list, ctx_one_empty):
        assert len(DirectOrbitEnumerator().orbit_elements(pf, ctx_one_empty)) == pf.orbit_size(type_sizes)

def test_orbit_enumerators_agree_with_empty_group_ss(dot, ctx_one_empty):
    """BruteForce and Direct agree for SS PairFlavours when a group is empty."""
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    fo_list = repS(ctx_one_empty, [dot])
    for pf in canonical_pair_flavours(fo_list, ctx_one_empty):
        assert set(bf.orbit_elements(pf, ctx_one_empty)) == \
               set(direct.orbit_elements(pf, ctx_one_empty))

def test_orbit_enumerators_agree_with_empty_group(dot, eps3, ctx_one_empty):
    """BruteForce and Direct produce identical element sets when a group is empty."""
    bf     = BruteForceOrbitEnumerator()
    direct = DirectOrbitEnumerator()
    fo_list = repS(ctx_one_empty, [dot, eps3])
    for pf in canonical_pair_flavours(fo_list, ctx_one_empty):
        assert set(bf.orbit_elements(pf, ctx_one_empty)) == \
               set(direct.orbit_elements(pf, ctx_one_empty))

