"""Tests for encode.py: the top-level encode() function."""
import pytest
import numpy as np
from symatom import (
    ArgumentSymmetry, VectorType, Context, Plan,
    repS, canonical_pair_flavours,
)
from symcoder import EvaluableOperation, encode, encode_brute, describe_encoding
from symcoder.pairs import eval_pair_orbit, eval_pair_orbit_positive
from symcoder.encode import _embed_compressed


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def dot():
    return EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )

@pytest.fixture
def eps3():
    return EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )

@pytest.fixture
def electrons():
    return VectorType("electrons", ("a", "b", "c", "d"))

@pytest.fixture
def ctx(electrons):
    return Context(types=(electrons,))

@pytest.fixture
def plan(dot, eps3, ctx):
    return Plan(context=ctx, operations=(dot, eps3))

@pytest.fixture
def event_3d():
    return {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "d": np.array([1.0, 1.0, 1.0]) / 3.0**0.5,
    }


# ---------------------------------------------------------------------------
# Basic shape and type checks
# ---------------------------------------------------------------------------

_SKIP_NO_REGISTRY = pytest.mark.skip(reason="Pre-registry path removed; will be replaced by registry-based tests")


@_SKIP_NO_REGISTRY
def test_encode_returns_ndarray(dot, plan, ctx, event_3d):
    result = encode(plan, event_3d)
    assert isinstance(result, np.ndarray)

@_SKIP_NO_REGISTRY
def test_encode_returns_real(dot, plan, ctx, event_3d):
    result = encode(plan, event_3d)
    assert result.dtype == np.float64

@_SKIP_NO_REGISTRY
def test_encode_length_matches_describe(dot, eps3, plan, ctx, event_3d):
    """encode() output length matches sum of segment lengths from describe_encoding()."""
    result = encode(plan, event_3d)
    assert len(result) == sum(s.length for s in describe_encoding(plan))

@_SKIP_NO_REGISTRY
def test_encode_dot_only_length(dot, ctx, event_3d):
    """Length check for dot-only plan matches describe_encoding()."""
    plan = Plan(context=ctx, operations=(dot,))
    result = encode(plan, event_3d)
    assert len(result) == sum(s.length for s in describe_encoding(plan))

def test_encode_brute_length_matches_sum_of_counts(dot, eps3, plan, ctx, event_3d):
    """encode_brute() gives 2*sum(pf.count()) reals — n complex coeffs × 2 reals each."""
    fo_list = repS(ctx, (dot, eps3))
    type_sizes = tuple(g.size for g in ctx.types)
    expected_len = 2 * sum(
        pf.count(type_sizes)
        for pf in canonical_pair_flavours(fo_list, ctx)
    )
    assert len(encode_brute(plan, event_3d)) == expected_len


# ---------------------------------------------------------------------------
# Permutation invariance
# ---------------------------------------------------------------------------

@_SKIP_NO_REGISTRY
def test_encode_invariant_under_label_permutation(dot, ctx):
    """
    Swapping two particle labels must not change the encoding.
    encode(E) == encode(E') where E' is E with labels a and b swapped.
    """
    plan = Plan(context=ctx, operations=(dot,))
    event = {
        "a": np.array([1.0, 2.0, 3.0]),
        "b": np.array([4.0, 5.0, 6.0]),
        "c": np.array([7.0, 0.0, 1.0]),
        "d": np.array([0.0, 3.0, 2.0]),
    }
    event_swapped = dict(event)
    event_swapped["a"], event_swapped["b"] = event["b"], event["a"]

    result       = encode(plan, event)
    result_swapped = encode(plan, event_swapped)
    np.testing.assert_array_almost_equal(result, result_swapped)

@_SKIP_NO_REGISTRY
def test_encode_invariant_under_cyclic_permutation(dot, ctx):
    """Cyclic permutation a→b→c→d→a leaves the encoding unchanged."""
    plan = Plan(context=ctx, operations=(dot,))
    event = {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "d": np.array([1.0, 1.0, 0.0]) / 2.0**0.5,
    }
    event_cycled = {
        "a": event["b"],
        "b": event["c"],
        "c": event["d"],
        "d": event["a"],
    }
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_cycled),
    )


# ---------------------------------------------------------------------------
# Determinism and reproducibility
# ---------------------------------------------------------------------------

@_SKIP_NO_REGISTRY
def test_encode_deterministic(dot, plan, ctx, event_3d):
    """encode() called twice on the same event returns the same vector."""
    r1 = encode(plan, event_3d)
    r2 = encode(plan, event_3d)
    np.testing.assert_array_equal(r1, r2)


# ---------------------------------------------------------------------------
# Two-group context
# ---------------------------------------------------------------------------

@_SKIP_NO_REGISTRY
def test_encode_two_groups(dot, eps3):
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(dot, eps3))
    event = {
        "a": np.array([1.0, 0.0, 0.0]),
        "b": np.array([0.0, 1.0, 0.0]),
        "c": np.array([0.0, 0.0, 1.0]),
        "p": np.array([1.0, 1.0, 0.0]) / 2.0**0.5,
        "q": np.array([0.0, 1.0, 1.0]) / 2.0**0.5,
    }
    result = encode(plan, event)
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    assert len(result) > 0

@_SKIP_NO_REGISTRY
def test_encode_two_groups_permutation_invariant(dot):
    """Swapping muon labels p↔q leaves the encoding unchanged."""
    electrons = VectorType("electrons", ("a", "b"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(dot,))
    event = {
        "a": np.array([1.0, 2.0, 0.0]),
        "b": np.array([3.0, 0.0, 1.0]),
        "p": np.array([0.0, 1.0, 4.0]),
        "q": np.array([2.0, 0.0, 1.0]),
    }
    event_swapped = dict(event)
    event_swapped["p"], event_swapped["q"] = event["q"], event["p"]
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_swapped),
    )


# ---------------------------------------------------------------------------
# Optimization 5a: compressed polynomial embedding
# ---------------------------------------------------------------------------

def _make_ops_and_event():
    """Shared setup: four 3-D vectors, all four symmetry-class operation pairs."""
    labels = ("a", "b", "c", "d")
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    event = {
        "a": np.array([1.0, 0.2, 0.3]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.3, 1.0]),
        "d": np.array([0.2, 0.7, 0.1]),
    }
    return dot, eps3, labels, event


def test_eval_pair_orbit_positive_count(dot, eps3, ctx):
    """eval_pair_orbit_positive always returns exactly pf.count() values."""
    fo_list = repS(ctx, (dot, eps3))
    type_sizes = tuple(g.size for g in ctx.types)
    for pf in canonical_pair_flavours(fo_list, ctx):
        plan = Plan(context=ctx, operations=(dot, eps3))
        event = {l: np.random.randn(3) for l in ctx.all_labels}
        pos = eval_pair_orbit_positive(pf, plan, event)
        assert len(pos) == pf.count(type_sizes), (
            f"Expected {pf.count(type_sizes)}, got {len(pos)} for {pf!r}"
        )


def test_eval_pair_orbit_positive_count_unusual_label_order():
    """Count invariant holds even with unusual alphabetical label ordering."""
    # "toast" < "apple" < "zebra" under default str ordering is False,
    # but "apple" < "toast" < "zebra" — the point is the Atom constructor
    # sorts internally and absorbs parity into sign, so count is always n.
    labels = ("apple", "toast", "zebra")
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    g = VectorType("particles", labels)
    ctx = Context((g,))
    plan = Plan(context=ctx, operations=(dot, eps3))
    event = {l: np.random.randn(3) for l in labels}
    fo_list = repS(ctx, (dot, eps3))
    type_sizes = tuple(g.size for g in ctx.types)
    for pf in canonical_pair_flavours(fo_list, ctx):
        pos = eval_pair_orbit_positive(pf, plan, event)
        assert len(pos) == pf.count(type_sizes), (
            f"Expected {pf.count(type_sizes)}, got {len(pos)} for {pf!r}"
        )


@pytest.mark.parametrize("perm", [
    {"a": "b", "b": "a", "c": "c", "d": "d"},
    {"a": "b", "b": "c", "c": "d", "d": "a"},
    {"a": "d", "b": "c", "c": "b", "d": "a"},
])
@_SKIP_NO_REGISTRY
def test_compressed_encoding_permutation_invariant_dot(perm, ctx):
    """Compressed encoding of dot-only plan is invariant under label permutation."""
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    plan = Plan(context=ctx, operations=(dot,))
    event = {
        "a": np.array([1.0, 0.2, 0.3]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.3, 1.0]),
        "d": np.array([0.2, 0.7, 0.1]),
    }
    event_permuted = {perm[k]: v for k, v in event.items()}
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_permuted),
    )


@pytest.mark.parametrize("perm", [
    {"a": "b", "b": "a", "c": "c", "d": "d"},
    {"a": "b", "b": "c", "c": "d", "d": "a"},
])
@_SKIP_NO_REGISTRY
def test_compressed_encoding_permutation_invariant_eps3(perm, ctx):
    # TODO Make this a much better test that looks at cases that can't be comressed too, e.e. eps3(a,p,v)
    """Compressed encoding with ANTISYMMETRIC eps3 is permutation-invariant."""
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    plan = Plan(context=ctx, operations=(eps3,))
    event = {
        "a": np.array([1.0, 0.2, 0.3]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.3, 1.0]),
        "d": np.array([0.2, 0.7, 0.1]),
    }
    event_permuted = {perm[k]: v for k, v in event.items()}
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_permuted),
    )


def test_compressed_length_sym_sym(dot, ctx, event_3d):
    """SYM×SYM: compressed length == orbit_size (no compression)."""
    plan = Plan(context=ctx, operations=(dot,))
    fo_list = repS(ctx, (dot,))
    type_sizes = tuple(g.size for g in ctx.types)
    for pf in canonical_pair_flavours(fo_list, ctx):
        assert pf.count(type_sizes) == pf.orbit_size(type_sizes)


def test_compressed_length_antisym_antisym(ctx, event_3d):
    # TODO FIXME : This is wrong, referring to ANTISYMM etc rather than TYPE_11,12,21,22
    """ANTISYM×ANTISYM: eval_pair_orbit_positive returns exactly count elements.

    For AA PairFlavours, DirectOrbitEnumerator generates all four sign-combinations
    per label assignment: (+u,+v), (-u,+v), (+u,-v), (-u,-v).  eval_pair_orbit_positive
    keeps only the (+,+) variant, so its length equals pf.count.  This is the property
    the _embed_compressed encoder relies on for AA pairs.

    Note: orbit_size is now the true G-orbit size (|G|/|Stab|), which equals 2*count
    for the non-zero-overlap AA case that arises in the 4-electron context.  The old
    assertion pf.count == orbit_size // 4 was based on an incorrect orbit_size formula
    and has been replaced with the encoding-relevant invariant here.
    """
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    plan = Plan(context=ctx, operations=(eps3,))
    fo_list = repS(ctx, (eps3,))
    type_sizes = tuple(g.size for g in ctx.types)
    for pf in canonical_pair_flavours(fo_list, ctx):
        pos_values = eval_pair_orbit_positive(pf, plan, event_3d)
        assert len(pos_values) == pf.count(type_sizes)


def test_compressed_length_sym_antisym(ctx):
    # TODO FIXME : This is wrong, referring to ANTISYMM etc rather than TYPE_11,12,21,22
    """SYM×ANTISYM or ANTISYM×SYM: compressed length == orbit_size / 2."""
    dot = EvaluableOperation(
        name="dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    eps3 = EvaluableOperation(
        name="eps3", rank=3, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], np.cross(vecs[1], vecs[2]))),
    )
    plan = Plan(context=ctx, operations=(dot, eps3))
    fo_list = repS(ctx, (dot, eps3))
    type_sizes = tuple(g.size for g in ctx.types)
    mixed_pfs = [
        pf for pf in canonical_pair_flavours(fo_list, ctx)
        if (pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
        != (pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
    ]
    assert len(mixed_pfs) > 0, "Need at least one SYM×ANTISYM pair for this test"
    for pf in mixed_pfs:
        assert pf.count(type_sizes) == pf.orbit_size(type_sizes) // 2


# ---------------------------------------------------------------------------
# TYPE_NEG: correlated-negation embedding in a 2-group context
# ---------------------------------------------------------------------------

def _make_eps2():
    return EvaluableOperation(
        name="eps2", rank=2, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda vecs: float(vecs[0][0]*vecs[1][1] - vecs[0][1]*vecs[1][0]),
    )


def test_embed_compressed_type_neg_permutation_invariant():
    # TODO: I don't like TYPE_NEG intruding on TYPES 11,12,21,22. That's not what they were for.
    # But TYPE_NEG is not all bad. It does represent a case where a special form of compression
    # can be used. It is just measuring a different thing to the 12/21/12/22 TYPES.
    """
    TYPE_NEG pair (AA, full per-group overlap in electron group): the encoding
    must be invariant under label permutations that negate BOTH atom signs.

    eps2[(2,0)] × eps2[(2,0)] with overlap=(2,0) in a 2-group context
    (electrons={a,b,c}, muons={p,q}).  Swapping labels a↔b negates both
    atom eval values simultaneously (TYPE_NEG orbit = {z_k, -z_k}).
    """
    from symcoder.encode import _embed_compressed
    from symcoder.pairs import eval_pair_orbit_positive
    eps2 = _make_eps2()
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(eps2,))

    np.random.seed(42)
    event = {
        "a": np.array([1.0, 0.3, 0.2]),
        "b": np.array([0.1, 1.0, 0.4]),
        "c": np.array([0.5, 0.2, 1.0]),
        "p": np.array([0.7, 0.1, 0.3]),
        "q": np.array([0.2, 0.9, 0.1]),
    }
    event_swapped = dict(event)
    event_swapped["a"], event_swapped["b"] = event["b"], event["a"]

    fo_list = repS(ctx, (eps2,))
    type_sizes = tuple(g.size for g in ctx.types)

    # Find a TYPE_NEG PairFlavour (AA, full electron overlap)
    from symcoder.encode import _sign_correlation_type_from_pf
    from symatom.group import SignCorrelationType
    type_neg_found = False
    for pf in canonical_pair_flavours(fo_list, ctx):
        sct = _sign_correlation_type_from_pf(pf)
        z_pos         = np.array(eval_pair_orbit_positive(pf, plan, event),         dtype=complex)
        z_pos_swapped = np.array(eval_pair_orbit_positive(pf, plan, event_swapped), dtype=complex)
        c1 = _embed_compressed(z_pos, pf)
        c2 = _embed_compressed(z_pos_swapped, pf)
        np.testing.assert_array_almost_equal(
            c1, c2,
            err_msg=f"_embed_compressed not invariant under a↔b swap for {pf!r} (type={sct})"
        )
        if sct == SignCorrelationType.TYPE_NEG:
            type_neg_found = True
    assert type_neg_found, "No TYPE_NEG PairFlavour found — test setup may be wrong"


def test_embed_compressed_type_neg_distinguishes_events():
    """
    TYPE_NEG embedding must be injective: two events with genuinely different
    orbit values must produce different encoded vectors.

    This test checks that the TYPE_NEG branch (embed z_pos²) preserves
    Im(z²) = 2·Re(z)·Im(z), which would be lost if we incorrectly used
    the TYPE_22 branch (embed {z², conj(z²)} → real polynomial → loses Im(z²)).
    """
    from symcoder.encode import _embed_compressed, _sign_correlation_type_from_pf
    from symatom.group import SignCorrelationType
    eps2 = _make_eps2()
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(eps2,))

    # Two events designed to yield z_k values with opposite Im(z_k²)
    # eps2(a,b) = a[0]*b[1] - a[1]*b[0] (the 2D cross product component)
    # Construct so that Re(z) = eps2 eval, Im(z) = another eps2 eval, with Im(z²) different
    event1 = {
        "a": np.array([2.0, 0.0]), "b": np.array([0.0, 1.0]),  # eps2(a,b) = 2
        "c": np.array([0.0, 0.0]),
        "p": np.array([1.0, 0.0]), "q": np.array([0.0, 1.0]),  # eps2(a,p) context
    }
    event2 = {
        "a": np.array([2.0, 0.0]), "b": np.array([0.0, -1.0]),  # eps2(a,b) = -2
        "c": np.array([0.0, 0.0]),
        "p": np.array([1.0, 0.0]), "q": np.array([0.0, 1.0]),
    }

    fo_list = repS(ctx, (eps2,))
    type_sizes = tuple(g.size for g in ctx.types)

    type_neg_pfs = [
        pf for pf in canonical_pair_flavours(fo_list, ctx)
        if _sign_correlation_type_from_pf(pf) == SignCorrelationType.TYPE_NEG
    ]
    assert len(type_neg_pfs) > 0, "No TYPE_NEG PairFlavour found — test setup may be wrong"

    for pf in type_neg_pfs:
        z1 = np.array(eval_pair_orbit_positive(pf, plan, event1), dtype=complex)
        z2 = np.array(eval_pair_orbit_positive(pf, plan, event2), dtype=complex)
        # Only test pairs where the z values actually differ
        if np.allclose(z1, z2) or np.allclose(z1**2, z2**2):
            continue
        c1 = _embed_compressed(z1, pf)
        c2 = _embed_compressed(z2, pf)
        # The embeddings must differ (TYPE_NEG should preserve z² information)
        assert not np.allclose(c1, c2), (
            f"TYPE_NEG embedding failed to distinguish events for {pf!r}: "
            f"z1={z1}, z2={z2}, z1²={z1**2}, z2²={z2**2}"
        )


@_SKIP_NO_REGISTRY
def test_encode_two_groups_permutation_invariant_eps2():
    """encode() with an ANTISYMMETRIC rank-2 op in a 2-group context is permutation-invariant.

    This exercises the TYPE_NEG branch (AA pairs with full electron overlap) and
    ensures the full encode pipeline handles it correctly.
    """
    eps2 = _make_eps2()
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx  = Context(types=(electrons, muons))
    plan = Plan(context=ctx, operations=(eps2,))
    event = {
        "a": np.array([1.3, 0.2]),
        "b": np.array([0.7, 1.1]),
        "c": np.array([0.4, 0.9]),
        "p": np.array([0.8, 0.3]),
        "q": np.array([0.1, 1.4]),
    }
    # Swap electrons a↔b: a TYPE_NEG orbit has z_k → -z_k, so encoding must be stable
    event_swapped = dict(event)
    event_swapped["a"], event_swapped["b"] = event["b"], event["a"]
    np.testing.assert_array_almost_equal(
        encode(plan, event),
        encode(plan, event_swapped),
        err_msg="encode() not invariant under a↔b swap in 2-group eps2 context"
    )
