"""
Tests for σ-compressed self-pairing encoders.

SelfPairType11Encoder (method "11_sigma") and SelfPairNegEncoder (method "neg_sigma")
exploit the swap-symmetry of self-pairing blocks (op_u==op_v, F_u==F_v) to halve
the encoding length from 2n to n reals.

The tests verify:
  1. TYPE_NEG self-pairing reps are always form-1 (+u, +v)
  2. output_dim is n (not 2n) for both σ-encoders
  3. Permutation invariance for both
  4. Injectivity (the "bold claim": σ-compression loses no information)
  5. Total encoding is shorter than without σ-factories
"""
import pytest
import numpy as np
from symatom import ArgumentSymmetry, VectorType, Context, Plan
from symcoder import EvaluableOperation, encode, describe_encoding
from symcoder.encoders import (
    OrbitEncoderFactory, SortEncoderFactory, HalfSortEncoderFactory,
    standard_row_pair_factories, OverlapBlockEncoderFactory, Phase2EncoderFactory,
)
from symcoder.encoders.row_pair_encoders import (
    _pair_orbit_atoms, _try_neg_partition, _is_self_pairing_block,
    SelfPairType11Encoder, SelfPairNegEncoder,
)


# ---------------------------------------------------------------------------
# Fixtures
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
def ctx():
    return Context(types=(VectorType("electrons", ("a", "b", "c", "d")),))

@pytest.fixture
def orbit_factory():
    return OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])

@pytest.fixture
def phase2_factory():
    return Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories())
    ])

@pytest.fixture
def phase2_factory_no_sigma():
    """Phase2 factory with σ-factories removed, for comparison."""
    from symcoder.encoders.row_pair_encoders import (
        SelfPairEncoderFactory, NegPairEncoderFactory, Type22PairEncoderFactory,
        Type21PairEncoderFactory, Type12PairEncoderFactory, Type11PairEncoderFactory,
    )
    no_sigma = [
        SelfPairEncoderFactory(),
        NegPairEncoderFactory(),
        Type22PairEncoderFactory(),
        Type21PairEncoderFactory(),
        Type12PairEncoderFactory(),
        Type11PairEncoderFactory(),
    ]
    return Phase2EncoderFactory([OverlapBlockEncoderFactory(no_sigma)])

def _random_event(ctx, seed):
    rng = np.random.default_rng(seed)
    return {label: rng.standard_normal(3) for label in ctx.all_labels}

def _all_assoc_segs(plan, orbit_factory, phase2_factory):
    return [s for s in describe_encoding(plan, orbit_factory, phase2_factory)
            if s.kind == "ASSOC"]


# ---------------------------------------------------------------------------
# 1. TYPE_NEG self-pairing reps: form-1 and form-2 both occur; encoding works for both
# ---------------------------------------------------------------------------

def test_neg_self_pairing_reps_u_sign_is_positive(eps3, ctx):
    """_try_neg_partition always selects u.sign=+1, but v.sign can be ±1 (form-1 or form-2).

    Form-1: {(+u,+v),(-u,-v)} → rep (+u,+v), v.sign=+1.
    Form-2: {(+u,-v),(-u,+v)} → rep (+u,-v), v.sign=-1.

    Both forms arise in practice.  The encoding uses z² which is τ-closed for
    both forms, so SelfPairNegEncoder is correct in either case.
    """
    plan = Plan(context=ctx, operations=(eps3,))
    from symatom import repS, canonical_pair_flavours
    fo_list = repS(plan.context, plan.operations)
    pfs = canonical_pair_flavours(fo_list, plan.context)
    self_pair_block_pfs = [pf for pf in pfs if _is_self_pairing_block(pf)]
    assert self_pair_block_pfs, "Expected at least one self-pairing block for eps3"
    for pf in self_pair_block_pfs:
        pairs = _pair_orbit_atoms(pf, plan)
        reps = _try_neg_partition(pairs)
        if reps is not None:
            for u, v in reps:
                assert u.sign == +1, f"TYPE_NEG self-pairing rep has u.sign={u.sign} (must be +1)"
                # v.sign can be +1 (form-1) or -1 (form-2) — both are valid


# ---------------------------------------------------------------------------
# 2. output_dim checks
# ---------------------------------------------------------------------------

def test_sigma_type11_output_dim_is_n(dot, ctx, orbit_factory, phase2_factory):
    """SelfPairType11Encoder ASSOC segments have length n, not 2n."""
    plan = Plan(context=ctx, operations=(dot,))
    assoc_segs = _all_assoc_segs(plan, orbit_factory, phase2_factory)
    for s in assoc_segs:
        assert s.symmetry_class == "11_sigma", f"Expected 11_sigma, got {s.symmetry_class}"
        # notional_length = 2 * output_dim (NOTIONAL_FACTOR=2, notional = 2n, actual = n)
        assert s.notional_length == 2 * s.length, (
            f"Expected notional=2*length for 11_sigma, got length={s.length}, notional={s.notional_length}"
        )


def test_sigma_neg_output_dim_is_n(eps3, ctx, orbit_factory, phase2_factory):
    """SelfPairNegEncoder ASSOC segments have length n, not 2n."""
    plan = Plan(context=ctx, operations=(eps3,))
    assoc_segs = _all_assoc_segs(plan, orbit_factory, phase2_factory)
    for s in assoc_segs:
        assert s.symmetry_class == "neg_sigma", f"Expected neg_sigma, got {s.symmetry_class}"
        assert s.notional_length == 4 * s.length, (
            f"Expected notional=4*length for neg_sigma, got length={s.length}, notional={s.notional_length}"
        )


# ---------------------------------------------------------------------------
# 3. Permutation invariance
# ---------------------------------------------------------------------------

def test_sigma_type11_permutation_invariant(dot, ctx, orbit_factory, phase2_factory):
    """dot-only encoding is invariant under relabelling electrons."""
    plan = Plan(context=ctx, operations=(dot,))
    # Two events related by swapping labels a↔c
    labels = list(ctx.all_labels)
    rng = np.random.default_rng(0)
    vecs = {l: rng.standard_normal(3) for l in labels}
    event1 = dict(vecs)
    event2 = dict(vecs)
    event2["a"], event2["c"] = vecs["c"], vecs["a"]
    enc1 = encode(plan, event1, orbit_factory, phase2_factory)
    enc2 = encode(plan, event2, orbit_factory, phase2_factory)
    np.testing.assert_allclose(enc1, enc2, atol=1e-12,
                               err_msg="dot σ-encoding not invariant under a↔c swap")


def test_sigma_neg_permutation_invariant(eps3, ctx, orbit_factory, phase2_factory):
    """eps3-only encoding is invariant under relabelling electrons."""
    plan = Plan(context=ctx, operations=(eps3,))
    labels = list(ctx.all_labels)
    rng = np.random.default_rng(1)
    vecs = {l: rng.standard_normal(3) for l in labels}
    event1 = dict(vecs)
    event2 = dict(vecs)
    event2["a"], event2["c"] = vecs["c"], vecs["a"]
    enc1 = encode(plan, event1, orbit_factory, phase2_factory)
    enc2 = encode(plan, event2, orbit_factory, phase2_factory)
    np.testing.assert_allclose(enc1, enc2, atol=1e-12,
                               err_msg="eps3 σ-encoding not invariant under a↔c swap")


# ---------------------------------------------------------------------------
# 4. Injectivity (bold claim: σ-compression is lossless)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(100))
def test_sigma_type11_injectivity(dot, ctx, orbit_factory, phase2_factory, seed):
    """Two distinct random events give distinct dot-only encodings (100 trials)."""
    plan = Plan(context=ctx, operations=(dot,))
    event_a = _random_event(ctx, seed * 2)
    event_b = _random_event(ctx, seed * 2 + 1)
    enc_a = encode(plan, event_a, orbit_factory, phase2_factory)
    enc_b = encode(plan, event_b, orbit_factory, phase2_factory)
    # Encodings should differ (collision would indicate information loss)
    assert not np.allclose(enc_a, enc_b, atol=1e-10), (
        f"Seed {seed}: two distinct events gave equal dot σ-encodings"
    )


@pytest.mark.parametrize("seed", range(100))
def test_sigma_neg_injectivity(eps3, ctx, orbit_factory, phase2_factory, seed):
    """Two distinct random events give distinct eps3-only encodings (100 trials)."""
    plan = Plan(context=ctx, operations=(eps3,))
    event_a = _random_event(ctx, seed * 2)
    event_b = _random_event(ctx, seed * 2 + 1)
    enc_a = encode(plan, event_a, orbit_factory, phase2_factory)
    enc_b = encode(plan, event_b, orbit_factory, phase2_factory)
    assert not np.allclose(enc_a, enc_b, atol=1e-10), (
        f"Seed {seed}: two distinct events gave equal eps3 σ-encodings"
    )


# ---------------------------------------------------------------------------
# 5. σ-encoders reduce total encoding length vs non-σ factories
# ---------------------------------------------------------------------------

def test_sigma_output_dim_smaller_than_non_sigma(dot, eps3, ctx, orbit_factory,
                                                  phase2_factory, phase2_factory_no_sigma):
    """Total encoding length is shorter with σ-factories than without."""
    plan = Plan(context=ctx, operations=(dot, eps3))
    segs_sigma    = describe_encoding(plan, orbit_factory, phase2_factory)
    segs_no_sigma = describe_encoding(plan, orbit_factory, phase2_factory_no_sigma)
    len_sigma    = sum(s.length for s in segs_sigma)
    len_no_sigma = sum(s.length for s in segs_no_sigma)
    assert len_sigma < len_no_sigma, (
        f"σ-encoding ({len_sigma}) should be shorter than non-σ ({len_no_sigma})"
    )
