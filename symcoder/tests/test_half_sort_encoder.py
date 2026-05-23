"""
Tests for HalfSortEncoder and HalfSortEncoderFactory.

Concrete orbit-size facts (verified empirically):
  eps3 on 3-particle group: orbit_size=2, half_dim=1, sort_dim=2
  eps3 on 4-particle group: orbit_size=8, half_dim=4, sort_dim=8
  dot  on 3-particle group: orbit_size=3, half_dim=N/A (not pairable)
  dot  on 4-particle group: orbit_size=6, half_dim=N/A (not pairable)
"""
from __future__ import annotations

import numpy as np
import pytest

from symatom import repS, Context, Plan, VectorType, ArgumentSymmetry
from symcoder import EvaluableOperation
from symcoder.encoders import (
    HalfSortEncoder,
    HalfSortEncoderFactory,
    OrbitSpec,
    OrbitEncoderFactory,
    SortEncoder,
    SortEncoderFactory,
)


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def eps3():
    return EvaluableOperation(
        "eps3", rank=3, odd_parity=True,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        eval_fn=lambda v: float(np.dot(v[0], np.cross(v[1], v[2]))),
    )


@pytest.fixture
def dot():
    return EvaluableOperation(
        "dot", rank=2, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda v: float(np.dot(v[0], v[1])),
    )


@pytest.fixture
def mag():
    return EvaluableOperation(
        "mag", rank=1, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))),
    )


@pytest.fixture
def ctx3():
    return Context((VectorType("g", ("a", "b", "c")),))


@pytest.fixture
def ctx4():
    return Context((VectorType("g", ("a", "b", "c", "d")),))


def _make_plan(ctx, *ops):
    return Plan(context=ctx, operations=ops)


def _first_spec_and_orbit(plan):
    """Return (OrbitSpec, orbit_list) for the first FlavouredOperator in the plan."""
    fo_list = repS(plan.context, plan.operations)
    fo = fo_list[0]
    spec = OrbitSpec.from_flavoured_operator(fo)
    orbit = plan.context.the_group.orbit(fo.canonical_representative())
    return spec, orbit


def _event_3d(labels):
    """Random 3-vector event for a list of labels."""
    rng = np.random.default_rng(42)
    return {l: rng.standard_normal(3) for l in labels}


# ---------------------------------------------------------------------------
# HalfSortEncoderFactory.assess — eligibility
# ---------------------------------------------------------------------------

class TestHalfSortEligibility:

    def test_assess_returns_half_sort_encoder_for_eps3_3p(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        result = HalfSortEncoderFactory().assess(spec, plan)
        assert len(result) == 1
        assert isinstance(result[0], HalfSortEncoder)

    def test_assess_returns_half_sort_encoder_for_eps3_4p(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        result = HalfSortEncoderFactory().assess(spec, plan)
        assert len(result) == 1
        assert isinstance(result[0], HalfSortEncoder)

    def test_assess_returns_empty_for_dot_3p(self, dot, ctx3):
        plan = _make_plan(ctx3, dot)
        spec, _ = _first_spec_and_orbit(plan)
        assert HalfSortEncoderFactory().assess(spec, plan) == []

    def test_assess_returns_empty_for_dot_4p(self, dot, ctx4):
        plan = _make_plan(ctx4, dot)
        spec, _ = _first_spec_and_orbit(plan)
        assert HalfSortEncoderFactory().assess(spec, plan) == []

    def test_assess_returns_empty_for_mag(self, mag, ctx3):
        plan = _make_plan(ctx3, mag)
        spec, _ = _first_spec_and_orbit(plan)
        assert HalfSortEncoderFactory().assess(spec, plan) == []


# ---------------------------------------------------------------------------
# HalfSortEncoder.output_dim — compression ratio
# ---------------------------------------------------------------------------

class TestHalfSortOutputDim:

    def test_output_dim_is_half_orbit_size_eps3_3p(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, orbit = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        assert enc.output_dim == len(orbit) // 2  # orbit_size=2 → dim=1

    def test_output_dim_is_half_orbit_size_eps3_4p(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, orbit = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        assert enc.output_dim == len(orbit) // 2  # orbit_size=8 → dim=4

    def test_output_dim_concrete_eps3_3p(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        assert enc.output_dim == 1

    def test_output_dim_concrete_eps3_4p(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        assert enc.output_dim == 4

    def test_half_sort_output_dim_less_than_sort_output_dim(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        half_enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        sort_enc  = SortEncoderFactory().assess(spec, plan)[0]
        assert half_enc.output_dim < sort_enc.output_dim


# ---------------------------------------------------------------------------
# HalfSortEncoder.method_name and priority
# ---------------------------------------------------------------------------

class TestHalfSortMetadata:

    def test_method_name_is_half_sort(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        assert enc.method_name == "half_sort"

    def test_priority_higher_than_sort_encoder(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        half_enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        sort_enc  = SortEncoderFactory().assess(spec, plan)[0]
        assert half_enc.priority > sort_enc.priority


# ---------------------------------------------------------------------------
# HalfSortEncoder.encode — correctness
# ---------------------------------------------------------------------------

class TestHalfSortEncoding:

    def test_output_length_matches_output_dim(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        event = _event_3d(["a", "b", "c"])
        result = enc.encode(event)
        assert len(result.values) == enc.output_dim

    def test_output_is_float64(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        result = enc.encode(_event_3d(["a", "b", "c"]))
        assert result.values.dtype == np.float64

    def test_output_is_sorted(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        result = enc.encode(_event_3d(["a", "b", "c", "d"]))
        assert np.all(result.values[:-1] <= result.values[1:])

    def test_output_is_non_negative(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        result = enc.encode(_event_3d(["a", "b", "c", "d"]))
        assert np.all(result.values >= 0.0)

    def test_permutation_invariant_eps3_3p(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        rng = np.random.default_rng(7)
        vecs = {l: rng.standard_normal(3) for l in ["a", "b", "c"]}
        out1 = enc.encode(vecs).values
        # Permute labels in the event
        swapped = {"a": vecs["b"], "b": vecs["c"], "c": vecs["a"]}
        out2 = enc.encode(swapped).values
        np.testing.assert_allclose(out1, out2, atol=1e-12)

    def test_permutation_invariant_eps3_4p(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        rng = np.random.default_rng(13)
        labels = ["a", "b", "c", "d"]
        vecs = {l: rng.standard_normal(3) for l in labels}
        out_orig = enc.encode(vecs).values
        # Reverse label order
        reversed_vecs = dict(zip(labels, [vecs[l] for l in reversed(labels)]))
        out_rev = enc.encode(reversed_vecs).values
        np.testing.assert_allclose(out_orig, out_rev, atol=1e-12)

    def test_different_events_give_different_outputs(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        spec, _ = _first_spec_and_orbit(plan)
        enc = HalfSortEncoderFactory().assess(spec, plan)[0]
        rng = np.random.default_rng(99)
        labels = ["a", "b", "c", "d"]
        ev1 = {l: rng.standard_normal(3) for l in labels}
        ev2 = {l: rng.standard_normal(3) for l in labels}
        out1 = enc.encode(ev1).values
        out2 = enc.encode(ev2).values
        assert not np.allclose(out1, out2)


# ---------------------------------------------------------------------------
# SortEncoder as fallback — eps3 absent, dot-only plan
# ---------------------------------------------------------------------------

class TestSortFallback:

    def test_sort_encoder_used_when_half_sort_inapplicable(self, dot, ctx3):
        plan = _make_plan(ctx3, dot)
        spec, _ = _first_spec_and_orbit(plan)
        assert HalfSortEncoderFactory().assess(spec, plan) == []
        result = SortEncoderFactory().assess(spec, plan)
        assert len(result) == 1
        assert isinstance(result[0], SortEncoder)

    def test_orbit_factory_uses_sort_for_dot(self, dot, ctx3):
        plan = _make_plan(ctx3, dot)
        factory = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
        orbit_enc = factory.build(plan)
        segs = orbit_enc.describe()
        assert all(s.method_name == "sort" for s in segs)

    def test_orbit_factory_uses_half_sort_for_eps3(self, eps3, ctx3):
        plan = _make_plan(ctx3, eps3)
        factory = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
        orbit_enc = factory.build(plan)
        segs = orbit_enc.describe()
        assert all(s.method_name == "half_sort" for s in segs)


# ---------------------------------------------------------------------------
# OrbitEncoderFactory — min(output_dim) selection
# ---------------------------------------------------------------------------

class TestOrbitFactorySelection:

    def test_half_sort_wins_over_sort_for_eps3(self, eps3, ctx4):
        plan = _make_plan(ctx4, eps3)
        factory = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
        orbit_enc = factory.build(plan)
        # orbit_size=8 for eps3/4p, half_sort gives dim=4, sort gives dim=8
        assert orbit_enc.output_dim == 4

    def test_sort_wins_when_half_sort_inapplicable(self, dot, ctx4):
        plan = _make_plan(ctx4, dot)
        factory = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
        orbit_enc = factory.build(plan)
        # orbit_size=6 for dot/4p, only sort applies
        assert orbit_enc.output_dim == 6

    def test_mixed_plan_uses_correct_encoder_per_orbit(self, eps3, dot, ctx3):
        plan = _make_plan(ctx3, eps3, dot)
        factory = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
        orbit_enc = factory.build(plan)
        segs = orbit_enc.describe()
        by_op = {s.op_u: s.method_name for s in segs}
        assert by_op["eps3"] == "half_sort"
        assert by_op["dot"] == "sort"
