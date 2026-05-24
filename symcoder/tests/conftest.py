"""Shared pytest fixtures for symcoder tests."""
import numpy as np
import pytest
from symatom import ArgumentSymmetry, Operation, VectorType, Context
from symcoder.encoders import (
    OrbitEncoderFactory,
    SortEncoderFactory,
    HalfSortEncoderFactory,
    standard_row_pair_factories,
    OverlapBlockEncoderFactory,
    Phase2EncoderFactory,
)

# ---------------------------------------------------------------------------
# Integer-mode vectors
# ---------------------------------------------------------------------------
# Pre-generated integer vectors for the integer-arithmetic test mode.
# Values in [-4, 4]: mixed signs and zeros keep polynomial coefficient
# magnitudes small.  Zeros are frequent in real-world inputs and must
# be handled correctly — they are NOT edge cases.
_INT_VECS: np.ndarray = (
    np.random.default_rng(seed=0).integers(-4, 5, size=(20, 3)).astype(float)
)


# ---------------------------------------------------------------------------
# Event-mode parametrisation
# ---------------------------------------------------------------------------

@pytest.fixture(params=["float", "integer"])
def event_mode(request):
    """Parametrise every test over float and integer event modes.

    float:
        Normally-distributed random vectors (N(0,1), seed 42).
        atol = 1e-10.
    integer:
        Exact-integer vectors from [-4, 4] (seed 0, pre-generated).
        Zeros are included — they are frequent in real-world inputs.
        atol = 1e-6 for unpolished outputs (double-root splitting from
        Phase 1 collisions gives ~2e-7 errors in root-finding).
        atol_exact = 0.0 for polished outputs (bitwise exact).
    """
    return request.param


# ---------------------------------------------------------------------------
# Context, operations, event
# ---------------------------------------------------------------------------

@pytest.fixture
def ctx():
    return Context((
        VectorType('electrons', ('a', 'b', 'c')),
        VectorType('muons',     ('p', 'q')),
    ))


@pytest.fixture
def ops():
    mag  = Operation('mag',  rank=1, odd_parity=False,
                              argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                              eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))))
    dot  = Operation('dot',  rank=2, odd_parity=False,
                              argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                              eval_fn=lambda v: float(np.dot(v[0], v[1])))
    eps3 = Operation('eps3', rank=3, odd_parity=True,
                              argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                              eval_fn=lambda v: float(np.dot(v[0], np.cross(v[1], v[2]))))
    return mag, dot, eps3


@pytest.fixture
def event(ctx, event_mode):
    """Test event: label → vector mapping.

    float mode:   normally-distributed random vectors (seed 42).
    integer mode: small exact-integer vectors from [-4, 4] (see _INT_VECS).
                  dot/eps3 evaluations are exact integers; polished decoded
                  values match Phase 1 values bitwise.  Zeros are included.
    """
    if event_mode == "integer":
        labels = list(ctx.all_labels)
        return {lbl: _INT_VECS[i] for i, lbl in enumerate(labels)}
    np.random.seed(42)
    return {lbl: np.random.randn(3) for lbl in ctx.all_labels}


# ---------------------------------------------------------------------------
# Tolerance fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def atol(event_mode):
    """Absolute tolerance for unpolished (leaf / block) decoder tests.

    Float mode:   1e-10.
    Integer mode: 1e-6 — double-root splitting from Phase 1 collisions
                  (two atom-pairs evaluating to the same value) reaches
                  ~2e-7 in unpolished root-finding.
    """
    return 1e-6 if event_mode == "integer" else 1e-10


@pytest.fixture
def atol_exact(event_mode):
    """Absolute tolerance for polished-output tests (alignment decoder etc.).

    Float mode:   1e-10.
    Integer mode: 0.0 — polishing replaces noisy roots with exact Phase 1
                  values, so the comparison is bitwise exact.
    """
    return 0.0 if event_mode == "integer" else 1e-10


# ---------------------------------------------------------------------------
# Encoder factories
# ---------------------------------------------------------------------------

@pytest.fixture
def orbit_factory():
    return OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])


@pytest.fixture
def phase2_factory():
    return Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories())
    ])


@pytest.fixture
def phase2_factory_no_drop():
    """Like phase2_factory but with complementarity drop disabled.

    Needed for encoder types (Type22, NegPairEncoder) that are always
    comp-dropped in the standard plan — they contribute no output under the
    complementarity drop rule, so they can only be roundtrip-tested when the
    drop is switched off.
    """
    return Phase2EncoderFactory([
        OverlapBlockEncoderFactory(
            standard_row_pair_factories(), use_complementarity_drop=False
        )
    ])
