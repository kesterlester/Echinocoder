"""Shared pytest fixtures for symcoder tests."""
import pytest
from symcoder.encoders import (
    AtomOrbitEncoderRegistry,
    SortEncoderFactory,
    standard_row_pair_factories,
    OverlapBlockEncoderFactory,
    Phase2EncoderFactory,
)


@pytest.fixture
def registry():
    r = AtomOrbitEncoderRegistry()
    r.register(SortEncoderFactory())
    return r


@pytest.fixture
def phase2_factory():
    return Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories())
    ])
