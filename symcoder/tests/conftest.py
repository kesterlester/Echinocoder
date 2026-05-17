"""Shared pytest fixtures for symcoder tests."""
import pytest
from symcoder.encoders import (
    OrbitEncoderFactory,
    SortEncoderFactory,
    standard_row_pair_factories,
    OverlapBlockEncoderFactory,
    Phase2EncoderFactory,
)


@pytest.fixture
def orbit_factory():
    return OrbitEncoderFactory([SortEncoderFactory()])


@pytest.fixture
def phase2_factory():
    return Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories())
    ])
