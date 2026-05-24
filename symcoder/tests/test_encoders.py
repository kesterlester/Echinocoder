"""
Tests for the symcoder.encoders registry infrastructure.

These tests exercise the structural plumbing — OrbitSpec construction,
registry registration, query_all filtering, best_for selection, encode
dispatch — using lightweight mock factories and encoders.  They do NOT
test the concrete SortEncoder / PolyEncoder stubs.
"""
from __future__ import annotations

import numpy as np
import pytest

from symcoder.encoders import (
    AtomOrbitEncoder,
    AtomOrbitEncoderFactory,
    AtomOrbitEncoderRegistry,
    EncodingResult,
    OrbitSpec,
    OrbitSpecForm,
)
from symatom.atoms import Atom, Operation, ArgumentSymmetry


# ---------------------------------------------------------------------------
# Minimal mock encoders and factories
# ---------------------------------------------------------------------------

class _MockEncoder(AtomOrbitEncoder):
    """A ready-to-use encoder with configurable properties."""
    def __init__(self, priority: float = 1.0, output_dim: int = 4,
                 method: str = "mock"):
        self._priority   = priority
        self._output_dim = output_dim
        self._method     = method

    @property
    def output_dim(self) -> int:
        return self._output_dim

    @property
    def priority(self) -> float:
        return self._priority

    @property
    def method_name(self) -> str:
        return self._method

    def encode(self, event: dict) -> EncodingResult:
        return EncodingResult(
            values=np.zeros(self._output_dim, dtype=np.float64),
            metadata={"method": self._method},
        )


class _CapableFactory(AtomOrbitEncoderFactory):
    """Always returns one encoder for any spec."""
    def __init__(self, priority: float = 1.0, output_dim: int = 4,
                 method: str = "mock"):
        self._priority   = priority
        self._output_dim = output_dim
        self._method     = method

    def assess(self, spec, plan) -> list[AtomOrbitEncoder]:
        return [_MockEncoder(self._priority, self._output_dim, self._method)]


class _IncapableFactory(AtomOrbitEncoderFactory):
    """Never returns any encoder."""
    def assess(self, spec, plan) -> list[AtomOrbitEncoder]:
        return []


class _SelectiveFactory(AtomOrbitEncoderFactory):
    """Only capable for FLAVOURED_OPERATOR specs."""
    def assess(self, spec, plan) -> list[AtomOrbitEncoder]:
        if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
            return [_MockEncoder(priority=1.5, output_dim=8, method="selective")]
        return []


class _MultiFactory(AtomOrbitEncoderFactory):
    """Returns two encoders for any spec (ascending and descending variants)."""
    def assess(self, spec, plan) -> list[AtomOrbitEncoder]:
        return [
            _MockEncoder(priority=0.4, output_dim=3, method="asc"),
            _MockEncoder(priority=0.6, output_dim=3, method="desc"),
        ]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _dummy_spec(form: OrbitSpecForm = OrbitSpecForm.EXPLICIT_ORBIT) -> OrbitSpec:
    return OrbitSpec(form=form, payload=[])


# ---------------------------------------------------------------------------
# OrbitSpec construction
# ---------------------------------------------------------------------------

class TestOrbitSpec:

    def test_from_atom_sets_form(self):
        op   = Operation("dot", rank=2, odd_parity=False,
                         argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                         eval_fn=lambda v: 0.0)
        atom = Atom(op, ("a", "b"), sign=1)
        spec = OrbitSpec.from_atom(atom)
        assert spec.form == OrbitSpecForm.REPRESENTATIVE_ATOM
        assert spec.payload is atom

    def test_from_explicit_orbit_materialises_iterable(self):
        op   = Operation("dot", rank=2, odd_parity=False,
                         argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                         eval_fn=lambda v: 0.0)
        atom = Atom(op, ("a", "b"), sign=1)
        spec = OrbitSpec.from_explicit_orbit(iter([atom]))  # iterator, not list
        assert spec.form == OrbitSpecForm.EXPLICIT_ORBIT
        assert isinstance(spec.payload, list)
        assert len(spec.payload) == 1

    def test_from_flavoured_operator(self):
        sentinel = object()
        spec     = OrbitSpec.from_flavoured_operator(sentinel)
        assert spec.form == OrbitSpecForm.FLAVOURED_OPERATOR
        assert spec.payload is sentinel

    def test_repr_is_informative(self):
        spec = _dummy_spec()
        assert "EXPLICIT_ORBIT" in repr(spec)


# ---------------------------------------------------------------------------
# Registry: registration
# ---------------------------------------------------------------------------

class TestRegistration:

    def test_empty_registry_has_zero_length(self):
        assert len(AtomOrbitEncoderRegistry()) == 0

    def test_register_increases_length(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory())
        assert len(reg) == 1
        reg.register(_IncapableFactory())
        assert len(reg) == 2

    def test_iteration_yields_factories_in_registration_order(self):
        reg = AtomOrbitEncoderRegistry()
        a   = _CapableFactory()
        b   = _IncapableFactory()
        reg.register(a)
        reg.register(b)
        listed = list(reg)
        assert listed[0] is a
        assert listed[1] is b


# ---------------------------------------------------------------------------
# Registry: query_all
# ---------------------------------------------------------------------------

class TestQueryAll:

    def test_returns_encoders_from_capable_factories_only(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(method="yes"))
        reg.register(_IncapableFactory())
        results = reg.query_all(_dummy_spec(), plan=None)
        assert len(results) == 1
        assert results[0].method_name == "yes"

    def test_empty_result_when_no_factory_capable(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_IncapableFactory())
        assert reg.query_all(_dummy_spec(), plan=None) == []

    def test_results_in_factory_registration_order(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(priority=1.0, method="first"))
        reg.register(_CapableFactory(priority=2.0, method="second"))
        results = reg.query_all(_dummy_spec(), plan=None)
        assert results[0].method_name == "first"
        assert results[1].method_name == "second"

    def test_multi_factory_returns_all_its_encoders(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_MultiFactory())
        results = reg.query_all(_dummy_spec(), plan=None)
        assert len(results) == 2
        methods = {e.method_name for e in results}
        assert methods == {"asc", "desc"}

    def test_selective_factory_filtered_by_form(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_SelectiveFactory())
        assert reg.query_all(_dummy_spec(OrbitSpecForm.EXPLICIT_ORBIT), None) == []
        results = reg.query_all(_dummy_spec(OrbitSpecForm.FLAVOURED_OPERATOR), None)
        assert len(results) == 1
        assert results[0].method_name == "selective"

    def test_encoder_properties_are_accessible(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(priority=3.7, output_dim=12, method="foo"))
        enc = reg.query_all(_dummy_spec(), None)[0]
        assert enc.output_dim  == 12
        assert enc.priority    == 3.7
        assert enc.method_name == "foo"

    def test_flattened_across_multiple_factories(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_MultiFactory())   # returns 2
        reg.register(_CapableFactory()) # returns 1
        results = reg.query_all(_dummy_spec(), plan=None)
        assert len(results) == 3


# ---------------------------------------------------------------------------
# Registry: best_for
# ---------------------------------------------------------------------------

class TestBestFor:

    def test_returns_none_when_no_factory_capable(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_IncapableFactory())
        assert reg.best_for(_dummy_spec(), plan=None) is None

    def test_returns_highest_priority_encoder(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(priority=0.5, output_dim=2, method="low"))
        reg.register(_CapableFactory(priority=2.0, output_dim=6, method="high"))
        enc = reg.best_for(_dummy_spec(), plan=None)
        assert enc.priority    == 2.0
        assert enc.output_dim  == 6
        assert enc.method_name == "high"

    def test_registration_order_breaks_ties(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(priority=1.0, method="first"))
        reg.register(_CapableFactory(priority=1.0, method="second"))
        enc = reg.best_for(_dummy_spec(), plan=None)
        assert enc.method_name == "first"

    def test_ignores_incapable_factory(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(priority=0.1, method="only"))
        reg.register(_IncapableFactory())
        enc = reg.best_for(_dummy_spec(), plan=None)
        assert enc.method_name == "only"

    def test_picks_best_across_multi_factory(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_MultiFactory())  # asc(0.4) and desc(0.6)
        enc = reg.best_for(_dummy_spec(), plan=None)
        assert enc.method_name == "desc"  # higher priority within MultiFactory


# ---------------------------------------------------------------------------
# Registry: encode (convenience dispatch)
# ---------------------------------------------------------------------------

class TestEncode:

    def test_raises_when_no_factory_capable(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_IncapableFactory())
        with pytest.raises(RuntimeError, match="No registered factory"):
            reg.encode(_dummy_spec(), event={}, plan=None)

    def test_dispatches_to_best_encoder(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(priority=0.5, output_dim=2))
        reg.register(_CapableFactory(priority=2.0, output_dim=6))
        result = reg.encode(_dummy_spec(), event={}, plan=None)
        assert result.values.shape == (6,)  # from the high-priority encoder

    def test_result_values_are_float64(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory(output_dim=4))
        result = reg.encode(_dummy_spec(), event={}, plan=None)
        assert result.values.dtype == np.float64

    def test_result_metadata_is_dict(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableFactory())
        result = reg.encode(_dummy_spec(), event={}, plan=None)
        assert isinstance(result.metadata, dict)


# ---------------------------------------------------------------------------
# EncodingResult defaults
# ---------------------------------------------------------------------------

class TestEncodingResultDefaults:

    def test_metadata_defaults_to_empty_dict(self):
        res = EncodingResult(values=np.array([1.0, 2.0]))
        assert res.metadata == {}

    def test_values_preserved(self):
        arr = np.array([3.0, 1.0, 2.0], dtype=np.float64)
        res = EncodingResult(values=arr)
        assert np.array_equal(res.values, arr)
