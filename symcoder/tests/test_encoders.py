"""
Tests for the symcoder.encoders registry infrastructure.

These tests exercise the structural plumbing — OrbitSpec construction,
registry registration, query_all filtering, best_for selection, encode
dispatch — using lightweight mock encoders.  They do NOT test the concrete
SortEncoder / PolyEncoder stubs (those raise NotImplementedError until
implemented).
"""
from __future__ import annotations

import numpy as np
import pytest

from symcoder.encoders import (
    AtomOrbitEncoder,
    AtomOrbitEncoderRegistry,
    EncodingCapability,
    EncodingResult,
    OrbitSpec,
    OrbitSpecForm,
)
from symatom.atoms import Atom, Operation, ArgumentSymmetry


# ---------------------------------------------------------------------------
# Minimal mock encoders
# ---------------------------------------------------------------------------

class _CapableEncoder(AtomOrbitEncoder):
    """Always reports can_encode=True with configurable priority and output_dim."""
    def __init__(self, priority: float = 1.0, output_dim: int = 4,
                 method: str = "mock"):
        self._priority = priority
        self._output_dim = output_dim
        self._method = method

    def assess(self, spec, plan):
        return EncodingCapability(
            can_encode=True,
            output_dim=self._output_dim,
            method_name=self._method,
            priority=self._priority,
        )

    def encode(self, spec, event, plan):
        return EncodingResult(
            values=np.zeros(self._output_dim, dtype=np.float64),
            metadata={"method": self._method},
        )


class _IncapableEncoder(AtomOrbitEncoder):
    """Always reports can_encode=False."""
    def assess(self, spec, plan):
        return EncodingCapability(
            can_encode=False, output_dim=None, method_name=None, priority=0.0
        )

    def encode(self, spec, event, plan):
        raise AssertionError("encode() must never be called on an incapable encoder")


class _SelectiveEncoder(AtomOrbitEncoder):
    """Only capable for FLAVOURED_OPERATOR specs."""
    def assess(self, spec, plan):
        if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
            return EncodingCapability(
                can_encode=True, output_dim=8, method_name="selective", priority=1.5
            )
        return EncodingCapability(
            can_encode=False, output_dim=None, method_name=None, priority=0.0
        )

    def encode(self, spec, event, plan):
        return EncodingResult(
            values=np.ones(8, dtype=np.float64),
            metadata={"method": "selective"},
        )


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
        op = Operation("dot", rank=2, parity=1,
                       argument_symmetry=ArgumentSymmetry.SYMMETRIC)
        atom = Atom(op, ("a", "b"), sign=1)
        spec = OrbitSpec.from_atom(atom)
        assert spec.form == OrbitSpecForm.REPRESENTATIVE_ATOM
        assert spec.payload is atom

    def test_from_explicit_orbit_materialises_iterable(self):
        op = Operation("dot", rank=2, parity=1,
                       argument_symmetry=ArgumentSymmetry.SYMMETRIC)
        atom = Atom(op, ("a", "b"), sign=1)
        spec = OrbitSpec.from_explicit_orbit(iter([atom]))  # iterator, not list
        assert spec.form == OrbitSpecForm.EXPLICIT_ORBIT
        assert isinstance(spec.payload, list)
        assert len(spec.payload) == 1

    def test_from_flavoured_operator(self):
        sentinel = object()
        spec = OrbitSpec.from_flavoured_operator(sentinel)
        assert spec.form == OrbitSpecForm.FLAVOURED_OPERATOR
        assert spec.payload is sentinel

    def test_repr_is_informative(self):
        spec = _dummy_spec()
        r = repr(spec)
        assert "EXPLICIT_ORBIT" in r


# ---------------------------------------------------------------------------
# Registry: registration
# ---------------------------------------------------------------------------

class TestRegistration:

    def test_empty_registry_has_zero_length(self):
        reg = AtomOrbitEncoderRegistry()
        assert len(reg) == 0

    def test_register_increases_length(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableEncoder())
        assert len(reg) == 1
        reg.register(_IncapableEncoder())
        assert len(reg) == 2

    def test_iteration_yields_encoders_in_registration_order(self):
        reg = AtomOrbitEncoderRegistry()
        a = _CapableEncoder()
        b = _IncapableEncoder()
        reg.register(a)
        reg.register(b)
        listed = list(reg)
        assert listed[0] is a
        assert listed[1] is b


# ---------------------------------------------------------------------------
# Registry: query_all
# ---------------------------------------------------------------------------

class TestQueryAll:

    def test_returns_capable_encoders_only(self):
        reg = AtomOrbitEncoderRegistry()
        enc_yes = _CapableEncoder()
        enc_no  = _IncapableEncoder()
        reg.register(enc_yes)
        reg.register(enc_no)
        results = reg.query_all(_dummy_spec(), plan=None)
        assert len(results) == 1
        enc, cap = results[0]
        assert enc is enc_yes
        assert cap.can_encode is True

    def test_empty_result_when_none_capable(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_IncapableEncoder())
        assert reg.query_all(_dummy_spec(), plan=None) == []

    def test_results_in_registration_order(self):
        reg = AtomOrbitEncoderRegistry()
        a = _CapableEncoder(priority=1.0)
        b = _CapableEncoder(priority=2.0)
        reg.register(a)
        reg.register(b)
        results = reg.query_all(_dummy_spec(), plan=None)
        assert results[0][0] is a
        assert results[1][0] is b

    def test_selective_encoder_filtered_by_form(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_SelectiveEncoder())
        # EXPLICIT_ORBIT form → not capable
        assert reg.query_all(_dummy_spec(OrbitSpecForm.EXPLICIT_ORBIT), None) == []
        # FLAVOURED_OPERATOR form → capable
        result = reg.query_all(_dummy_spec(OrbitSpecForm.FLAVOURED_OPERATOR), None)
        assert len(result) == 1
        assert result[0][1].method_name == "selective"

    def test_capability_fields_are_passed_through(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableEncoder(priority=3.7, output_dim=12, method="foo"))
        _, cap = reg.query_all(_dummy_spec(), None)[0]
        assert cap.output_dim == 12
        assert cap.priority == 3.7
        assert cap.method_name == "foo"


# ---------------------------------------------------------------------------
# Registry: best_for
# ---------------------------------------------------------------------------

class TestBestFor:

    def test_returns_none_when_no_capable(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_IncapableEncoder())
        assert reg.best_for(_dummy_spec(), plan=None) is None

    def test_returns_highest_priority(self):
        reg = AtomOrbitEncoderRegistry()
        low  = _CapableEncoder(priority=0.5, output_dim=2)
        high = _CapableEncoder(priority=2.0, output_dim=6)
        reg.register(low)
        reg.register(high)
        enc, cap = reg.best_for(_dummy_spec(), plan=None)
        assert enc is high
        assert cap.priority == 2.0

    def test_registration_order_breaks_ties(self):
        reg = AtomOrbitEncoderRegistry()
        first  = _CapableEncoder(priority=1.0, method="first")
        second = _CapableEncoder(priority=1.0, method="second")
        reg.register(first)
        reg.register(second)
        enc, _ = reg.best_for(_dummy_spec(), plan=None)
        # Python's max is stable; first registered wins on equal priority
        assert enc is first

    def test_ignores_incapable_encoder_even_if_registered_last(self):
        reg = AtomOrbitEncoderRegistry()
        capable = _CapableEncoder(priority=0.1)
        reg.register(capable)
        reg.register(_IncapableEncoder())
        enc, _ = reg.best_for(_dummy_spec(), plan=None)
        assert enc is capable


# ---------------------------------------------------------------------------
# Registry: encode (convenience dispatch)
# ---------------------------------------------------------------------------

class TestEncode:

    def test_raises_when_no_encoder_capable(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_IncapableEncoder())
        with pytest.raises(RuntimeError, match="No registered encoder"):
            reg.encode(_dummy_spec(), event={}, plan=None)

    def test_dispatches_to_best_encoder(self):
        reg = AtomOrbitEncoderRegistry()
        low  = _CapableEncoder(priority=0.5, output_dim=2)
        high = _CapableEncoder(priority=2.0, output_dim=6)
        reg.register(low)
        reg.register(high)
        result = reg.encode(_dummy_spec(), event={}, plan=None)
        assert result.values.shape == (6,)  # came from high-priority encoder

    def test_result_values_are_float64(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableEncoder(output_dim=4))
        result = reg.encode(_dummy_spec(), event={}, plan=None)
        assert result.values.dtype == np.float64

    def test_result_metadata_is_dict(self):
        reg = AtomOrbitEncoderRegistry()
        reg.register(_CapableEncoder())
        result = reg.encode(_dummy_spec(), event={}, plan=None)
        assert isinstance(result.metadata, dict)


# ---------------------------------------------------------------------------
# EncodingCapability / EncodingResult defaults
# ---------------------------------------------------------------------------

class TestDataclassDefaults:

    def test_capability_metadata_defaults_to_empty_dict(self):
        cap = EncodingCapability(can_encode=True, output_dim=4,
                                 method_name="x", priority=1.0)
        assert cap.metadata == {}

    def test_result_metadata_defaults_to_empty_dict(self):
        res = EncodingResult(values=np.array([1.0, 2.0]))
        assert res.metadata == {}

    def test_two_capabilities_with_same_args_are_equal(self):
        a = EncodingCapability(can_encode=True, output_dim=4,
                               method_name="x", priority=1.0)
        b = EncodingCapability(can_encode=True, output_dim=4,
                               method_name="x", priority=1.0)
        assert a == b
