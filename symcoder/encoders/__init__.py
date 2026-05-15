"""
symcoder.encoders — pluggable registry of atom orbit encoders.

Public API
----------
OrbitSpecForm               — enum: REPRESENTATIVE_ATOM | FLAVOURED_OPERATOR | EXPLICIT_ORBIT
OrbitSpec                   — wraps one orbit specification (use class-method constructors)
EncodingResult              — embedding output: float64 numpy array + metadata dict
AtomOrbitEncoder            — ABC for a ready-to-use encoder bound to one orbit spec
AtomOrbitEncoderFactory     — ABC for a factory that creates AtomOrbitEncoders via assess()
AtomOrbitEncoderRegistry    — holds registered factories; query_all / best_for / encode
SortEncoder                 — ready-to-use sort encoder (produced by SortEncoderFactory)
SortEncoderFactory          — factory: all orbit forms, low priority (0.5)
PolyEncoder                 — ready-to-use poly encoder stub (produced by PolyEncoderFactory)
PolyEncoderFactory          — factory stub: ANTISYM orbits, high priority (1.0)

See encoders/README.md for a full concept overview and usage walkthrough.
"""
from ._base import (
    OrbitSpecForm,
    OrbitSpec,
    EncodingResult,
    AtomOrbitEncoder,
    AtomOrbitEncoderFactory,
)
from ._registry import AtomOrbitEncoderRegistry
from .sort_encoder import SortEncoder, SortEncoderFactory
from .poly_encoder import PolyEncoder, PolyEncoderFactory

__all__ = [
    "OrbitSpecForm",
    "OrbitSpec",
    "EncodingResult",
    "AtomOrbitEncoder",
    "AtomOrbitEncoderFactory",
    "AtomOrbitEncoderRegistry",
    "SortEncoder",
    "SortEncoderFactory",
    "PolyEncoder",
    "PolyEncoderFactory",
]
