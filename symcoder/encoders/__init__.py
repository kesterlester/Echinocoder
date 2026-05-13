"""
symcoder.encoders — pluggable registry of atom orbit encoders.

Public API
----------
OrbitSpecForm           — enum: REPRESENTATIVE_ATOM | FLAVOURED_OPERATOR | EXPLICIT_ORBIT
OrbitSpec               — wraps one orbit specification (use class-method constructors)
EncodingCapability      — what an encoder says it can do (returned by assess())
EncodingResult          — embedding output: float64 numpy array + metadata dict
AtomOrbitEncoder        — ABC; subclass this to write a new encoder
AtomOrbitEncoderRegistry — holds registered encoders; query_all / best_for / encode
SortEncoder             — stub: sorted-evaluation encoder (all orbit forms, low priority)
PolyEncoder             — stub: polynomial-compressed encoder (ANTISYM orbits, hi priority)

See encoders/README.md for a full concept overview and usage walkthrough.
"""
from ._base import (
    OrbitSpecForm,
    OrbitSpec,
    EncodingCapability,
    EncodingResult,
    AtomOrbitEncoder,
)
from ._registry import AtomOrbitEncoderRegistry
from .sort_encoder import SortEncoder
from .poly_encoder import PolyEncoder

__all__ = [
    "OrbitSpecForm",
    "OrbitSpec",
    "EncodingCapability",
    "EncodingResult",
    "AtomOrbitEncoder",
    "AtomOrbitEncoderRegistry",
    "SortEncoder",
    "PolyEncoder",
]
