"""
symcoder.encoders — pluggable registry of atom orbit encoders (Phase 1) and
                    pair orbit encoders (Phase 2).

Phase 1 — Atom Orbit Encoders
------------------------------
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

Phase 2 — Pair Orbit Encoders
-------------------------------
PairOrbitSpec               — wraps a PairFlavour (the pair-encoding unit)
PairOrbitEncoder            — ABC for a ready-to-use encoder bound to one pair spec
PairOrbitEncoderFactory     — ABC for a factory that creates PairOrbitEncoders via assess()
PairOrbitEncoderRegistry    — holds registered factories; query_all / best_for / encode
EmbedCompressedEncoder      — concrete encoder wrapping _embed_compressed pipeline
EmbedCompressedEncoderFactory — factory: non-self pairs, priority 0.5

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
from .pair_base import (
    PairOrbitSpec,
    PairOrbitEncoder,
    PairOrbitEncoderFactory,
    PairOrbitEncoderRegistry,
)
from .embed_compressed_encoder import EmbedCompressedEncoder, EmbedCompressedEncoderFactory

__all__ = [
    # Phase 1
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
    # Phase 2
    "PairOrbitSpec",
    "PairOrbitEncoder",
    "PairOrbitEncoderFactory",
    "PairOrbitEncoderRegistry",
    "EmbedCompressedEncoder",
    "EmbedCompressedEncoderFactory",
]
