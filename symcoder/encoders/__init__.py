"""
symcoder.encoders — pluggable encoder hierarchies for Phase 1 and Phase 2.

Phase 1 — Atom Orbit Encoders
------------------------------
OrbitSpecForm               — enum: REPRESENTATIVE_ATOM | FLAVOURED_OPERATOR | EXPLICIT_ORBIT
OrbitSpec                   — wraps one orbit specification (use class-method constructors)
EncodingResult              — embedding output: float64 numpy array + metadata dict
AtomOrbitEncoder            — ABC for a ready-to-use encoder bound to one orbit spec
AtomOrbitEncoderFactory     — ABC for a factory that creates AtomOrbitEncoders via assess()
AtomOrbitEncoderRegistry    — holds registered factories; query_all / best_for / encode
SortEncoder / SortEncoderFactory       — sort-based Phase 1 encoder (fallback)
PolyEncoder / PolyEncoderFactory       — polynomial Phase 1 encoder (ANTISYM, higher priority)

Phase 2 — Row-pair Encoders (lowest level)
--------------------------------------------
PairOrbitSpec               — wraps a PairFlavour (the pair-encoding unit)
PairOrbitEncoder            — ABC for a ready-to-use encoder bound to one pair spec
PairOrbitEncoderFactory     — ABC for a factory that creates PairOrbitEncoders

NullPairEncoder / SelfPairEncoderFactory        — null encoding for self-pairs (output_dim=0)
Type11PairEncoder / Type11PairEncoderFactory    — TYPE_11 pairs
Type12PairEncoder / Type12PairEncoderFactory    — TYPE_12 pairs
Type21PairEncoder / Type21PairEncoderFactory    — TYPE_21 pairs
Type22PairEncoder / Type22PairEncoderFactory    — TYPE_22 pairs
NegPairEncoder / NegPairEncoderFactory          — TYPE_NEG pairs
standard_row_pair_factories()                   — convenience: all six factories in order

Phase 2 — Overlap Block Encoders (middle level)
-------------------------------------------------
OverlapBlockSpec              — wraps the PairFlavour list for one overlap block
OverlapBlockEncoder           — ready-to-use encoder for one block
OverlapBlockEncoderFactory    — factory: takes row-pair factories; manages selection
                                and complementarity drop (use_complementarity_drop kwarg)

Phase 2 — Top Level (all blocks)
----------------------------------
Phase2Encoder                 — encodes all overlap blocks; delegates to block encoders
Phase2EncoderFactory          — factory: takes block factories; build(plan) → Phase2Encoder

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
from .sort_encoder import SortEncoder, SortEncoderFactory, HalfSortEncoder, HalfSortEncoderFactory
from .poly_encoder import PolyEncoder, PolyEncoderFactory
from .pair_base import (
    PairOrbitSpec,
    PairOrbitEncoder,
    PairOrbitEncoderFactory,
)
from .row_pair_encoders import (
    NullPairEncoder,        SelfPairEncoderFactory,
    Type11PairEncoder,      Type11PairEncoderFactory,
    Type12PairEncoder,      Type12PairEncoderFactory,
    Type21PairEncoder,      Type21PairEncoderFactory,
    Type22PairEncoder,      Type22PairEncoderFactory,
    NegPairEncoder,         NegPairEncoderFactory,
    standard_row_pair_factories,
)
from .overlap_block import (
    OverlapBlockSpec,
    OverlapBlockEncoder,
    OverlapBlockEncoderFactory,
)
from .phase2_encoder import Phase2Encoder, Phase2EncoderFactory
from .orbit_encoder import OrbitEncoder, OrbitEncoderFactory

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
    "HalfSortEncoder",
    "HalfSortEncoderFactory",
    "PolyEncoder",
    "PolyEncoderFactory",
    # Phase 2 row-pair level
    "PairOrbitSpec",
    "PairOrbitEncoder",
    "PairOrbitEncoderFactory",
    "NullPairEncoder",        "SelfPairEncoderFactory",
    "Type11PairEncoder",      "Type11PairEncoderFactory",
    "Type12PairEncoder",      "Type12PairEncoderFactory",
    "Type21PairEncoder",      "Type21PairEncoderFactory",
    "Type22PairEncoder",      "Type22PairEncoderFactory",
    "NegPairEncoder",         "NegPairEncoderFactory",
    "standard_row_pair_factories",
    # Phase 2 block level
    "OverlapBlockSpec",
    "OverlapBlockEncoder",
    "OverlapBlockEncoderFactory",
    # Phase 2 top level
    "Phase2Encoder",
    "Phase2EncoderFactory",
    # Phase 1 top level
    "OrbitEncoder",
    "OrbitEncoderFactory",
]
