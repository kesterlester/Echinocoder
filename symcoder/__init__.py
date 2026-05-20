"""
symcoder — numerical encoding layer built on top of symatom.

Takes symatom's symbolic atom/orbit machinery and produces concrete
invariant embeddings of physics events via polynomial zipping.

Public API:
  eval         : EvaluableOperation, evaluate
  encode       : encode, encode_and_describe, describe_encoding
  describe     : SegmentInfo
  encoders     : OrbitEncoderFactory (Phase 1 top level)
                 SortEncoderFactory, PolyEncoderFactory
                 standard_row_pair_factories (Phase 2 row-pair convenience)
                 OverlapBlockEncoderFactory  (Phase 2 block level)
                 Phase2EncoderFactory        (Phase 2 top level)
"""
from .eval     import EvaluableOperation, evaluate
from .encode   import encode, encode_and_describe, describe_encoding
from .describe import SegmentInfo, Phase1Tree, OverlapBlockNode, Phase2Tree, EncodingTree
from .encoders import (
    OrbitEncoderFactory,
    SortEncoderFactory, HalfSortEncoderFactory, PolyEncoderFactory,
    standard_row_pair_factories,
    OverlapBlockEncoderFactory,
    Phase2Encoder, Phase2EncoderFactory,
)

__all__ = [
    "EvaluableOperation",
    "evaluate",
    "encode",
    "encode_and_describe",
    "describe_encoding",
    "SegmentInfo",
    "Phase1Tree",
    "OverlapBlockNode",
    "Phase2Tree",
    "EncodingTree",
    "OrbitEncoderFactory",
    "SortEncoderFactory",
    "HalfSortEncoderFactory",
    "PolyEncoderFactory",
    "standard_row_pair_factories",
    "OverlapBlockEncoderFactory",
    "Phase2Encoder",
    "Phase2EncoderFactory",
]
