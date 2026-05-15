"""
symcoder — numerical encoding layer built on top of symatom.

Takes symatom's symbolic atom/orbit machinery and produces concrete
invariant embeddings of physics events via polynomial zipping.

Public API:
  eval     : EvaluableOperation, evaluate
  encode   : encode
  describe : describe_encoding, SegmentInfo
"""
from .eval     import EvaluableOperation, evaluate
from .encode   import encode, encode_brute, encode_and_describe, describe_encoding
from .describe import SegmentInfo

__all__ = [
    "EvaluableOperation",
    "evaluate",
    "encode",
    "encode_brute",
    "encode_and_describe",
    "describe_encoding",
    "SegmentInfo",
]
