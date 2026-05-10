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
from .encode   import encode
from .describe import describe_encoding, SegmentInfo

__all__ = [
    "EvaluableOperation",
    "evaluate",
    "encode",
    "describe_encoding",
    "SegmentInfo",
]
