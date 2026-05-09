"""
symcoder — numerical encoding layer built on top of symatom.

Takes symatom's symbolic atom/orbit machinery and produces concrete
invariant embeddings of physics events via polynomial zipping.

Public API:
  eval   : EvaluableOperation, evaluate
  encode : encode
"""
from .eval   import EvaluableOperation, evaluate
from .encode import encode

__all__ = [
    "EvaluableOperation",
    "evaluate",
    "encode",
]
