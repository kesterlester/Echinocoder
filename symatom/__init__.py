"""
symatom — symbolic atom and orbit machinery for permutation + continuous
symmetry encodings.

Public API:
  atoms   : ArgumentSymmetry, Operation, VectorGroup, Atom, are_negatives
  context : Context, Plan
  canon   : Canonicaliser (Protocol), SimpleCanonicaliser
  orbits  : orbit, stabiliser_size, orbit_and_stabiliser_size
"""
from .atoms   import ArgumentSymmetry, Operation, VectorGroup, Atom, are_negatives
from .context import Context, Plan
from .canon   import Canonicaliser, SimpleCanonicaliser
from .orbits  import orbit, stabiliser_size, orbit_and_stabiliser_size
from .rep     import Flavour, FlavouredOperator, repL, repS

__all__ = [
    "ArgumentSymmetry",
    "Operation",
    "VectorGroup",
    "Atom",
    "are_negatives",
    "Context",
    "Plan",
    "Canonicaliser",
    "SimpleCanonicaliser",
    "orbit",
    "stabiliser_size",
    "orbit_and_stabiliser_size",
    "Flavour",
    "FlavouredOperator",
    "repL",
    "repS",
]
