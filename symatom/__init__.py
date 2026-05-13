"""
symatom — symbolic atom and orbit machinery for permutation + continuous
symmetry encodings.

Public API:
  atoms      : ArgumentSymmetry, Operation, VectorGroup, Atom, are_negatives
  context    : Context, Plan
  canon      : Canonicaliser (Protocol), SimpleCanonicaliser, DirectCanonicaliser
  orbit_enum : OrbitEnumerator, BruteForceOrbitEnumerator, DirectOrbitEnumerator
  orbits     : orbit, stabiliser_size, orbit_and_stabiliser_size
  rep        : Flavour, FlavouredOperator, repS
  group      : TheGroup, SignCorrelationType
"""
from .atoms      import ArgumentSymmetry, Operation, VectorGroup, Atom, are_negatives
from .context    import Context, Plan
from .canon      import Canonicaliser, SimpleCanonicaliser, DirectCanonicaliser
from .orbit_enum import OrbitEnumerator, BruteForceOrbitEnumerator, DirectOrbitEnumerator
from .orbits     import orbit, stabiliser_size, orbit_and_stabiliser_size
from .rep        import (Flavour, FlavouredOperator, repS,
                         PairFlavour, pair_flavour_of,
                         canonical_pair_flavours,
                         brute_force_canonical_pair_flavours)
from .group      import TheGroup, SignCorrelationType

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
    "DirectCanonicaliser",
    "OrbitEnumerator",
    "BruteForceOrbitEnumerator",
    "DirectOrbitEnumerator",
    "orbit",
    "stabiliser_size",
    "orbit_and_stabiliser_size",
    "Flavour",
    "FlavouredOperator",
    "repS",
    "PairFlavour",
    "pair_flavour_of",
    "canonical_pair_flavours",
    "brute_force_canonical_pair_flavours",
    "TheGroup",
    "SignCorrelationType",
]
