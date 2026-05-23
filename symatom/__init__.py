"""
symatom — symbolic atom and orbit machinery for permutation + continuous
symmetry encodings.

Public API:
  atoms      : ArgumentSymmetry, Operation, VectorType, Atom, are_negatives
  context    : Context, Plan
  orbit_enum : OrbitEnumerator, BruteForceOrbitEnumerator, DirectOrbitEnumerator
  rep        : Flavour, FlavouredOperator, repS
  group      : TheGroup, GroupElement
"""
from .atoms      import ArgumentSymmetry, Operation, VectorType, Atom, are_negatives
from .context    import Context, Plan
from .orbit_enum import OrbitEnumerator, BruteForceOrbitEnumerator, DirectOrbitEnumerator
from .rep        import (Flavour, FlavouredOperator, repS,
                         PairFlavour, pair_flavour_of,
                         canonical_pair_flavours,
                         brute_force_canonical_pair_flavours)
from .group      import TheGroup, GroupElement

__all__ = [
    "ArgumentSymmetry",
    "Operation",
    "VectorType",
    "Atom",
    "are_negatives",
    "Context",
    "Plan",
    "OrbitEnumerator",
    "BruteForceOrbitEnumerator",
    "DirectOrbitEnumerator",
    "Flavour",
    "FlavouredOperator",
    "repS",
    "PairFlavour",
    "pair_flavour_of",
    "canonical_pair_flavours",
    "brute_force_canonical_pair_flavours",
    "TheGroup",
    "GroupElement",
]
