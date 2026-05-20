"""
symcoder.decoded_types
=======================
Output types for Step-A decoding (invariant-space reconstruction).

Every decode() call returns one of these annotated-multiset types.  The value
list is semantically a multiset: callers must not rely on element order.  The
atom / atom-pair reference records which algebraic objects produced the values,
where the bijection between values and atoms is unknown (G-ambiguity).

Values are stored as plain Python lists rather than sets so that floating-point
near-duplicates are never silently merged or dropped.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class AnnotatedMultisetOfReals:
    """
    A multiset of real values produced by a Phase 1 orbit decoder.

    values    : list of floats — the full decoded orbit, fully decompressed.
                Treat as a multiset; do not rely on order.
    atoms     : list of Atom objects — the algebraic atoms whose evaluations
                produced these values, in an unknown correspondence.
    """
    values: list   # list[float]
    atoms:  list   # list[Atom]

    def __repr__(self) -> str:
        return (f"AnnotatedMultisetOfReals("
                f"{len(self.values)} values, "
                f"{len(self.atoms)} atoms)")


@dataclass
class AnnotatedMultisetOfRealPairs:
    """
    A multiset of (u_value, v_value) real pairs produced by a Phase 2
    row-pair decoder.

    pairs      : list of (float, float) tuples — the full decoded orbit of
                 pairs, fully decompressed.  Treat as a multiset; do not rely
                 on order.
    atom_pairs : list of (Atom, Atom) tuples — the algebraic atom-pairs whose
                 evaluations produced the pairs, in an unknown correspondence.
    """
    pairs:      list   # list[tuple[float, float]]
    atom_pairs: list   # list[tuple[Atom, Atom]]

    def __repr__(self) -> str:
        return (f"AnnotatedMultisetOfRealPairs("
                f"{len(self.pairs)} pairs, "
                f"{len(self.atom_pairs)} atom-pairs)")
