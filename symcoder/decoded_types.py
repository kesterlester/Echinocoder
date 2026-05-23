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


@dataclass
class AnnotatedMultisetOfRepSEvalVectors:
    """
    The G-orbit of eval(repS, E): the output of the alignment decoder.

    This is the richest Step-A decoded object.  It contains all G-invariant
    information present in the encoding — the full correlational structure of
    atom evaluations across the discrete group orbit — and no more.

    vectors : list of tuples of floats — the multiset of column vectors, one
              per element of the G-orbit of repS at event E.  Each tuple has
              length |repS|; position r holds eval(g · repS[r], E) for the
              group element g corresponding to that column.  Treat as a
              multiset; do not rely on order.  The multiset has
              |G| / |stab(repS, E)| elements; for generic events (trivial
              stabiliser) this equals |G|.
    atoms   : list of Atom objects — repS, the canonical ordered list of all
              atoms across all FlavouredOperators.  The same sequence underlies
              every column vector; it provides the algebraic row labels of the
              alignment table.  The bijection between vector positions and group
              elements is unknown (G-ambiguity, as always).
    """
    vectors: list   # list[tuple[float, ...]]
    atoms:   list   # list[Atom]  (repS)

    def __repr__(self) -> str:
        return (f"AnnotatedMultisetOfRepSEvalVectors("
                f"{len(self.vectors)} vectors, "
                f"{len(self.atoms)} atoms)")
