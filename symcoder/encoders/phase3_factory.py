"""
symcoder.encoders.phase3_factory
=================================
Stacking wrapper for Phase 3 encoders.

Phase 3 has multiple alternative one-shot whole-table encoders (simplicial,
Vandermonde-Cinf, Vandermonde-C0).  Conceptually they all encode the same
object — the alignment table T as a multiset of |G| tuples in R^{|rep|} —
but with different trade-offs:

  =================  ===============  ===========================  ==========
  Encoder            Smoothness        Output size                  Notes
  =================  ===============  ===========================  ==========
  Simplicial         C0 (spiky)        2nk + 1 - k                  Smallest
  Vandermonde Cinf   C-infinity        ((k-1)n + 1) · n             Smooth
  Vandermonde C0     C0 (sort)         ((k-1)n + 1) · n             Dim-OK
  =================  ===============  ===========================  ==========

The user picks any subset by passing a list of factories to
``Phase3EncoderFactory``.  All chosen encoders are run and their outputs are
concatenated end-to-end, just as Phase 1 stacks ``SortEncoderFactory`` and
``HalfSortEncoderFactory`` and Phase 2 stacks block factories.

Concatenation rather than min-selection?
----------------------------------------
Phase 1 picks the *minimum*-output-dim encoder per orbit (min_assess
policy), because the encoders are exact alternatives — having more than one
adds redundancy.  Phase 3 *concatenates* by default because the user
explicitly told us they want to be able to "use them all" — different
encoders here have different smoothness / scale / faithfulness trade-offs
and may be useful simultaneously to downstream consumers.  Passing a
single-element list to this factory recovers single-encoder behaviour.

Examples
--------
Default (simplicial only)::

    Phase3EncoderFactory([Phase3SimplicialEncoderFactory()])

Both Vandermonde variants together::

    Phase3EncoderFactory([
        Phase3VandermondeEncoderFactory(mode="Cinf"),
        Phase3VandermondeEncoderFactory(mode="C0"),
    ])

All three side by side::

    Phase3EncoderFactory([
        Phase3SimplicialEncoderFactory(),
        Phase3VandermondeEncoderFactory(mode="Cinf"),
        Phase3VandermondeEncoderFactory(mode="C0"),
    ])
"""
from __future__ import annotations

from typing import Iterable, List

import numpy as np

from symcoder.encoders._base import EncodingResult


def _Phase3Tree_class():
    from symcoder.describe import Phase3Tree
    return Phase3Tree


class Phase3CompositeEncoder:
    """Composite encoder that concatenates the output of multiple Phase 3
    encoders.

    Built by ``Phase3EncoderFactory.build(plan)``.  Holds the list of
    individual encoders and their cumulative output-dim layout.
    """

    def __init__(self, encoders: List) -> None:
        self._encoders = list(encoders)

    @property
    def output_dim(self) -> int:
        return sum(e.output_dim for e in self._encoders)

    def encode(self, event: dict) -> EncodingResult:
        if not self._encoders:
            return EncodingResult(values=np.array([], dtype=np.float64))
        parts = [e.encode(event).values for e in self._encoders]
        values = np.concatenate(parts) if parts else np.array([], dtype=np.float64)
        return EncodingResult(values=values)

    def describe(self, start_offset: int = 0):
        Phase3Tree = _Phase3Tree_class()
        all_segments = []
        cursor = start_offset
        # Use first encoder's (n, k) for the descriptor — all encoders should
        # report identical n, k because they all consume the same alignment
        # table.  An empty composite reports (n, k) = (None, None).
        n_val, k_val = None, None
        for enc in self._encoders:
            tree = enc.describe(start_offset=cursor)
            all_segments.extend(tree.segments)
            if n_val is None:
                n_val = getattr(tree, "n", None)
                k_val = getattr(tree, "k", None)
            cursor += enc.output_dim
        return Phase3Tree(segments=all_segments, n=n_val, k=k_val)


class Phase3EncoderFactory:
    """Stacking factory for Phase 3 encoders.

    Parameters
    ----------
    factories : iterable
        Iterable of Phase-3 sub-factories.  Each must have a ``build(plan)``
        method returning an encoder with ``encode(event)``, ``output_dim``,
        and ``describe(start_offset)``.

    See module docstring for usage examples.
    """

    def __init__(self, factories: Iterable) -> None:
        self._factories = list(factories)

    def build(self, plan) -> Phase3CompositeEncoder:
        encoders = [f.build(plan) for f in self._factories]
        return Phase3CompositeEncoder(encoders)
