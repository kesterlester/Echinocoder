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

Architectural relationship to Phase 1 / Phase 2
-----------------------------------------------
Phase 3 mirrors the Phase-1 and Phase-2 factory pattern in spirit (sub-
factories expose ``assess`` to report whether they can offer a bound
encoder; the stacking wrapper collects them and raises if everyone
refuses) but diverges in two deliberate ways:

1. **No per-orbit "spec" argument.**  Phase 1's
   ``AtomOrbitEncoderFactory.assess(spec, plan)`` is called once per
   FlavouredOperator with the per-orbit spec.  Phase 2's
   ``OverlapBlockEncoderFactory.assess(spec, plan)`` is called once per
   overlap block.  Phase 3 is a *whole-table* encoder — there is no
   per-orbit or per-block iteration to do — so its sub-factory
   ``assess(plan)`` takes the plan alone and is called once.  Adding a
   vacuous ``spec`` would be cargo-cult.

2. **Concatenate, not min-select.**  Phase 1 picks
   ``min(offered, key=output_dim)`` per orbit because its alternative
   encoders are *interchangeable* — they all losslessly encode the same
   orbit, so storing two would be redundant.  Phase 3's sub-encoders are
   *complementary*: simplicial vs Vandermonde-Cinf vs Vandermonde-C0
   trade off smoothness, output magnitude, and (in principle) faithfulness
   against each other in non-overlapping ways, and may legitimately be
   useful side-by-side to downstream consumers.  Concatenating is the
   honest default; users who want only one pass a one-element factory
   list.

In short: same ``assess``-then-collect skeleton, different selection
policy and a one-shot rather than per-orbit dispatch.

Refusal semantics
-----------------
A sub-factory's ``assess(plan)`` may return ``[]`` when the user's
construction parameters cannot be honoured on the given plan — for
instance, a ``Phase3VandermondeEncoderFactory(scale=λ)`` applied to a
plan whose operations don't all declare a ``mass_dimension`` returns
``[]``.  The stacking wrapper silently skips refusing sub-factories,
but raises a clear ``RuntimeError`` if *every* sub-factory refuses
(producing a zero-length composite is almost never what the caller
wanted).

Examples
--------
Default (simplicial only)::

    Phase3EncoderFactory([Phase3SimplicialEncoderFactory()])

Both Vandermonde variants together::

    Phase3EncoderFactory([
        Phase3VandermondeEncoderFactory(mode="Cinf"),
        Phase3VandermondeEncoderFactory(mode="C0"),
    ])

All three side by side, with the Cinf variant dimensionally rescaled to
a ~1 GeV reference scale::

    Phase3EncoderFactory([
        Phase3SimplicialEncoderFactory(),
        Phase3VandermondeEncoderFactory(mode="Cinf", scale=1.0),
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
        Iterable of Phase-3 sub-factories.  Each must have an
        ``assess(plan)`` method returning a list of bound encoders (often a
        singleton, sometimes empty — see below) and a ``build(plan)``
        method that returns a single bound encoder or raises.

    Plan compatibility and refusal
    ------------------------------
    Each sub-factory's ``assess(plan)`` may legitimately return an empty
    list, indicating that the sub-factory cannot honour the user's
    construction parameters on this plan (e.g.\ a Vandermonde factory
    with ``scale=λ`` applied to a plan whose operations don't declare
    ``mass_dimension``).  Such sub-factories are silently skipped; the
    composite encoder is built from whichever sub-factories *do* offer.

    If *every* sub-factory refuses, ``build(plan)`` raises a clear error
    rather than silently returning an empty composite — Phase 3 was
    explicitly requested by the user and producing a zero-length
    encoding is almost certainly not what they wanted.

    See module docstring for usage examples.
    """

    def __init__(self, factories: Iterable) -> None:
        self._factories = list(factories)

    def build(self, plan) -> Phase3CompositeEncoder:
        encoders = []
        refusals = []
        for f in self._factories:
            offered = f.assess(plan) if hasattr(f, "assess") else [f.build(plan)]
            if not offered:
                refusals.append(type(f).__name__)
            encoders.extend(offered)
        if not encoders and self._factories:
            raise RuntimeError(
                f"Phase3EncoderFactory.build: every sub-factory refused to "
                f"offer a bound encoder for this plan.  Refusing factories: "
                f"{refusals}.  Common cause: a Vandermonde factory was "
                f"constructed with scale set, but the plan contains "
                f"operation(s) without a declared mass_dimension; either "
                f"drop the scale (scale=None) or declare a mass_dimension "
                f"on every Operation in the plan."
            )
        return Phase3CompositeEncoder(encoders)
