"""
symcoder.encoders.pair_base
============================
Core abstractions for the Pair Orbit Encoder registry framework — the Phase 2
analogue of the Atom Orbit Encoder framework (_base.py / _registry.py).

A *pair orbit* is the G-orbit of a pair of atoms (u, v) under the symmetry
group.  In the current encoding scheme these are the ASSOC segments produced
in Phase 2 of encode_and_describe().

Public names
------------
PairOrbitSpec           — wraps a PairFlavour (the pair-encoding unit)
PairOrbitEncoder        — ABC for a ready-to-use encoder bound to one pair spec
PairOrbitEncoderFactory — ABC for a factory that creates PairOrbitEncoders via assess()
PairOrbitEncoderRegistry — holds registered factories; query_all / best_for / encode

Design
------
Mirrors the Atom Orbit Encoder hierarchy exactly:

  PairOrbitEncoderFactory
      Long-lived.  assess(spec, plan) returns a list of PairOrbitEncoder
      instances — one per encoding the factory will offer for that spec.
      An empty list means "I cannot handle this spec".

  PairOrbitEncoder
      Short-lived, ready-to-use.  Created by a factory's assess(); carries all
      structural info (output_dim, priority, method_name) and exposes
      encode(event) to do the actual work.  No spec or plan needed at encode
      time — everything was captured at construction.

  PairOrbitEncoderRegistry
      Holds factories.  query_all(spec, plan) calls assess on every factory
      and flattens the results.  best_for() picks the highest-priority encoder.

Drop decisions (NULL_SELF, NULL_COMP) remain in the encode_and_describe()
orchestrator; the registry is consulted only for non-dropped ASSOC pairs.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from symatom.rep import PairFlavour
    from symatom.context import Plan


# ---------------------------------------------------------------------------
# EncodingResult  (re-exported from _base for convenience; same type)
# ---------------------------------------------------------------------------
# We re-use the same EncodingResult dataclass from the Atom side rather than
# duplicating it, since the semantics are identical.
from ._base import EncodingResult  # noqa: E402  (after TYPE_CHECKING block)


# ---------------------------------------------------------------------------
# PairOrbitSpec
# ---------------------------------------------------------------------------

@dataclass
class PairOrbitSpec:
    """
    Wraps a PairFlavour to communicate one pair orbit to a factory.

    Always use the class-method constructor from_pair_flavour() so that intent
    is clear at the call site.

    Fields
    ------
    pf : PairFlavour
        The pair flavour describing this orbit's structure.
    """
    pf: Any  # symatom.rep.PairFlavour — Any to avoid hard circular import

    @classmethod
    def from_pair_flavour(cls, pf: PairFlavour) -> PairOrbitSpec:
        """Specify the pair orbit by its PairFlavour."""
        return cls(pf=pf)

    def __repr__(self) -> str:
        return f"PairOrbitSpec({self.pf!r})"


# ---------------------------------------------------------------------------
# PairOrbitEncoder  (abstract base — the ready-to-use encoder)
# ---------------------------------------------------------------------------

class PairOrbitEncoder(ABC):
    """
    A ready-to-use encoder bound to one specific pair orbit.

    Instances are created by PairOrbitEncoderFactory.assess() and carry all
    structural information computed during that call.  encode(event) requires
    only the event data — no spec or plan is supplied again.

    Properties (inspected by managers before selecting an encoder)
    --------------------------------------------------------------
    output_dim  : int    — number of reals the embedding will produce
    priority    : float  — self-assigned preference score; higher wins
    method_name : str    — short human-readable label for the strategy used
    """

    @property
    @abstractmethod
    def output_dim(self) -> int:
        """Number of real-valued outputs the embedding produces."""
        ...

    @property
    @abstractmethod
    def priority(self) -> float:
        """Self-assigned preference score.  Higher means preferred."""
        ...

    @property
    def method_name(self) -> str | None:
        """Short label for the encoding strategy.  Override to provide one."""
        return None

    @abstractmethod
    def encode(self, event: dict) -> EncodingResult:
        """
        Perform the embedding for the given event.

        event maps label strings to numpy arrays (one per vector in the context).
        The returned values array must have dtype float64 and length output_dim.
        """
        ...


# ---------------------------------------------------------------------------
# PairOrbitEncoderFactory  (abstract base — the factory)
# ---------------------------------------------------------------------------

class PairOrbitEncoderFactory(ABC):
    """
    Abstract base for objects that create PairOrbitEncoders for pair orbits they support.

    assess(spec, plan) returns a list of PairOrbitEncoder instances — one per
    distinct encoding the factory will offer.  An empty list signals "I cannot
    handle this spec".
    """

    @abstractmethod
    def assess(self, spec: PairOrbitSpec, plan: Plan) -> list[PairOrbitEncoder]:
        """
        Return a list of ready-to-use encoders for this spec, or [] if unsupported.

        Contract
        --------
        - Must never raise for a well-formed spec; return [] instead.
        - Should be cheap relative to a full encode: no evaluation of event data.
        - Every returned encoder's output_dim must be correctly set.
        """
        ...


# ---------------------------------------------------------------------------
# PairOrbitEncoderRegistry
# ---------------------------------------------------------------------------

class PairOrbitEncoderRegistry:
    """
    A collection of PairOrbitEncoderFactory objects that can be queried and
    dispatched to.

    Typical usage
    -------------
        pair_registry = PairOrbitEncoderRegistry()
        pair_registry.register(EmbedCompressedEncoderFactory())

        spec = PairOrbitSpec.from_pair_flavour(pf)

        # Survey all factories and get back ready-to-use encoders:
        encoders = pair_registry.query_all(spec, plan)
        chosen   = min(encoders, key=lambda e: e.output_dim)  # manager selects
        result   = chosen.encode(event)

        # Convenience: use the highest-priority encoder directly:
        result = pair_registry.encode(spec, event, plan)

    Registration order
    ------------------
    Factories are queried in registration order.  query_all() returns results in
    that same order within each factory's output.  When two capable encoders have
    equal priority, best_for() returns the one that appears first in query_all()
    (stable sort).
    """

    def __init__(self) -> None:
        self._factories: list[PairOrbitEncoderFactory] = []

    # ------------------------------------------------------------------
    # Registration
    # ------------------------------------------------------------------

    def register(self, factory: PairOrbitEncoderFactory) -> None:
        """Add a factory to the registry (appended to the end)."""
        self._factories.append(factory)

    def __len__(self) -> int:
        """Return the number of registered factories."""
        return len(self._factories)

    def __iter__(self):
        """Iterate over registered factories in registration order."""
        return iter(self._factories)

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def query_all(
        self,
        spec: PairOrbitSpec,
        plan: Plan,
    ) -> list[PairOrbitEncoder]:
        """
        Ask every registered factory for encoders that can handle spec.

        Returns a flat list of PairOrbitEncoder instances in factory-registration
        order, with each factory's own encoders in the order it returned them.
        Factories that returned [] are silently skipped.
        """
        results: list[PairOrbitEncoder] = []
        for factory in self._factories:
            results.extend(factory.assess(spec, plan))
        return results

    def best_for(
        self,
        spec: PairOrbitSpec,
        plan: Plan,
    ) -> PairOrbitEncoder | None:
        """
        Return the highest-priority encoder among all candidates, or None.

        Ties are broken by position in query_all() (earlier wins).
        """
        candidates = self.query_all(spec, plan)
        if not candidates:
            return None
        return max(candidates, key=lambda e: e.priority)

    # ------------------------------------------------------------------
    # Convenience dispatch
    # ------------------------------------------------------------------

    def encode(
        self,
        spec: PairOrbitSpec,
        event: dict,
        plan: Plan,
    ) -> EncodingResult:
        """
        Encode spec using the highest-priority available encoder.

        Raises RuntimeError if no registered factory can handle spec.
        """
        enc = self.best_for(spec, plan)
        if enc is None:
            raise RuntimeError(
                f"No registered pair factory produced an encoder for {spec!r}. "
                f"({len(self._factories)} factory/factories registered.)"
            )
        return enc.encode(event)
