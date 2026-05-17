"""
symcoder.encoders.pair_base
============================
Core abstractions for the row-pair encoder level of the Phase 2 hierarchy.

A *pair orbit* is the G-orbit of a pair of atoms (u, v) under the symmetry
group.  In the encoding scheme these are the ASSOC segments produced in Phase 2.

Public names
------------
PairOrbitSpec           — wraps a PairFlavour (the pair-encoding unit)
PairOrbitEncoder        — ABC for a ready-to-use encoder bound to one pair spec
PairOrbitEncoderFactory — ABC for a factory that creates PairOrbitEncoders

Design (constructor-injection, no external registry)
-----------------------------------------------------
There is no PairOrbitEncoderRegistry class.  Instead, higher-level factories
(OverlapBlockEncoderFactory) receive a list of PairOrbitEncoderFactory instances
at construction time.  This makes the composition explicit and avoids a separate
mutable registry object.

  PairOrbitEncoderFactory
      assess(spec, plan) returns a list of PairOrbitEncoder instances — one per
      encoding the factory offers for that spec.  [] means "I cannot handle this".

  PairOrbitEncoder
      encode(event) performs the embedding.
      describe() returns SegmentInfo metadata (start=0 placeholder; caller adjusts).
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

    @abstractmethod
    def describe(self) -> list:
        """
        Return a list of SegmentInfo objects describing this encoder's output.

        start values are 0-relative placeholders; the calling overlap-block
        encoder adjusts them to correct positions within the block.
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

