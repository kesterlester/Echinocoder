"""
symcoder.encoders._base
=======================
Core abstractions for the Atom Orbit Encoder registry framework.

Public names
------------
OrbitSpecForm           — enum discriminator for the three orbit-specification forms
OrbitSpec               — thin wrapper carrying one orbit specification
EncodingResult          — what an encoder returns after embedding (returned by encode())
AtomOrbitEncoder        — ABC for a ready-to-use encoder bound to one orbit spec
AtomOrbitEncoderFactory — ABC for a factory that creates AtomOrbitEncoders via assess()

Design overview
---------------
The framework uses a two-type pattern:

  AtomOrbitEncoderFactory
      A long-lived object (typically one per encoding strategy, e.g. sort or
      polynomial compression).  Its only job is assess(spec, plan), which returns
      a (possibly empty) list of AtomOrbitEncoder instances — one per distinct
      encoding it is willing to offer for that spec.  An empty list means "I
      cannot handle this spec".  Factories never hold event data.

  AtomOrbitEncoder
      A short-lived, ready-to-use encoder returned by a factory's assess().  It
      carries all structural information computed during assess (orbit, output
      dimension, etc.) and exposes encode(event) to do the actual work.  No
      spec or plan is needed at encode time — everything was captured at
      construction.  Inspectable properties (output_dim, priority, method_name)
      let a manager rank and select among competing encoders before committing.

  AtomOrbitEncoderRegistry
      Holds a list of factories.  query_all(spec, plan) calls assess on every
      factory and flattens the results.  best_for() picks the highest-priority
      encoder.  encode() is a convenience that calls best_for then .encode(event).
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from symatom.atoms import Atom
    from symatom.rep import FlavouredOperator
    from symatom.context import Plan


# ---------------------------------------------------------------------------
# OrbitSpec
# ---------------------------------------------------------------------------

class OrbitSpecForm(Enum):
    """
    Discriminator for the three ways an atom orbit can be communicated to a factory.

    REPRESENTATIVE_ATOM
        A single Atom that is a member of the orbit.  The factory uses the Atom's
        operation and label composition — together with the Plan's context — to
        deduce all structural properties of the orbit (size, sign-correlation type,
        etc.) without enumerating every member.

    FLAVOURED_OPERATOR
        A symatom.rep.FlavouredOperator, i.e. a (operation, flavour, context)
        triple that acts as a compact structural description of the orbit.
        No specific atom instance is required; the FlavouredOperator already
        encodes everything needed for counting and enumeration.

    EXPLICIT_ORBIT
        An already-enumerated list of Atom objects forming the orbit.  Useful
        when the caller has already done the enumeration for other reasons, or
        when testing with hand-crafted orbits.
    """
    REPRESENTATIVE_ATOM = "representative_atom"
    FLAVOURED_OPERATOR  = "flavoured_operator"
    EXPLICIT_ORBIT      = "explicit_orbit"


@dataclass
class OrbitSpec:
    """
    Wraps one of three representations of an atom orbit.

    Always use the class-method constructors (from_atom, from_flavoured_operator,
    from_explicit_orbit) rather than the raw constructor, so that intent is clear
    and payload types are documented at the call site.

    Factories inspect spec.form to decide whether they can handle the spec, then
    access spec.payload with the appropriate cast.
    """
    form:    OrbitSpecForm
    payload: Any  # Atom | FlavouredOperator | list[Atom]  (see OrbitSpecForm docs)

    @classmethod
    def from_atom(cls, atom: Atom) -> OrbitSpec:
        """Specify the orbit by a single representative Atom."""
        return cls(form=OrbitSpecForm.REPRESENTATIVE_ATOM, payload=atom)

    @classmethod
    def from_flavoured_operator(cls, fo: FlavouredOperator) -> OrbitSpec:
        """Specify the orbit by a FlavouredOperator (structural description)."""
        return cls(form=OrbitSpecForm.FLAVOURED_OPERATOR, payload=fo)

    @classmethod
    def from_explicit_orbit(cls, atoms) -> OrbitSpec:
        """Specify the orbit by providing all of its Atom members as an iterable."""
        return cls(form=OrbitSpecForm.EXPLICIT_ORBIT, payload=list(atoms))

    def __repr__(self) -> str:
        return f"OrbitSpec({self.form.name}, {self.payload!r})"


# ---------------------------------------------------------------------------
# EncodingResult  (returned by AtomOrbitEncoder.encode())
# ---------------------------------------------------------------------------

@dataclass
class EncodingResult:
    """
    What an AtomOrbitEncoder returns after performing an embedding.

    Fields
    ------
    values : np.ndarray
        The ordered real-number embedding, dtype float64, shape (output_dim,).

    metadata : dict[str, Any]
        Diagnostic information about this particular embedding.  May be empty.
    """
    values:   np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# AtomOrbitEncoder  (abstract base — the ready-to-use encoder)
# ---------------------------------------------------------------------------

class AtomOrbitEncoder(ABC):
    """
    A ready-to-use encoder bound to one specific orbit.

    Instances are created by AtomOrbitEncoderFactory.assess() and carry all
    structural information computed during that call.  encode(event) requires
    only the event data — no spec or plan is supplied again.

    Properties (inspected by managers before selecting an encoder)
    --------------------------------------------------------------
    output_dim  : int    — number of reals the embedding will produce
    priority    : float  — self-assigned preference score; higher wins;
                           only used to break ties when checks on other more
                           important thigns (e.g. output_dim) do not select
                           a clear winner
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
# AtomOrbitEncoderFactory  (abstract base — the factory)
# ---------------------------------------------------------------------------

class AtomOrbitEncoderFactory(ABC):
    """
    Abstract base for objects that create AtomOrbitEncoders for orbits they support.

    A factory is a long-lived, stateless-ish object (one per encoding strategy).
    assess(spec, plan) is the only required method.  It returns a list of
    AtomOrbitEncoder instances — one per distinct encoding the factory wishes to
    offer.  An empty list signals "I cannot handle this spec".

    A factory may return more than one encoder for a single spec when it offers
    multiple variants (e.g. sort-ascending and sort-descending, or different
    polynomial compression levels).

    Internal delegation
    -------------------
    A factory may hold another factory or encoder internally.  For example, a
    polynomial factory might hold a SortEncoderFactory and call its assess() to
    obtain the sorted evaluations, then wrap the result with polynomial compression.
    """

    @abstractmethod
    def assess(self, spec: OrbitSpec, plan: Plan) -> list[AtomOrbitEncoder]:
        """
        Return a list of ready-to-use encoders for this spec, or [] if unsupported.

        Contract
        --------
        - Must never raise for a well-formed spec; return [] instead.
        - Should be cheap relative to a full encode: no evaluation of event data.
        - Every returned encoder's output_dim must be correctly set.
        """
        ...
