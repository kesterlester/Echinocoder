"""
symcoder.encoders._base
=======================
Core abstractions for the Atom Orbit Encoder registry framework.

Public names
------------
OrbitSpecForm       — enum discriminator for the three orbit-specification forms
OrbitSpec           — thin wrapper carrying one orbit specification
EncodingCapability  — what an encoder says it can do (returned by assess())
EncodingResult      — what an encoder returns after embedding (returned by encode())
AtomOrbitEncoder    — abstract base class; every concrete encoder subclasses this
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    # Avoid importing at runtime for callers that do not use these types directly.
    from symatom.atoms import Atom
    from symatom.rep import FlavouredOperator
    from symatom.context import Plan


# ---------------------------------------------------------------------------
# OrbitSpec
# ---------------------------------------------------------------------------

class OrbitSpecForm(Enum):
    """
    Discriminator for the three ways an atom orbit can be communicated to an encoder.

    REPRESENTATIVE_ATOM
        A single Atom that is a member of the orbit.  The encoder uses the Atom's
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

    Encoders inspect spec.form to decide whether they can handle the spec, then
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
# EncodingCapability  (returned by assess())
# ---------------------------------------------------------------------------

@dataclass
class EncodingCapability:
    """
    What an AtomOrbitEncoder says it can (or cannot) do for a given OrbitSpec.

    Returned by AtomOrbitEncoder.assess(); used by AtomOrbitEncoderRegistry
    to filter and rank encoders.

    Fields
    ------
    can_encode : bool
        True iff this encoder is able to embed the specified orbit.
        If False, all other fields may be None / meaningless.

    output_dim : int | None
        Number of real-valued outputs the embedding would produce.
        Set to None when can_encode is False.

    method_name : str | None
        Short human-readable label for the embedding method, e.g. "sort" or
        "poly_compressed".  Useful for logging and diagnostics.

    priority : float
        Self-assigned quality / preference score.  Higher means the encoder
        considers itself a better fit for this orbit.  The high-level manager
        will typically call query_all() and apply its own selection logic;
        priority is used only by the registry's convenience best_for() helper.

    metadata : dict[str, Any]
        Any additional information the encoder wishes to expose, e.g. estimated
        compute cost, whether the embedding is injective, precision notes, etc.
    """
    can_encode:  bool
    output_dim:  int | None = None
    method_name: str | None = ""
    priority:    float = 0
    metadata:    dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# EncodingResult  (returned by encode())
# ---------------------------------------------------------------------------

@dataclass
class EncodingResult:
    """
    What an AtomOrbitEncoder returns after performing an embedding.

    Fields
    ------
    values : np.ndarray
        The ordered real-number embedding, dtype float64, shape (output_dim,).
        numpy arrays are used (rather than plain lists) because the encoding
        pipeline is designed to run over many events at speed.

    metadata : dict[str, Any]
        Diagnostic information about this particular embedding, e.g. which
        encoder was used, the orbit size encountered, normalisation applied.
        May be empty.
    """
    values:   np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# AtomOrbitEncoder  (abstract base class)
# ---------------------------------------------------------------------------

class AtomOrbitEncoder(ABC):
    """
    Abstract base for objects that embed a single-atom orbit into a vector of reals.

    Design overview
    ---------------
    Each concrete subclass specialises in certain kinds of orbit (e.g. orbits of
    symmetric operations, orbits with a specific sign-correlation structure, or
    orbits in contexts with particular group structure).  The two-step protocol
    — assess() then encode() — lets a manager survey all registered encoders
    before committing to one:

        registry = AtomOrbitEncoderRegistry()
        registry.register(SortEncoder())
        registry.register(PolyEncoder())

        spec = OrbitSpec.from_flavoured_operator(fo)
        candidates = registry.query_all(spec, plan)
        # candidates: [(encoder, capability), ...] for all capable encoders
        # high-level manager picks one and calls encoder.encode(spec, event, plan)

    Lifecycle of a single encoder instance
    ----------------------------------------
    1. Instantiate — encoder may set up internal state or hold a sub-encoder.
    2. assess(spec, plan)  — cheap probe; returns EncodingCapability.
    3. encode(spec, event, plan)  — only called after a True assess(); does the work.

    Orbit specification forms
    -------------------------
    spec is an OrbitSpec wrapping one of:
      - a representative Atom            (REPRESENTATIVE_ATOM form)
      - a FlavouredOperator              (FLAVOURED_OPERATOR form)
      - an already-enumerated list[Atom] (EXPLICIT_ORBIT form)

    Concrete encoders return can_encode=False for forms they do not support.

    Internal delegation
    -------------------
    An encoder may hold another encoder instance and delegate part of its work.
    For example, a polynomial encoder might hold a SortEncoder to obtain Phase-1
    sorted evaluations, then apply polynomial compression on top:

        class PolyEncoder(AtomOrbitEncoder):
            def __init__(self):
                self._sort_enc = SortEncoder()
            def encode(self, spec, event, plan):
                sorted_vals = self._sort_enc.encode(spec, event, plan).values
                # ... compress further ...

    Extending to higher-level objects
    -----------------------------------
    The same assess / encode / OrbitSpec / Capability / Result / Registry pattern
    can be replicated for higher-level objects (atom pairs, molecules, composite
    objects that are not themselves Atoms but contain Atoms).  Create a parallel
    ABC (e.g. MoleculeOrbitEncoder) and a corresponding registry class.  The
    AtomOrbitEncoderRegistry has no knowledge of the concrete encoder types, so
    it can be replicated without modification.
    """

    @abstractmethod
    def assess(self, spec: OrbitSpec, plan: Plan) -> EncodingCapability:
        """
        Return capability information without performing any embedding.

        Contract
        --------
        - Must never raise for a well-formed spec; return can_encode=False instead.
        - Should be cheap: no evaluation of event data, minimal combinatorics.
        - If can_encode is True, output_dim must be set to the correct value.
        """
        ...

    @abstractmethod
    def encode(self, spec: OrbitSpec, event: dict, plan: Plan) -> EncodingResult:
        """
        Perform the embedding and return a vector of reals.

        Contract
        --------
        - Only called when assess() returned can_encode=True for the same spec.
        - event maps label strings to numpy arrays (one per vector in the context).
        - The returned values array must have dtype float64 and shape (output_dim,)
          where output_dim matches what assess() reported.
        """
        ...
