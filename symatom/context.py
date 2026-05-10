from __future__ import annotations
from dataclasses import dataclass, field
from .atoms import VectorGroup, Operation, Atom
from .orbit_enum import OrbitEnumerator, BruteForceOrbitEnumerator, DirectOrbitEnumerator
from .canon import DirectCanonicaliser


@dataclass(frozen=True)
class Context:
    """
    The set of VectorGroups in scope for a computation.  Immutable once created.
    The full symmetry group is the direct product of the S_n for each group.
    """
    groups: tuple   # tuple of VectorGroup

    def __post_init__(self):
        if not isinstance(self.groups, tuple):
            raise TypeError(f"groups must be a tuple, got {type(self.groups)}")
        names = [g.name for g in self.groups]
        if len(set(names)) != len(names):
            raise ValueError(f"VectorGroup names must be distinct, got {names!r}")
        all_labels = []
        for g in self.groups:
            all_labels.extend(g.labels)
        if len(set(all_labels)) != len(all_labels):
            raise ValueError("Labels must be distinct across all VectorGroups")

    def group_of(self, label) -> VectorGroup:
        """Return the VectorGroup that contains label, or raise KeyError."""
        for g in self.groups:
            if label in g.labels:
                return g
        raise KeyError(f"Label {label!r} not found in any VectorGroup")

    @property
    def all_labels(self) -> tuple:
        result = []
        for g in self.groups:
            result.extend(g.labels)
        return tuple(result)

    def check_atom(self, atom: Atom) -> None:
        """Raise ValueError if any label in atom is not in this context."""
        for label in atom.labels:
            self.group_of(label)   # raises KeyError if absent

    def __repr__(self):
        parts = ", ".join(f"{g.name}[{g.size}]" for g in self.groups)
        return f"Context({parts})"


@dataclass
class Plan:
    """
    Bundles the local configuration for a computation: which vector groups are
    in scope, which canonicalisation implementation to use, which operations
    are registered, and which orbit enumeration strategy to use.

    Plans are passed explicitly; there is no global plan.
    """
    context:          Context
    canonicaliser:    object           = field(default_factory=DirectCanonicaliser)
    operations:       tuple            = field(default_factory=tuple)
    orbit_enumerator: OrbitEnumerator  = field(default_factory=DirectOrbitEnumerator)

    def canonicalise(self, atom_tuple: tuple) -> tuple:
        """Delegate to the plan's canonicaliser."""
        return self.canonicaliser.canonicalise(atom_tuple, self.context)

    def get_operation(self, name: str) -> Operation:
        for op in self.operations:
            if op.name == name:
                return op
        raise KeyError(f"Operation {name!r} is not registered in this plan")

    def __repr__(self):
        ops = ", ".join(op.name for op in self.operations)
        return f"Plan(context={self.context!r}, ops=({ops}))"
