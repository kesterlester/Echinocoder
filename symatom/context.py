from __future__ import annotations
from dataclasses import dataclass, field
from .atoms import VectorType, Operation, Atom
from .group import TheGroup
from .orbit_enum import OrbitEnumerator, BruteForceOrbitEnumerator, DirectOrbitEnumerator


@dataclass(frozen=True)
class Context:
    """
    The set of VectorTypes in scope for a computation.  Immutable once created.
    The full symmetry group is the direct product of the S_n for each group.
    """
    types:     tuple     # tuple of VectorType
    the_group: TheGroup  = field(init=False, repr=False)

    def __post_init__(self):
        if not isinstance(self.types, tuple):
            raise TypeError(f"types must be a tuple, got {type(self.types)}")
        names = [g.name for g in self.types]
        if len(set(names)) != len(names):
            raise ValueError(f"VectorType names must be distinct, got {names!r}")
        all_labels = []
        for g in self.types:
            all_labels.extend(g.labels)
        if len(set(all_labels)) != len(all_labels):
            raise ValueError("Labels must be distinct across all VectorTypes")
        object.__setattr__(self, 'the_group', TheGroup(self.types))

    def type_of(self, label) -> VectorType:
        """Return the VectorType that contains label, or raise KeyError."""
        for g in self.types:
            if label in g.labels:
                return g
        raise KeyError(f"Label {label!r} not found in any VectorType")

    @property
    def all_labels(self) -> tuple:
        result = []
        for g in self.types:
            result.extend(g.labels)
        return tuple(result)

    def check_atom(self, atom: Atom) -> None:
        """Raise ValueError if any label in atom is not in this context."""
        for label in atom.labels:
            self.type_of(label)   # raises KeyError if absent

    def __repr__(self):
        parts = ", ".join(f"{g.name}[{g.size}]" for g in self.types)
        return f"Context({parts})"


@dataclass
class Plan:
    """
    Bundles the local configuration for a computation: which vector types are
    in scope, which operations are registered, and which orbit enumeration
    strategy to use.

    Plans are passed explicitly; there is no global plan.
    """
    context:          Context
    operations:       tuple            = field(default_factory=tuple)
    orbit_enumerator: OrbitEnumerator  = field(default_factory=DirectOrbitEnumerator)

    def get_operation(self, name: str) -> Operation:
        for op in self.operations:
            if op.name == name:
                return op
        raise KeyError(f"Operation {name!r} is not registered in this plan")

    def __repr__(self):
        ops = ", ".join(op.name for op in self.operations)
        return f"Plan(context={self.context!r}, ops=({ops}))"
