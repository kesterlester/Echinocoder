from __future__ import annotations
from dataclasses import dataclass
from enum import Enum


class ArgumentSymmetry(Enum):
    """
    Declares how an operation's value transforms under permutations of its own
    arguments.  Designed as a subclassable Enum so that a MixedSymmetry variant
    can be added later without breaking existing code.
    """
    SYMMETRIC     = "symmetric"
    ANTISYMMETRIC = "antisymmetric"
    UNSTRUCTURED  = "unstructured"


@dataclass(frozen=True)
class Operation:
    """
    A named, reusable description of a type of multilinear contraction.
    The framework is agnostic about what the operation *is*; it only uses
    these declared properties.
    """
    name:               str
    rank:               int             # number of vector arguments (>= 1)
    parity:             int             # +1 or -1 under spatial inversion
    argument_symmetry:  ArgumentSymmetry

    def __post_init__(self):
        if self.rank < 1:
            raise ValueError(f"rank must be >= 1, got {self.rank!r}")
        if self.parity not in (+1, -1):
            raise ValueError(f"parity must be +1 or -1, got {self.parity!r}")


@dataclass(frozen=True)
class VectorGroup:
    """
    A named, ordered collection of vector labels all of the same physical type.
    The symmetric group S_n acts on it by permuting its labels.
    """
    name:   str
    labels: tuple   # tuple of hashable labels, typically strings

    def __post_init__(self):
        if not isinstance(self.labels, tuple):
            raise TypeError(f"labels must be a tuple, got {type(self.labels)}")
        if len(self.labels) == 0:
            raise ValueError("VectorGroup must have at least one label")
        if len(set(self.labels)) != len(self.labels):
            raise ValueError(f"VectorGroup labels must be distinct, got {self.labels!r}")

    @property
    def size(self) -> int:
        return len(self.labels)


@dataclass(frozen=True)
class Atom:
    """
    An operation applied to an ordered tuple of vector labels, with an overall
    sign.  The sign is the *only* place a sign lives; canonicalisation absorbs
    all sign changes into this field and never returns an external sign.

    Well-formedness rules (enforced at construction):
      1. len(labels) == operation.rank
      2. sign in {+1, -1}
      3. sign == +1  whenever  operation.argument_symmetry != ANTISYMMETRIC
         (a dot product can never legitimately carry sign = -1)
    Rule 2 (label membership in a context) is checked separately by Context.
    """
    operation:  Operation
    labels:     tuple   # ordered tuple of vector labels
    sign:       int     # +1 or -1

    def __post_init__(self):
        if not isinstance(self.labels, tuple):
            raise TypeError(f"labels must be a tuple, got {type(self.labels)}")
        if len(self.labels) != self.operation.rank:
            raise ValueError(
                f"Operation '{self.operation.name}' has rank {self.operation.rank} "
                f"but {len(self.labels)} label(s) were supplied: {self.labels!r}"
            )
        if self.sign not in (+1, -1):
            raise ValueError(f"sign must be +1 or -1, got {self.sign!r}")
        if (self.sign == -1 and
                self.operation.argument_symmetry != ArgumentSymmetry.ANTISYMMETRIC):
            raise ValueError(
                f"sign=-1 is not valid for operation '{self.operation.name}' whose "
                f"argument_symmetry is {self.operation.argument_symmetry}; "
                f"only ANTISYMMETRIC operations may carry sign=-1"
            )


def are_negatives(a: Atom, b: Atom) -> bool:
    """
    Return True iff a and b are negatives of each other: same operation, same
    labels, opposite signs.  Always False for non-ANTISYMMETRIC operations
    (consistent with the well-formedness rule that such atoms carry sign=+1).
    """
    return (
        a.operation == b.operation
        and a.labels == b.labels
        and a.sign == -b.sign
        and a.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
    )
