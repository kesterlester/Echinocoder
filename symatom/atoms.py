from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable
from enum import Enum


def _sort_sign(labels: tuple) -> tuple[tuple, int]:
    """
    Return (sorted_labels, sign) where sign is the parity of the permutation
    that sorts the labels.  Uses bubble sort to count adjacent transpositions.
    Assumes labels are mutually comparable (e.g. all strings).
    """
    lst = list(labels)
    sign = 1
    for i in range(len(lst)):
        for j in range(len(lst) - 1 - i):
            if lst[j] > lst[j + 1]:
                lst[j], lst[j + 1] = lst[j + 1], lst[j]
                sign *= -1
    return tuple(lst), sign


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

    eval_fn parameter
    -----------------
    eval_fn(vectors) -> float
        vectors : list of arrays, one per label, in atom.labels order.
        Return value : the numerical value of the operation on those vectors,
        BEFORE the atom's sign is applied.  The sign is applied by evaluate().
        eval_fn is included in equality and hashing using Python's default
        function identity (f == g iff f is g).  This means two Operation
        objects with different callable objects are considered distinct even
        if they happen to compute the same value — which is the desired
        behaviour for library singletons (e.g. euclidean2.dot ≠ euclidean3.dot).

    Optional tex parameter
    ----------------------
    tex : str | None
        A LaTeX macro *body* using #1, #2, … as positional placeholders for
        the atom's labels, e.g. r"\\varepsilon(\\vc{#1},\\vc{#2},\\vc{#3})".
        When supplied, TeX emitters can define a \\newcommand with this body
        (rank arguments) and call it with the atom's labels.  LaTeX supports
        at most 9 macro parameters, so passing tex with rank > 9 raises
        ValueError at construction time.
        tex is excluded from equality and hashing — it is display metadata,
        not part of the symbolic identity of the operation.
    """
    name:               str
    rank:               int             # number of vector arguments (>= 1)
    odd_parity:         bool            # True → parity = -1 (pseudoscalar); False → parity = +1 (scalar)
    argument_symmetry:  ArgumentSymmetry
    eval_fn:            Callable = field(kw_only=True)
    tex:                str | None = field(default=None, kw_only=True,
                                           compare=False, hash=False)

    def __post_init__(self):
        if self.rank < 1:
            raise ValueError(f"rank must be >= 1, got {self.rank!r}")
        if self.tex is not None and self.rank > 9:
            raise ValueError(
                f"tex macro requires at most 9 LaTeX parameters; "
                f"operation '{self.name}' has rank {self.rank}"
            )

    @property
    def parity(self) -> int:
        """Return the spatial-parity sign as +1 or -1 (derived from odd_parity)."""
        return -1 if self.odd_parity else +1

    def __repr__(self):
        par = "+1" if not self.odd_parity else "-1"
        sym = self.argument_symmetry.name
        tex_part = f", tex={self.tex!r}" if self.tex is not None else ""
        return f"Operation('{self.name}', rank={self.rank}, parity={par}, {sym}{tex_part})"


@dataclass(frozen=True)
class VectorType:
    """
    A named, ordered collection of vector labels all of the same physical type.
    The symmetric group S_n acts on it by permuting its labels.
    """
    name:   str
    labels: tuple   # tuple of hashable labels, typically strings

    def __post_init__(self):
        if not isinstance(self.labels, tuple):
            raise TypeError(f"labels must be a tuple, got {type(self.labels)}")
        if len(set(self.labels)) != len(self.labels):
            raise ValueError(f"VectorType labels must be distinct, got {self.labels!r}")

    def __repr__(self):
        lbls = ", ".join(str(l) for l in self.labels)
        return f"VectorType('{self.name}': {lbls})"

    @property
    def size(self) -> int:
        return len(self.labels)


def _canonicalise_fields(
    operation: Operation, labels: tuple, sign: int, *, pin_sign: bool
) -> tuple[tuple, int]:
    """
    Validate and canonicalise (labels, sign) for an Atom.  Shared by both
    construction paths so that changing the canonicalisation rule only ever
    requires editing this one function.

    pin_sign=False  (normal path): for ANTISYMMETRIC operations the sign is
                    multiplied by the parity of the label-sorting permutation;
                    for SYMMETRIC and UNSTRUCTURED operations the sign is
                    stored exactly as given.
    pin_sign=True   (pinned path): the sign is stored exactly as given;
                    labels are still sorted but parity is never computed.

    All three ArgumentSymmetry values accept sign ∈ {+1, -1}.  The sign
    represents the algebraic sign of the atom's value (e.g. -dot(a,b) is a
    valid atom with SYMMETRIC operation and sign=-1).

    Returns (canonical_labels, final_sign).
    Raises TypeError / ValueError on malformed input.
    """
    if not isinstance(labels, tuple):
        raise TypeError(f"labels must be a tuple, got {type(labels)}")
    if len(labels) != operation.rank:
        raise ValueError(
            f"Operation '{operation.name}' has rank {operation.rank} "
            f"but {len(labels)} label(s) were supplied: {labels!r}"
        )
    if sign not in (+1, -1):
        raise ValueError(f"sign must be +1 or -1, got {sign!r}")
    if len(set(labels)) != len(labels):
        raise ValueError(f"Atom labels must be distinct, got {labels!r}")

    sym = operation.argument_symmetry
    if sym == ArgumentSymmetry.UNSTRUCTURED:
        return labels, sign

    if pin_sign:
        # Sort labels but ignore parity — caller's sign is stored as-is.
        sorted_labels, _ = _sort_sign(labels)
        return sorted_labels, sign
    else:
        sorted_labels, perm_sign = _sort_sign(labels)
        final_sign = sign * perm_sign if sym == ArgumentSymmetry.ANTISYMMETRIC else sign
        return sorted_labels, final_sign


@dataclass(frozen=True)
class Atom:
    """
    An operation applied to an ordered tuple of vector labels, with an overall
    sign.  The sign is the *only* place a sign lives; canonicalisation absorbs
    all sign changes into this field and never returns an external sign.

    Well-formedness rules (enforced at construction):
      1. len(labels) == operation.rank
      2. sign in {+1, -1}  — all ArgumentSymmetry values permit both signs.
      3. all labels supplied at initialisation must be distinct. So if forming dot(a,b) is valid then dot(a,a) is ill-formed;
         use a rank-1 lenSq taking just a if you want something equivalent to evaluating dot(a,a).

    Internal argument canonicalization (also enforced at construction):
      - SYMMETRIC:     labels are sorted into ascending order; sign is unchanged.
                       e.g. Atom(dot, ("b","a"), sign=-1) stores labels=("a","b"), sign=-1.
      - ANTISYMMETRIC: labels are sorted into ascending order; sign is multiplied
                       by the parity of the sorting permutation.
      - UNSTRUCTURED:  labels are stored in the order given; no reordering.

    As a result, Atom(dot, ("b","a"), +1) == Atom(dot, ("a","b"), +1), and
    are_negatives(Atom(eps, ("b","a"), +1), Atom(eps, ("a","b"), +1)) is True
    (the unsorted eps atom self-canonicalizes to the sorted form with sign -1).

    Rule 2 (label membership in a context) is checked separately by Context.

    Alternative constructor
    -----------------------
    Atom.with_pinned_sign(operation, labels, sign) sorts the labels into
    canonical order but preserves the given sign exactly. The parity of the
    sorting permutation is not folded in.  Use this when the sign is already
    known to be correct for the canonical label arrangement.
    """
    operation:  Operation
    labels:     tuple   # stored in canonical argument order (see above)
    sign:       int     # +1 or -1

    def __post_init__(self):
        labels, sign = _canonicalise_fields(
            self.operation, self.labels, self.sign, pin_sign=False
        )
        object.__setattr__(self, 'labels', labels)
        object.__setattr__(self, 'sign',   sign)

    @classmethod
    def with_pinned_sign(cls, operation: Operation, labels, sign: int = +1) -> "Atom":
        """
        Construct an Atom with canonical (sorted) label order, but preserving sign
        exactly as given; the parity of the sorting permutation is NOT folded in.

        Bypasses __init__ / __post_init__ entirely so no perm-parity is computed
        only to be discarded.  All validation still runs via _canonicalise_fields.
        """
        canon_labels, _ = _canonicalise_fields(
            operation, tuple(labels), sign, pin_sign=True
        )
        inst = object.__new__(cls)
        object.__setattr__(inst, 'operation', operation)
        object.__setattr__(inst, 'labels',    canon_labels)
        object.__setattr__(inst, 'sign',      sign)   # pin_sign=True guarantees no change
        return inst

    def __repr__(self):
        sign_str = "+" if self.sign == 1 else "-"
        args = ", ".join(str(l) for l in self.labels)
        return f"{sign_str}{self.operation.name}({args})"


def are_negatives(a: Atom, b: Atom) -> bool:
    """
    Return True iff a and b are negatives of each other: same operation, same
    canonical labels, opposite signs.  Works for all ArgumentSymmetry values.
    """
    return (
        a.operation == b.operation
        and a.labels == b.labels
        and a.sign == -b.sign
    )
