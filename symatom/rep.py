from __future__ import annotations
import math
from dataclasses import dataclass
from itertools import combinations as _combinations, product as _iproduct
from .atoms import Atom, ArgumentSymmetry, Operation
from .context import Context


# ---------------------------------------------------------------------------
# Flavour
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Flavour:
    """
    Label-composition tuple: how many arguments an operation draws from each
    vector group.  For a context with m groups, Flavour((k_1,...,k_m)) means
    k_i arguments come from group i.  The rank of the associated operation
    must equal sum(counts).

    Flavour is context-agnostic — it only carries counts; the context is held
    by FlavouredOperator.
    """
    counts: tuple   # non-negative ints, one per group in the context

    def __post_init__(self):
        if not isinstance(self.counts, tuple):
            raise TypeError(f"Flavour.counts must be a tuple, got {type(self.counts)}")
        for c in self.counts:
            if not isinstance(c, int) or c < 0:
                raise ValueError(
                    f"Flavour counts must be non-negative integers, got {self.counts!r}"
                )

    @property
    def rank(self) -> int:
        return sum(self.counts)

    def __iter__(self):
        return iter(self.counts)

    def __len__(self):
        return len(self.counts)


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _valid_flavours(groups: tuple, rank: int):
    """Yield every Flavour consistent with group sizes and the given rank."""
    def _rec(remaining, idx, acc):
        if idx == len(groups):
            if remaining == 0:
                yield Flavour(tuple(acc))
            return
        for k in range(min(remaining, groups[idx].size) + 1):
            yield from _rec(remaining - k, idx + 1, acc + [k])
    yield from _rec(rank, 0, [])


# ---------------------------------------------------------------------------
# FlavouredOperator
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FlavouredOperator:
    """
    A specific operation applied to a specific Flavour within a given context.
    Acts as a lazy recipe for the group of atoms sharing the same operation and
    label-composition.

    signed=True  → repL semantics: for ANTISYMMETRIC operations both sign=+1
                   and sign=-1 atoms are included for each label combination.
    signed=False → repS semantics: for ANTISYMMETRIC operations only the sign=+1
                   atom per combination is included.
    Symmetric and UNSTRUCTURED operations are identical under repL and repS.

    Mixed-symmetry operations are not yet supported and will raise
    NotImplementedError if encountered (future extension point).
    """
    operation:  Operation
    flavour:    Flavour
    context:    Context
    signed:     bool

    def __post_init__(self):
        if len(self.flavour) != len(self.context.groups):
            raise ValueError(
                f"Flavour has {len(self.flavour)} counts but context has "
                f"{len(self.context.groups)} groups"
            )
        if self.flavour.rank != self.operation.rank:
            raise ValueError(
                f"Flavour rank {self.flavour.rank} != operation rank {self.operation.rank}"
            )
        for group, k in zip(self.context.groups, self.flavour):
            if k > group.size:
                raise ValueError(
                    f"Flavour count {k} exceeds group '{group.name}' size {group.size}"
                )

    def count(self) -> int:
        """
        Number of atoms this FlavouredOperator will yield.  Pure combinatorics —
        no atoms are generated.
        """
        base = math.prod(
            math.comb(g.size, k)
            for g, k in zip(self.context.groups, self.flavour)
        )
        if self.signed and self.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC:
            return base * 2
        return base

    def atoms(self):
        """
        Lazily yield all atoms for this (operation, flavour) combination.

        Labels are generated via combinations (one per species group), so each
        combination is visited exactly once with labels in sorted order within
        each group.  For ANTISYMMETRIC operations with signed=True (repL), both
        sign=+1 and sign=-1 variants are yielded for every combination.
        """
        per_group = [
            list(_combinations(g.labels, k))
            for g, k in zip(self.context.groups, self.flavour)
        ]
        antisym_and_signed = (
            self.signed and
            self.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        )
        for combo_parts in _iproduct(*per_group):
            labels = tuple(lbl for part in combo_parts for lbl in part)
            yield Atom(self.operation, labels, sign=+1)
            if antisym_and_signed:
                yield Atom(self.operation, labels, sign=-1)

    def contains(self, atom: Atom) -> bool:
        """
        Return True if atom belongs to the vocabulary of this FlavouredOperator.

        Checks operation identity, that every label belongs to the correct group
        in the correct count (matching the flavour), and that the sign is
        consistent with the signed/unsigned semantics.  Does not require the
        atom's labels to be in any particular order.
        """
        if atom.operation != self.operation:
            return False
        label_group = {
            lbl: i
            for i, g in enumerate(self.context.groups)
            for lbl in g.labels
        }
        counts = [0] * len(self.context.groups)
        for lbl in atom.labels:
            idx = label_group.get(lbl)
            if idx is None:
                return False
            counts[idx] += 1
        if tuple(counts) != self.flavour.counts:
            return False
        if atom.sign == -1 and not self.signed:
            return False
        return True


# ---------------------------------------------------------------------------
# repL and repS
# ---------------------------------------------------------------------------

def repL(context: Context, operations) -> list:
    """
    Return a list of FlavouredOperator (signed=True) covering every valid
    (operation, Flavour) combination for the given context and operations.

    "Valid" means each species contributes no more labels than its pool allows.
    """
    return [
        FlavouredOperator(operation=op, flavour=fl, context=context, signed=True)
        for op in operations
        for fl in _valid_flavours(context.groups, op.rank)
    ]


def repS(context: Context, operations) -> list:
    """
    Return a list of FlavouredOperator (signed=False) covering every valid
    (operation, Flavour) combination for the given context and operations.
    """
    return [
        FlavouredOperator(operation=op, flavour=fl, context=context, signed=False)
        for op in operations
        for fl in _valid_flavours(context.groups, op.rank)
    ]
