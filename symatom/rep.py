from __future__ import annotations
import math
from dataclasses import dataclass
from itertools import combinations as _combinations, product as _iproduct
from .atoms import Atom, ArgumentSymmetry, Operation
from .context import Context


# ---------------------------------------------------------------------------
# Internal sort key used by PairFlavour canonical ordering
# ---------------------------------------------------------------------------

def _op_fl_key(op: Operation, fl: Flavour) -> tuple:
    """Comparable sort key for an (Operation, Flavour) pair."""
    return (op.name, op.rank, op.parity, op.argument_symmetry.value, fl.counts)


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

    def describe(self, group_names) -> str:
        """Return a human-readable string using group names, e.g. '2 electrons + 1 muon'."""
        parts = [f"{k} {name}" for k, name in zip(self.counts, group_names) if k > 0]
        return " + ".join(parts) if parts else "empty"

    def __repr__(self):
        return f"Flavour({', '.join(str(c) for c in self.counts)})"

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

    def __repr__(self):
        group_names = [g.name for g in self.context.groups]
        fl_str = self.flavour.describe(group_names)
        mode = "repL" if self.signed else "repS"
        return f"FO({self.operation.name}, {fl_str}, {mode})"

    def canonical_representative(self, canonicaliser) -> Atom:
        """Return the canonical orbit representative for this FlavouredOperator."""
        return canonicaliser.canonicalise((next(self.atoms()),), self.context)[0]

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


# ---------------------------------------------------------------------------
# PairFlavour
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PairFlavour:
    """
    The orbit-type of a pair of atoms (u, v) under the full symmetry group G.

    Two atom-pairs lie in the same G-orbit if and only if they have the same
    PairFlavour (provided all labels are distinct, as the well-formedness rules
    guarantee).  The distinct PairFlavours therefore directly index the distinct
    G-orbits of atom-pairs — equivalently, the non-redundant row-pairs for the
    encoding layer above symatom.

    Fields:
        op_u, flavour_u  — operation and flavour of atom u
        op_v, flavour_v  — operation and flavour of atom v
        overlap          — per-group count of labels shared between u and v

    Canonical ordering: the (op_u, flavour_u) side is always ≤ (op_v, flavour_v)
    under an internal sort key.  Construction automatically swaps the sides if
    needed, so PairFlavour(A, fA, B, fB, ov) == PairFlavour(B, fB, A, fA, ov).
    This mirrors the fact that SimpleCanonicaliser sorts atoms within a tuple,
    identifying the G-orbits of (u, v) and (v, u).

    Signs are not tracked.  Sign-related encoding redundancies — e.g. that
    zip(rowA, rowB) and zip(rowA, −rowB) carry the same information — are
    encoding-layer concerns handled above symatom.
    """
    op_u:      Operation
    flavour_u: Flavour
    op_v:      Operation
    flavour_v: Flavour
    overlap:   tuple        # one non-negative int per group

    def __post_init__(self):
        n_groups = len(self.overlap)
        if len(self.flavour_u.counts) != n_groups:
            raise ValueError(
                f"flavour_u has {len(self.flavour_u.counts)} entries but "
                f"overlap has {n_groups}"
            )
        if len(self.flavour_v.counts) != n_groups:
            raise ValueError(
                f"flavour_v has {len(self.flavour_v.counts)} entries but "
                f"overlap has {n_groups}"
            )
        for i, (ku, kv, s) in enumerate(
            zip(self.flavour_u.counts, self.flavour_v.counts, self.overlap)
        ):
            if not isinstance(s, int) or s < 0:
                raise ValueError(
                    f"overlap[{i}] must be a non-negative integer, got {s!r}"
                )
            if s > min(ku, kv):
                raise ValueError(
                    f"overlap[{i}]={s} exceeds "
                    f"min(flavour_u[{i}]={ku}, flavour_v[{i}]={kv})"
                )
        # Impose canonical ordering: (op_u, flavour_u) <= (op_v, flavour_v).
        # Overlap is symmetric so it needs no adjustment when the sides are swapped.
        key_u = _op_fl_key(self.op_u, self.flavour_u)
        key_v = _op_fl_key(self.op_v, self.flavour_v)
        if key_u > key_v:
            op_u, fl_u = self.op_u, self.flavour_u
            object.__setattr__(self, 'op_u',      self.op_v)
            object.__setattr__(self, 'flavour_u', self.flavour_v)
            object.__setattr__(self, 'op_v',      op_u)
            object.__setattr__(self, 'flavour_v', fl_u)

    def __repr__(self):
        def _side(op, fl):
            counts_str = ", ".join(str(c) for c in fl.counts)
            return f"{op.name}[{counts_str}]"
        ov_str = ", ".join(str(s) for s in self.overlap)
        return (
            f"PairFlavour({_side(self.op_u, self.flavour_u)}"
            f" × {_side(self.op_v, self.flavour_v)}, overlap=({ov_str}))"
        )

    def count(self, group_sizes: tuple) -> int:
        """
        Number of ordered atom-pairs (u, v) with this PairFlavour.

        For group i with size n_i, flavour counts k_u and k_v, overlap s:

            ways_i = C(n_i, s) * C(n_i-s, k_u-s) * C(n_i-k_u, k_v-s)

        Total = product over all groups.  Counts ordered pairs; when
        op_u == op_v and flavour_u == flavour_v, each unordered pair {u, v}
        with u != v appears twice.  Sign variants are not counted here.
        """
        if len(group_sizes) != len(self.overlap):
            raise ValueError(
                f"group_sizes has {len(group_sizes)} entries but "
                f"overlap has {len(self.overlap)}"
            )
        total = 1
        for n, ku, kv, s in zip(
            group_sizes, self.flavour_u.counts, self.flavour_v.counts, self.overlap
        ):
            if n - ku - kv + s < 0:
                return 0    # group too small; should not arise for valid PairFlavours
            total *= (
                math.comb(n, s)
                * math.comb(n - s,  ku - s)
                * math.comb(n - ku, kv - s)
            )
        return total


# ---------------------------------------------------------------------------
# Functions that work with PairFlavour
# ---------------------------------------------------------------------------

def pair_flavour_of(atom_u: Atom, atom_v: Atom, context: Context) -> PairFlavour:
    """
    Compute the PairFlavour of a concrete atom-pair given a context.

    For each group, counts how many of atom_u's labels and atom_v's labels
    belong to that group, and how many are shared between the two atoms.
    """
    labels_v_set = set(atom_v.labels)
    flu, flv, ov = [], [], []
    for g in context.groups:
        g_set = set(g.labels)
        ku = sum(1 for lbl in atom_u.labels if lbl in g_set)
        kv = sum(1 for lbl in atom_v.labels if lbl in g_set)
        s  = sum(1 for lbl in atom_u.labels if lbl in g_set and lbl in labels_v_set)
        flu.append(ku)
        flv.append(kv)
        ov.append(s)
    return PairFlavour(
        op_u=atom_u.operation, flavour_u=Flavour(tuple(flu)),
        op_v=atom_v.operation, flavour_v=Flavour(tuple(flv)),
        overlap=tuple(ov),
    )


def canonical_pair_flavours(fo_list, context: Context) -> list:
    """
    Return all distinct PairFlavours for the given FlavouredOperators and
    context, generated directly without materialising atom-pairs.

    For each ordered pair (fo_u, fo_v) of FlavouredOperators the valid overlap
    range for group i is:
        s ∈ [max(0, k_u_i + k_v_i − n_i),  min(k_u_i, k_v_i)]

    PairFlavour's canonical ordering ensures (fo_u, fo_v) and (fo_v, fo_u)
    contribute identical PairFlavours; a seen-set handles deduplication.

    Returns a deterministically sorted list.
    """
    group_sizes = tuple(g.size for g in context.groups)
    seen = set()
    for fo_u in fo_list:
        for fo_v in fo_list:
            per_group = []
            valid = True
            for n, ku, kv in zip(
                group_sizes, fo_u.flavour.counts, fo_v.flavour.counts
            ):
                s_lo = max(0, ku + kv - n)
                s_hi = min(ku, kv)
                if s_lo > s_hi:
                    valid = False
                    break
                per_group.append(range(s_lo, s_hi + 1))
            if not valid:
                continue
            for overlap in _iproduct(*per_group):
                pf = PairFlavour(
                    op_u=fo_u.operation, flavour_u=fo_u.flavour,
                    op_v=fo_v.operation, flavour_v=fo_v.flavour,
                    overlap=overlap,
                )
                seen.add(pf)
    return sorted(seen, key=lambda pf: (
        _op_fl_key(pf.op_u, pf.flavour_u),
        _op_fl_key(pf.op_v, pf.flavour_v),
        pf.overlap,
    ))


def brute_force_canonical_pair_flavours(fo_list, context: Context) -> set:
    """
    Reference implementation for cross-checking canonical_pair_flavours().

    Materialises every atom from every FlavouredOperator, forms all ordered
    atom-pairs, and collects their PairFlavours.  O(N²) in the total atom
    count N.  Not suitable for large contexts; use only in tests.
    """
    all_atoms = [atom for fo in fo_list for atom in fo.atoms()]
    return {
        pair_flavour_of(u, v, context)
        for u in all_atoms
        for v in all_atoms
    }
