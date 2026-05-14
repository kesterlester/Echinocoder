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
    vector type.  For a context with m vector types, Flavour((k_1,...,k_m)) means
    k_i arguments come from vector type i.  The rank of the associated operation
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

    def describe(self, type_names) -> str:
        """Return a human-readable string using group names, e.g. '2 electrons + 1 muon'."""
        parts = [f"{k} {name}" for k, name in zip(self.counts, type_names) if k > 0]
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

def _valid_flavours_rear_first(types: tuple, rank: int):
    """Yield every Flavour consistent with group sizes and the given rank.

    Enumerates rear-first: iterates k from 0 upward at each position, so
    later types tend to be filled before earlier ones.
    Example with types of sizes (2, 2, 2) and rank 2:
      (0,0,2), (0,1,1), (0,2,0), (1,0,1), (1,1,0), (2,0,0)
    """
    def _rec(remaining, idx, acc):
        if idx == len(types):
            if remaining == 0:
                yield Flavour(tuple(acc))
            return
        for k in range(min(remaining, types[idx].size) + 1):
            yield from _rec(remaining - k, idx + 1, acc + [k])
    yield from _rec(rank, 0, [])


def _valid_flavours_front_first(types: tuple, rank: int):
    """Yield every Flavour consistent with group sizes and the given rank.

    Enumerates front-first: iterates k from max downward at each position, so
    earlier types tend to be filled before later ones.
    Example with types of sizes (2, 2, 2) and rank 2:
      (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)
    """
    def _rec(remaining, idx, acc):
        if idx == len(types):
            if remaining == 0:
                yield Flavour(tuple(acc))
            return
        for k in range(min(remaining, types[idx].size), -1, -1):
            yield from _rec(remaining - k, idx + 1, acc + [k])
    yield from _rec(rank, 0, [])


# Convenience alias — default enumeration order is front-first.
_valid_flavours = _valid_flavours_front_first


# ---------------------------------------------------------------------------
# FlavouredOperator
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FlavouredOperator:
    """
    A specific operation applied to a specific Flavour within a given context.
    Acts as a lazy recipe for the group of atoms sharing the same operation and
    label-composition.

    Mixed-symmetry operations are not yet supported and will raise
    NotImplementedError if encountered (future extension point).
    """
    operation:  Operation
    flavour:    Flavour
    context:    Context

    def __post_init__(self):
        if len(self.flavour) != len(self.context.types):
            raise ValueError(
                f"Flavour has {len(self.flavour)} counts but context has "
                f"{len(self.context.types)} vector types"
            )
        if self.flavour.rank != self.operation.rank:
            raise ValueError(
                f"Flavour rank {self.flavour.rank} != operation rank {self.operation.rank}"
            )
        for vt, k in zip(self.context.types, self.flavour):
            if k > vt.size:
                raise ValueError(
                    f"Flavour count {k} exceeds vector type '{vt.name}' size {vt.size}"
                )

    def count_of_atoms_one_per_sign(self) -> int:
        """
        Number of atoms this FlavouredOperator will yield.  Pure combinatorics —
        no atoms are generated.
        """
        base = math.prod(
            math.comb(g.size, k)
            for g, k in zip(self.context.types, self.flavour)
        )
        return base

    def atoms_one_per_sign(self):
        """
        Lazily yield all atoms for this (operation, flavour) combination.

        Labels are generated via combinations (one per species group), so each
        combination is visited exactly once with labels in sorted order within
        each group.

        TODO: This should move to a lazy use of TheGroup.
        """
        per_type = [
            list(_combinations(g.labels, k))
            for g, k in zip(self.context.types, self.flavour)
        ]
        for combo_parts in _iproduct(*per_type):
            labels = tuple(lbl for part in combo_parts for lbl in part)
            yield Atom(self.operation, labels, sign=+1)

    def __repr__(self):
        type_names = [g.name for g in self.context.types]
        fl_str = self.flavour.describe(type_names)
        return f"FO({self.operation.name}, {fl_str})"

    def canonical_representative(self) -> Atom:
        """Return the first atom from atoms_one_per_sign() as the orbit representative."""
        return next(self.atoms_one_per_sign())

    def matches_ignoring_sign(self, atom: Atom) -> bool:
        """
        Return True if atom belongs to the vocabulary of this FlavouredOperator.

        Checks operation identity, that every label belongs to the correct group
        in the correct count (matching the flavour). Does not require the
        atom's labels to be in any particular order or that the sign matches.
        """
        if atom.operation != self.operation:
            return False
        label_type = {
            lbl: i
            for i, g in enumerate(self.context.types)
            for lbl in g.labels
        }
        counts = [0] * len(self.context.types)
        for lbl in atom.labels:
            idx = label_type.get(lbl)
            if idx is None:
                return False
            counts[idx] += 1
        if tuple(counts) != self.flavour.counts:
            return False
        return True


# ---------------------------------------------------------------------------
# repS
# ---------------------------------------------------------------------------


def repS(context: Context, operations) -> list:
    """
    Return a list of FlavouredOperator) covering every valid
    (operation, Flavour) combination for the given context and operations.
    """
    return [
        FlavouredOperator(operation=op, flavour=fl, context=context)
        for op in operations
        for fl in _valid_flavours(context.types, op.rank)
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
        n_types = len(self.overlap)
        if len(self.flavour_u.counts) != n_types:
            raise ValueError(
                f"flavour_u has {len(self.flavour_u.counts)} entries but "
                f"overlap has {n_types}"
            )
        if len(self.flavour_v.counts) != n_types:
            raise ValueError(
                f"flavour_v has {len(self.flavour_v.counts)} entries but "
                f"overlap has {n_types}"
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

    def count(self, type_sizes: tuple) -> int:
        """
        Number of ordered atom-pairs (u, v) with this PairFlavour.

        For group i with size n_i, flavour counts k_u and k_v, overlap s:

            ways_i = C(n_i, s) * C(n_i-s, k_u-s) * C(n_i-k_u, k_v-s)

        Total = product over all groups.  Counts ordered pairs; when
        op_u == op_v and flavour_u == flavour_v, each unordered pair {u, v}
        with u != v appears twice.  Sign variants are not counted here.
        """
        if len(type_sizes) != len(self.overlap):
            raise ValueError(
                f"type_sizes has {len(type_sizes)} entries but "
                f"overlap has {len(self.overlap)}"
            )
        total = 1
        for n, ku, kv, s in zip(
            type_sizes, self.flavour_u.counts, self.flavour_v.counts, self.overlap
        ):
            if n - ku - kv + s < 0:
                return 0    # group too small; should not arise for valid PairFlavours
            total *= (
                math.comb(n, s)
                * math.comb(n - s,  ku - s)
                * math.comb(n - ku, kv - s)
            )
        return total

    def orbit_size(self, type_sizes: tuple) -> int:
        """
        Number of (Atom, Atom) pairs in the G-orbit of the canonical
        representative (+u_c, +v_c) of this PairFlavour.

        Computed as |G| / |Stab(+u_c, +v_c)| where G = S_{n_1} × ... × S_{n_m}
        and the stabiliser counts group elements that preserve both atoms
        *including their signs*.

        For SYMMETRIC operations the sign is always +1, so the stabiliser
        is unconstrained on sign.  For ANTISYMMETRIC operations the sign is
        the parity of the permutation on the atom's labels; the stabiliser
        must therefore have even total parity on those labels to leave the
        atom's sign unchanged.

        Assumption: labels from different groups do not interleave in the
        global sort order (i.e. all labels of one group sort before all
        labels of the next).  This is satisfied for the typical naming
        conventions used in this codebase.  If labels interleave, the
        formula may be incorrect; use len(orbit_elements(context)) instead.
        """
        if self.count(type_sizes) == 0:
            return 0

        def E(k):
            """Number of even permutations of k elements (= k!/2 for k>=2, else 1)."""
            return math.factorial(k) // 2 if k >= 2 else 1

        def O(k):
            """Number of odd permutations of k elements (= k!/2 for k>=2, else 0)."""
            return math.factorial(k) // 2 if k >= 2 else 0

        antisym_u = self.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        antisym_v = self.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

        group_order = 1
        stab = 1
        for n, ku, kv, s in zip(type_sizes,
                                  self.flavour_u.counts,
                                  self.flavour_v.counts,
                                  self.overlap):
            a = ku - s          # u-only labels in this group
            b = kv - s          # v-only labels in this group
            r = n - ku - kv + s  # labels used by neither u nor v in this group
            group_order *= math.factorial(n)
            if antisym_u and antisym_v:
                # Both signs must remain +1: need parity(τ_s)·parity(τ_u) = +1
                # AND parity(τ_s)·parity(τ_v) = +1 simultaneously.
                stab_g = (E(s) * E(a) * E(b) + O(s) * O(a) * O(b)) * math.factorial(r)
            elif not antisym_u and antisym_v:
                # u sign always +1 (symmetric); need parity(τ_s)·parity(τ_v) = +1.
                stab_g = (E(s) * E(b) + O(s) * O(b)) * math.factorial(a) * math.factorial(r)
            elif antisym_u and not antisym_v:
                # v sign always +1 (symmetric); need parity(τ_s)·parity(τ_u) = +1.
                stab_g = (E(s) * E(a) + O(s) * O(a)) * math.factorial(b) * math.factorial(r)
            else:
                # Both symmetric: no parity constraints on the stabiliser.
                stab_g = (math.factorial(s) * math.factorial(a)
                          * math.factorial(b) * math.factorial(r))
            stab *= stab_g
        return group_order // stab

    def orbit_elements(self, context) -> list:
        """
        Return all (Atom, Atom) pairs in the G-orbit of the canonical
        representative (+u_c, +v_c) of this PairFlavour.

        Applies every element of G = S_{n_1} × ... × S_{n_m} simultaneously
        to the canonical pair, collecting unique results via a seen-set.
        The length of the returned list equals orbit_size(type_sizes).

        Complexity: O(∏_g n_g!) — suitable for small test contexts only.
        For encoding, use DirectOrbitEnumerator which generates the (+,+)
        positive-sign subset directly without materialising the full orbit.

        Note on sign correlations
        -------------------------
        For ANTISYMMETRIC operations with non-zero overlap, signs of u and v
        are correlated by the simultaneous group action and NOT independent.
        In particular, for full-overlap pairs (NULL_SELF) any permutation
        acts identically on both atoms, so only (++, −−) sign combinations
        appear — not (+−) or (−+), which belong to a separate G-orbit.
        """
        from itertools import permutations as _perms

        type_sizes = tuple(g.size for g in context.types)
        if self.count(type_sizes) == 0:
            return []

        # Canonical representative: first valid label assignment.
        u_labels, v_labels = [], []
        for g, ku, kv, s in zip(context.types, self.flavour_u.counts,
                                  self.flavour_v.counts, self.overlap):
            u_labels.extend(g.labels[:ku])               # first ku labels for u
            v_labels.extend(g.labels[:s])                # shared labels for v
            v_labels.extend(g.labels[ku:ku + kv - s])   # v-only labels

        u_canon = Atom(self.op_u, tuple(u_labels), sign=+1)
        v_canon = Atom(self.op_v, tuple(v_labels), sign=+1)

        # Apply every σ ∈ G simultaneously to (u_canon, v_canon).
        seen = set()
        result = []
        for combo in _iproduct(*[list(_perms(g.labels)) for g in context.types]):
            perm_map = {}
            for g, perm in zip(context.types, combo):
                for orig, new in zip(g.labels, perm):
                    perm_map[orig] = new
            new_u = Atom(self.op_u,
                         tuple(perm_map[l] for l in u_canon.labels),
                         sign=u_canon.sign)
            new_v = Atom(self.op_v,
                         tuple(perm_map[l] for l in v_canon.labels),
                         sign=v_canon.sign)
            pair = (new_u, new_v)
            if pair not in seen:
                seen.add(pair)
                result.append(pair)
        return result


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
    for g in context.types:
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

    For each ordered pair (fo_u, fo_v) of FlavouredOperators — including
    fo_u == fo_v (self-pairing is intentional: pairing a FO with itself captures
    correlations between different atoms of the same type, and at full overlap
    it degenerates to the single-atom orbit which is handled by Phase 1 encoding)
    — the valid overlap range for group i is:
        s ∈ [max(0, k_u_i + k_v_i − n_i),  min(k_u_i, k_v_i)]

    PairFlavour's canonical ordering ensures (fo_u, fo_v) and (fo_v, fo_u)
    contribute identical PairFlavours; a seen-set handles deduplication.

    Sort order guarantee
    --------------------
    The returned list is sorted by:
        (op_u_key, op_v_key, overlap)
    where op_key = (name, rank, parity, argument_symmetry, flavour.counts).

    Consequence: all PairFlavours sharing the same (op_u, flavour_u, op_v,
    flavour_v) — i.e. an OVERLAP BLOCK — appear as a contiguous run in the
    output, ordered by ascending overlap tuple within the block.  Code that
    groups PairFlavours into OVERLAP BLOCKS (e.g. encode() and
    describe_encoding()) relies on this contiguity and must continue to do so.
    Any alternative implementation of canonical_pair_flavours must preserve it.
    """
    type_sizes = tuple(g.size for g in context.types)
    seen = set()
    for fo_u in fo_list:
        for fo_v in fo_list:
            per_type = []
            valid = True
            for n, ku, kv in zip(
                type_sizes, fo_u.flavour.counts, fo_v.flavour.counts
            ):
                s_lo = max(0, ku + kv - n)
                s_hi = min(ku, kv)
                if s_lo > s_hi:
                    valid = False
                    break
                per_type.append(range(s_lo, s_hi + 1))
            if not valid:
                continue
            for overlap in _iproduct(*per_type):
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
    all_atoms = [atom for fo in fo_list for atom in fo.atoms_one_per_sign()]
    return {
        pair_flavour_of(u, v, context)
        for u in all_atoms
        for v in all_atoms
    }
