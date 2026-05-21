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

        Note that atom_one_per_sign does NOT (and is NOT intended to) generate an orbit.
        E.g. it could in principle yield:

            [+eps2(a,b), +eps2(a,c), +eps2(b,c)]

        which is not an orbit, but it could also yield:
        
            [+dot(a,b), +dot(a,c), +dot(b,c)]

        which is an orbit! It can even yield stranger things with mixes signs like:

            [ -eps2(apple, zebra), -eps2(toast, zebra), +eps2(apple, toast)]

        For example, the demo program below:

        ###################################################################################
        def demo_scrambled_labels():
    
            # Show that atoms_one_per_sign() can yield atoms with MIXED signs when the
            # VectorType's label tuple is not in sorted order.

            # Root cause: combinations() preserves the declaration order of g.labels.
            # If that order is not sorted, the Atom constructor's label-sorting step
            # requires an odd permutation for some choices, flipping their sign to -1.

            # Example: labels=("zebra","apple","toast") with a rank-2 ANTISYMMETRIC op.
            # combinations(("zebra","apple","toast"), 2) yields:
            #     ("zebra","apple")  →  Atom sorts to ("apple","zebra"), 1 swap → sign=-1
            #     ("zebra","toast")  →  Atom sorts to ("toast","zebra"), 1 swap → sign=-1
            #     ("apple","toast")  →  already sorted                         → sign=+1
            # Output: a mix of + and - atoms despite all being passed sign=+1.
    
            eps2 = Operation("eps2", rank=2, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)
            scrambled = VectorType("scrambled", labels=("zebra", "apple", "toast"))
            ctx_s = Context(types=(scrambled,))
            fo = FlavouredOperator(operation=eps2, flavour=Flavour((2,)), context=ctx_s)

            _section("atoms_one_per_sign with scrambled labels — mixed sign output")
            print("  VectorType labels (declaration order): ('zebra', 'apple', 'toast')")
            print("  combinations yield pairs in that order, not in sorted order.")
            print()
            atoms = list(fo.atoms_one_per_sign())
            for atom in atoms:
                print(f"  {atom!r}")
            signs = [a.sign for a in atoms]
            print(f"\n  signs: {signs}  ← mix of +1 and -1 despite all constructed with sign=+1")
        ###################################################################################

        yields as output:

        ###################################################################################
        ============================================================
        atoms_one_per_sign with scrambled labels — mixed sign output
        ============================================================
        VectorType labels (declaration order): ('zebra', 'apple', 'toast')
        combinations yield pairs in that order, not in sorted order.

        -eps2(apple, zebra)
        -eps2(toast, zebra)
        +eps2(apple, toast)

        signs: [-1, -1, 1]  ← mix of +1 and -1 despite all constructed with sign=+1
        ###################################################################################

        TODO: Think about (1) whether it would be better for this algorithm to simply make
        operators all with +signs, (it could use the Atom.pinned_sign constructor), and/or
        (2) whether it should lean more on tools in TheGroup even though it is NOT a
        genrator of orbits. Thought (1) may be easy if nothing relies upon the historically
        strange sign conventions seen in the output of this method.
        """

        # For each vector type g and its required label count k, enumerate all
        # ways to choose k labels from g.labels using combinations (so order
        # within each group follows the declaration order of g.labels, not
        # sorted order).  The result is one list-of-tuples per type.
        per_type = [
            list(_combinations(g.labels, k))
            for g, k in zip(self.context.types, self.flavour)
        ]
        # Take the Cartesian product across all types: one choice per type per
        # iteration.  Flatten each multi-type choice into a single label tuple
        # by concatenating the per-type label tuples in type order.
        for combo_parts in _iproduct(*per_type):
            labels = tuple(lbl for part in combo_parts for lbl in part)
            # Construct with sign=+1.  The Atom constructor will sort the labels
            # and, for ANTISYMMETRIC operations, multiply sign by the parity of
            # the sorting permutation.  If the labels arrived in a non-sorted
            # order (possible when g.labels itself is not sorted), the stored
            # sign may end up -1.  See demo_scrambled_labels() in demo.py.
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

    def orbit_size(self, context) -> int:
        """
        Number of (Atom, Atom) pairs in the G-orbit of the canonical
        representative of this PairFlavour.

        Delegates entirely to TheGroup, which is the sole authority on orbit
        sizes.  TheGroup.orbit_size_pair() uses the orbit-stabiliser theorem
        |orbit| = |G| / |Stab(u, v)| and computes the stabiliser algebraically
        without inspecting argument_symmetry directly.
        """
        type_sizes = tuple(g.size for g in context.types)
        if self.count(type_sizes) == 0:
            return 0
        u_can, v_can = self.canonical_pair(context)
        return context.the_group.orbit_size_pair(u_can, v_can)

    def canonical_pair(self, context) -> tuple:
        """
        Return the canonical representative atom-pair (u_canon, v_canon) for
        this PairFlavour in the given context.

        Both atoms are constructed with sign=+1; the Atom constructor handles
        label sorting and absorbs permutation parity into the sign for
        ANTISYMMETRIC operations.

        This is the single definition of the canonical representative used by
        orbit_elements() and by plan.orbit_enumerator implementations.  It
        lives here — on PairFlavour — because the label-selection rule is a
        property of the PairFlavour structure, not of the group.

        Returns a 2-tuple (u_canon, v_canon).  Callers should pass these atoms
        directly to TheGroup methods for orbit computation.
        """
        u_labels, v_labels = [], []
        for g, ku, kv, s in zip(context.types, self.flavour_u.counts,
                                  self.flavour_v.counts, self.overlap):
            u_labels.extend(g.labels[:ku])               # first ku labels for u
            v_labels.extend(g.labels[:s])                # shared labels for v
            v_labels.extend(g.labels[ku:ku + kv - s])   # v-only labels
        return (
            Atom(self.op_u, tuple(u_labels), sign=+1),
            Atom(self.op_v, tuple(v_labels), sign=+1),
        )

    def orbit_elements(self, context) -> list:
        """
        Return all (Atom, Atom) pairs in the G-orbit of the canonical
        representative of this PairFlavour.

        Delegates to context.the_group.orbit_brute_pair(), which is the
        single authoritative brute-force pair-orbit implementation inside
        TheGroup.  All orbit logic — permutation enumeration, sign absorption,
        sign-correlation notes — lives in TheGroup.orbit_brute_pair(); see
        that method for full documentation.

        Complexity: O(∏_g n_g!) — suitable for small contexts only.
        For production encoding use plan.orbit_enumerator.orbit_elements(),
        which may be backed by a faster implementation.
        """
        type_sizes = tuple(g.size for g in context.types)
        if self.count(type_sizes) == 0:
            return []
        u_canon, v_canon = self.canonical_pair(context)
        return context.the_group.orbit_brute_pair(u_canon, v_canon)


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
