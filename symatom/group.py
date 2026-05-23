"""
symatom.group — GroupElement and TheGroup: the discrete symmetry group acting on
atoms and atom-pairs.

TheGroup encapsulates the group G that acts simultaneously on atoms and on pairs
of atoms.  In the common case G = S_{n_1} × ... × S_{n_m} (the direct product of
symmetric groups, one per species), but the class is designed for extension to
additional symmetries such as spatial parity (ParityOp), BeamExchange, or other
discrete operations that may be incorporated in future without breaking the
interface.

Key objects
-----------
GroupElement
    One element of the discrete symmetry group.  Obtained from
    TheGroup.all_group_elements().  Internal representation is opaque to callers.
    Provides:
      apply(atom)            — return g·atom (a new Atom)
      apply_to_event(event)  — return a new event dict with labels permuted by g
      is_identity()          — True iff this is the identity element
      __repr__()             — cycle notation; identity prints as "1"
      __eq__, __hash__       — value semantics (two GroupElements are equal iff
                               they represent the same group action)

TheGroup
    Frozen dataclass.  Constructed from a Context.  Provides:

    Group structure:
      order()                       — |G|
      all_group_elements()          — all GroupElements; identity FIRST (contract)

    Single-atom methods (short names):
      orbit(u)                      — the G-orbit of a single atom
      stabiliser_size(u)            — |Stab_G(u)|  (algebraic, O(n))
      orbit_size(u)                 — |G| / |Stab_G(u)|  (algebraic, O(n))
      in_orbit(candidate, rep)      — single-atom membership test
      orbit_brute(u)                — same as orbit(), O(∏ n_g!)

    Atom-pair methods (_pair suffix):
      orbit_pair(u, v)              — the G-orbit of the atom-pair
      stabiliser_size_pair(u, v)    — |Stab_G(u, v)|  (algebraic, O(n))
      orbit_size_pair(u, v)         — |G| / |Stab_G(u, v)|  (algebraic, O(n))
      in_orbit_pair(candidate, rep) — pair membership test
      orbit_brute_pair(u, v)        — same as orbit_pair(), O(∏ n_g!)

    Atom-sequence methods (_sequence suffix):
      orbit_sequence(atoms)         — the G-orbit of an ordered sequence of atoms
      orbit_brute_sequence(atoms)   — same as orbit_sequence(), O(∏ n_g!)

    All *_brute methods are permanent O(∏ n_g!) reference implementations.
    They must not be removed; they serve as ground-truth for validating faster
    replacements.  Cross-validation tests in test_group.py verify that every
    primary method agrees with its *_brute counterpart.

Step history
------------
Step 2: DirectOrbitEnumerator rewritten to use TheGroup.orbit_pair() rather than
    the old OrbitUnion approach.
Step 3: orbit methods made algebraic O(n); *_brute variants added for
    all four orbit-related methods; in_orbit_pair/orbit_pair primary
    methods still delegate to their *_brute counterparts pending Step 4.
    sign_correlation_type() (formerly here) moved to symcoder/encode.py.
Step 4 (planned): replace orbit_pair(), in_orbit_pair() with O(orbit_size)
    direct combinatorial algorithms.
Step 5: GroupElement introduced as the public handle for one group action.
    _apply() removed from TheGroup; logic lives in GroupElement.apply().
    all_group_elements() added as the public enumeration API.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from itertools import permutations as _perms, product as _iproduct
from .atoms import Atom, ArgumentSymmetry, VectorType


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _E(k: int) -> int:
    """Number of even permutations of k elements: k!/2 for k>=2, else 1."""
    return math.factorial(k) // 2 if k >= 2 else 1


def _O(k: int) -> int:
    """Number of odd permutations of k elements: k!/2 for k>=2, else 0."""
    return math.factorial(k) // 2 if k >= 2 else 0


def _cycle_notation(perm_map: dict) -> str:
    """Return cycle notation for a label-permutation dict.

    Fixed points are omitted.  The identity (all fixed points) returns "1".
    Labels within each cycle are ordered by first encounter; cycles are
    ordered lexicographically by their first element.

    Examples
    --------
    {'a': 'b', 'b': 'a', 'p': 'p'} → "(a b)"
    {'a': 'a', 'b': 'b', 'p': 'p'} → "1"
    {'a': 'b', 'b': 'c', 'c': 'a'} → "(a b c)"
    """
    remaining = set(perm_map.keys())
    cycles = []
    for start in sorted(remaining):
        if start not in remaining:
            continue
        if perm_map[start] == start:
            remaining.discard(start)
            continue
        cycle = [start]
        remaining.discard(start)
        current = perm_map[start]
        while current != start:
            cycle.append(current)
            remaining.discard(current)
            current = perm_map[current]
        cycles.append("(" + " ".join(cycle) + ")")
    return "".join(cycles) if cycles else "1"


# ---------------------------------------------------------------------------
# GroupElement
# ---------------------------------------------------------------------------

class GroupElement:
    """
    One element of the discrete symmetry group G.

    A GroupElement encapsulates a single group action.  Its internal
    representation is private; callers should treat it as an opaque handle.
    GroupElements are obtained from TheGroup.all_group_elements() — TheGroup
    is the sole authority on what elements the group contains and how they act.

    At present, the internal representation is a label-permutation dict, but
    this is an implementation detail that may change (e.g. when parity flips or
    other discrete generators are added).  Callers must not rely on it.

    Methods
    -------
    apply(atom)           → Atom
        Return g·atom: the atom with its labels permuted (and sign adjusted for
        ANTISYMMETRIC operations) by this group element.

    apply_to_event(event) → dict
        Return a new event dict with particle labels permuted by this group
        element.  TheGroup is the sole authority on group action; callers must
        use this method rather than permuting event dicts by hand.

    is_identity()         → bool
        True iff this is the identity element e ∈ G.

    Comparison
    ----------
    Two GroupElements are equal iff they represent the same group action.
    GroupElements are hashable and may be used as dict keys or in sets.

    Representation
    --------------
    Printed in cycle notation (ANTISYMMETRIC-sign extensions will add a ± prefix
    when parity generators are introduced).  The identity prints as "1".
    """

    def __init__(self, perm_map: dict) -> None:
        # _perm_map is private.  Do not access it outside this class or TheGroup.
        self._perm_map: dict = perm_map

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def apply(self, atom: Atom) -> Atom:
        """Return g·atom: atom with labels permuted by this group element.

        For ANTISYMMETRIC operations the Atom constructor absorbs the
        permutation parity into the atom's sign automatically.
        """
        new_labels = tuple(self._perm_map[lbl] for lbl in atom.labels)
        return Atom(atom.operation, new_labels, sign=atom.sign)

    def apply_to_event(self, event: dict) -> dict:
        """Return a new event dict with particle labels permuted by this element.

        ``event`` maps label → vector (any array-like).  The returned dict maps
        ``g(label) → vector`` for every label in the event.

        TheGroup is the sole authority on what the group action is.  Callers
        must not permute event labels by hand; use this method so that future
        extensions (parity flips, beam exchange, …) are handled automatically.
        """
        return {self._perm_map[lbl]: vec for lbl, vec in event.items()}

    def is_identity(self) -> bool:
        """Return True iff this is the identity element (every label maps to itself)."""
        return all(v == k for k, v in self._perm_map.items())

    # ------------------------------------------------------------------
    # Python data-model
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return _cycle_notation(self._perm_map)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GroupElement):
            return NotImplemented
        return self._perm_map == other._perm_map

    def __hash__(self) -> int:
        return hash(frozenset(self._perm_map.items()))


# ---------------------------------------------------------------------------
# TheGroup
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TheGroup:
    """
    A discrete group acting on atoms and atom-pairs.

    In the common case the group is G = S_{n_1} × ... × S_{n_m}: each factor
    S_{n_g} permutes the labels of the g-th VectorType, and the permutation is
    applied simultaneously to one atom (single-atom methods) or to both atoms
    of a pair (_pair methods).  This is the physically natural symmetry for
    systems of indistinguishable particles.

    The class is designed to accommodate extensions without breaking the
    interface.  Future variants may add ParityOp (spatial inversion),
    BeamExchange (swapping beam-left and beam-right labels), or other discrete
    generators as optional constructor arguments.  The current implementation
    realises only the base permutation group.

    All methods are O(∏ n_g!) brute-force.  For n_g ~ 30 these will be
    replaced by algebraic shortcut methods (Step 3 of the architecture
    redesign).  The brute-force implementations serve as the permanent
    reference against which fast replacements will be validated.

    Atom-sequence methods act on an ordered list of atoms simultaneously,
    generalising orbit_brute_pair to arbitrary length.  Used for ground-truth
    computation of the G-orbit of repS evaluations in alignment-decoder tests.

    Parameters
    ----------
    types : tuple[VectorType]
        The VectorTypes whose labels this group acts on.  Each VectorType
        contributes one symmetric-group factor S_{n_i}, permuting its labels.
        TheGroup stores the label sets directly rather than a full Context so
        that Context can hold a TheGroup without a circular dependency.

    Long-run design note
    --------------------
    The algebraic methods (stabiliser_size, order) only need the per-type
    label *sets* to partition an atom's labels into per-group counts.  In
    principle those methods are label-agnostic: they could work from
    pre-computed per-type counts without storing VectorType objects at all.
    The brute-force enumeration methods (_all_perm_maps, orbit_brute) need the
    actual label sequences to generate concrete permutations.  A future step
    could decouple these two concerns and reduce TheGroup to TheGroup(sizes).
    """
    types: tuple  # tuple[VectorType]

    # ------------------------------------------------------------------
    # Group structure
    # ------------------------------------------------------------------

    def order(self) -> int:
        """Return |G| = ∏_g n_g!  (total number of group elements)."""
        result = 1
        for g in self.types:
            result *= math.factorial(g.size)
        return result

    def all_group_elements(self) -> tuple:
        """Return all elements of the group as GroupElement objects.

        **Contract**: the identity element is always first.  This is a
        guaranteed API promise, not an implementation detail.  Callers may
        rely on ``all_group_elements()[0].is_identity()`` always being True.

        The identity-first ordering is a structural consequence of how
        _all_perm_maps() iterates (see its docstring).  The assert below
        guards against any future change to _all_perm_maps() that would
        silently break the contract; it is compiled away under ``python -O``.

        TheGroup is the sole authority on group membership and ordering.
        """
        result = tuple(GroupElement(pm) for pm in self._all_perm_maps())
        assert result[0].is_identity(), (
            "_all_perm_maps() must yield the identity first — "
            "all_group_elements() identity-first contract violated"
        )
        return result

    # ------------------------------------------------------------------
    # Internal machinery
    # ------------------------------------------------------------------

    def _all_perm_maps(self):
        """Yield every σ ∈ G as a {label → label} substitution dict.

        Private.  Used only inside all_group_elements() to construct
        GroupElement objects.  All other code should iterate via
        all_group_elements() and call g.apply() / g.apply_to_event().

        Implementation note: the identity permutation is always yielded
        first.  This is a structural guarantee: itertools.permutations()
        starts with the natural (identity) ordering, and itertools.product()
        of sequences that each begin with their identity element produces the
        all-identity combination first.  all_group_elements() relies on this.
        """
        group_perm_lists = [list(_perms(g.labels)) for g in self.types]
        for combo in _iproduct(*group_perm_lists):
            perm_map = {}
            for g, perm in zip(self.types, combo):
                for orig, new in zip(g.labels, perm):
                    perm_map[orig] = new
            yield perm_map

    # ------------------------------------------------------------------
    # Single-atom orbit and stabiliser
    # ------------------------------------------------------------------

    def orbit(self, u: Atom) -> list[Atom]:
        """
        Return all distinct Atoms in the G-orbit of u.

        Currently delegates to orbit_brute().  A direct O(orbit_size)
        combinatorial implementation is planned for Step 4.   TODO !!
        """
        return self.orbit_brute(u)

    def orbit_brute(self, u: Atom) -> list[Atom]:
        """
        Brute-force O(∏ n_g!) single-atom orbit enumeration.  Permanent reference.

        Every g ∈ G is applied to u.  The Atom constructor handles label sorting
        and absorbs permutation parity into the sign for ANTISYMMETRIC operations,
        so the result is always a list of canonical (sorted-labels) atoms.

        Returns a list with no duplicates in the same order as first encountered
        while iterating over G.  The first element is always u itself (from the
        identity element, which all_group_elements() guarantees comes first).
        """
        seen: set[Atom] = set()
        result: list[Atom] = []
        for g in self.all_group_elements():
            a = g.apply(u)
            if a not in seen:
                seen.add(a)
                result.append(a)
        return result

    def stabiliser_size(self, u: Atom) -> int:
        """
        Return |Stab_G(u)| — the number of σ ∈ G with σ·u = u.

        Computed algebraically.  For each group g let ku_g be the number of
        u's labels drawn from g, and r = n_g - ku_g the remaining labels.

          SYMMETRIC:     stab_g = ku_g! · r!
          ANTISYMMETRIC: stab_g = E(ku_g) · r!

        For ANTISYMMETRIC atoms only even permutations of u's labels within
        each group preserve the sign, so only those elements fix u.
        The total stabiliser size is the product across all groups.
        """
        antisym = u.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        u_set = set(u.labels)
        stab = 1
        for g in self.types:
            ku = len(u_set & set(g.labels))
            r = g.size - ku
            stab_g = (_E(ku) if antisym else math.factorial(ku)) * math.factorial(r)
            stab *= stab_g
        return stab

    def orbit_size(self, u: Atom) -> int:
        """Return |orbit(u)| = |G| / |Stab_G(u)|."""
        return self.order() // self.stabiliser_size(u)

    # ------------------------------------------------------------------
    # Single-atom membership
    # ------------------------------------------------------------------

    def in_orbit(self, candidate: Atom, representative: Atom) -> bool:
        """
        Return True iff candidate ∈ G · representative.

        Currently delegates to in_orbit_brute().
        """
        return self.in_orbit_brute(candidate, representative)

    def in_orbit_brute(self, candidate: Atom, representative: Atom) -> bool:
        """
        Brute-force O(∏ n_g!) single-atom membership test.  Permanent reference.
        """
        return candidate in self.orbit_brute(representative)

    # ------------------------------------------------------------------
    # Atom-pair orbit and stabiliser
    # ------------------------------------------------------------------

    def orbit_pair(self, u: Atom, v: Atom) -> list[tuple[Atom, Atom]]:
        """
        Return all distinct (Atom, Atom) pairs in the G-orbit of (u, v).

        Currently delegates to orbit_brute_pair().  A direct O(orbit_size)
        combinatorial implementation is planned for Step 4.   TODO !!
        """
        return self.orbit_brute_pair(u, v)

    def orbit_brute_pair(self, u: Atom, v: Atom) -> list[tuple[Atom, Atom]]:
        """
        Brute-force O(∏ n_g!) pair-orbit enumeration.  Permanent reference.

        Every g ∈ G is applied simultaneously to both atoms.  The Atom
        constructor handles label sorting and absorbs permutation parity into
        the sign for ANTISYMMETRIC operations, so the result is always a
        canonical (sorted-labels) pair.

        Returns a list with no duplicates in the same order as they are first
        encountered while iterating over G.  The first element is always
        (u, v) itself (from the identity element).
        """
        seen: set[tuple[Atom, Atom]] = set()
        result: list[tuple[Atom, Atom]] = []
        for g in self.all_group_elements():
            pair = (g.apply(u), g.apply(v))
            if pair not in seen:
                seen.add(pair)
                result.append(pair)
        return result

    def stabiliser_size_pair(self, u: Atom, v: Atom) -> int:
        """
        Return |Stab_G(u, v)| — the number of σ ∈ G with σ·(u, v) = (u, v).

        Computed algebraically using the E(k)/O(k) formula.  For each group g
        the labels used by u and v in that group are partitioned into
        (shared s, u-only a, v-only b, unused r).  The per-group stabiliser
        count depends on whether each operation is ANTISYMMETRIC:

          SS:  s! · a! · b! · r!
          AS:  (E(s)·E(a) + O(s)·O(a)) · b! · r!
          SA:  (E(s)·E(b) + O(s)·O(b)) · a! · r!
          AA:  (E(s)·E(a)·E(b) + O(s)·O(a)·O(b)) · r!

        The total stabiliser size is the product across all groups.
        """
        antisym_u = u.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        antisym_v = v.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        u_label_set = set(u.labels)
        v_label_set = set(v.labels)
        stab = 1
        for g in self.types:
            g_set = set(g.labels)
            u_g = u_label_set & g_set
            v_g = v_label_set & g_set
            shared_g = u_g & v_g
            u_only_g = u_g - shared_g
            v_only_g = v_g - shared_g
            s = len(shared_g)
            a = len(u_only_g)
            b = len(v_only_g)
            r = g.size - s - a - b
            if antisym_u and antisym_v:
                stab_g = (
                    _E(s) * _E(a) * _E(b) + _O(s) * _O(a) * _O(b)
                ) * math.factorial(r)
            elif antisym_u and not antisym_v:
                stab_g = (
                    _E(s) * _E(a) + _O(s) * _O(a)
                ) * math.factorial(b) * math.factorial(r)
            elif not antisym_u and antisym_v:
                stab_g = (
                    _E(s) * _E(b) + _O(s) * _O(b)
                ) * math.factorial(a) * math.factorial(r)
            else:  # SS
                stab_g = (
                    math.factorial(s) * math.factorial(a)
                    * math.factorial(b) * math.factorial(r)
                )
            stab *= stab_g
        return stab

    def orbit_size_pair(self, u: Atom, v: Atom) -> int:
        """Return |orbit_pair(u, v)| = |G| / |Stab_G(u, v)|."""
        return self.order() // self.stabiliser_size_pair(u, v)

    # ------------------------------------------------------------------
    # Atom-pair membership
    # ------------------------------------------------------------------

    def in_orbit_pair(
        self,
        candidate: tuple[Atom, Atom],
        representative: tuple[Atom, Atom],
    ) -> bool:
        """
        Return True iff candidate ∈ G · representative.

        Currently delegates to in_orbit_brute_pair().  An algebraic O(n)
        implementation is planned for Step 4.
        """
        return self.in_orbit_brute_pair(candidate, representative)

    def in_orbit_brute_pair(
        self,
        candidate: tuple[Atom, Atom],
        representative: tuple[Atom, Atom],
    ) -> bool:
        """
        Brute-force O(∏ n_g!) pair membership test.  Permanent reference.

        Equivalent to asking whether candidate and representative are in the
        same G-orbit.
        """
        orb_set = self.orbit_brute_pair(representative[0], representative[1])
        return candidate in orb_set

    # ------------------------------------------------------------------
    # Atom-sequence orbit
    # ------------------------------------------------------------------

    def orbit_sequence(self, atoms: list) -> list:
        """
        Return all distinct tuples in the G-orbit of an ordered sequence of atoms.

        Currently delegates to orbit_brute_sequence().
        """
        return self.orbit_brute_sequence(atoms)

    def orbit_brute_sequence(self, atoms: list) -> list:
        """
        Brute-force O(∏ n_g!) orbit enumeration for an ordered sequence of atoms.
        Permanent reference.

        Every g ∈ G is applied simultaneously to all atoms in the sequence.
        The Atom constructor handles label sorting and sign absorption for
        ANTISYMMETRIC operations, so all results are in canonical form.

        Returns a list of tuples with no duplicates, in the same order as first
        encountered while iterating over G.  The first element is always
        tuple(atoms) itself (from the identity element).

        This is the natural generalisation of orbit_brute_pair to sequences of
        arbitrary length.  TheGroup is the sole authority on orbit construction;
        callers must not attempt to build orbits of sequences by other means.
        """
        seen: set = set()
        result: list = []
        for g in self.all_group_elements():
            seq = tuple(g.apply(a) for a in atoms)
            if seq not in seen:
                seen.add(seq)
                result.append(seq)
        return result
