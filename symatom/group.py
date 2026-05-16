"""
symatom.group — TheGroup: the discrete symmetry group acting on atoms and atom-pairs.

TheGroup encapsulates the group G that acts simultaneously on atoms and on pairs
of atoms.  In the common case G = S_{n_1} × ... × S_{n_m} (the direct product of
symmetric groups, one per species), but the class is designed for extension to
additional symmetries such as spatial parity (ParityOp), BeamExchange, or other
discrete operations that may be incorporated in future without breaking the
interface.

Key objects
-----------
TheGroup
    Frozen dataclass.  Constructed from a Context.  Provides:

    Single-atom methods (short names):
      orbit(u)                      — the G-orbit of a single atom
      stabiliser_size(u)            — |Stab_G(u)|  (algebraic, O(n))
      orbit_size(u)                 — |G| / |Stab_G(u)|  (algebraic, O(n))
      in_orbit(candidate, rep)      — single-atom membership test
      orbit_brute(u)                — same as orbit(), O(∏ n_g!)
      in_orbit_brute(candidate, rep)— same as in_orbit(), O(∏ n_g!)

    Atom-pair methods (_pair suffix):
      orbit_pair(u, v)              — the G-orbit of the atom-pair
      stabiliser_size_pair(u, v)    — |Stab_G(u, v)|  (algebraic, O(n))
      orbit_size_pair(u, v)         — |G| / |Stab_G(u, v)|  (algebraic, O(n))
      in_orbit_pair(candidate, rep) — pair membership test
      orbit_brute_pair(u, v)        — same as orbit_pair(), O(∏ n_g!)
      in_orbit_brute_pair(c, rep)   — same as in_orbit_pair(), O(∏ n_g!)

    Group structure:
      order()                       — |G|

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

    # ------------------------------------------------------------------
    # Internal machinery
    # ------------------------------------------------------------------

    def _apply(self, perm_map: dict, atom: Atom) -> Atom:
        """Return σ·atom for a given label-substitution map σ."""
        new_labels = tuple(perm_map[lbl] for lbl in atom.labels)
        return Atom(atom.operation, new_labels, sign=atom.sign)

    def _all_perm_maps(self):
        """Yield every σ ∈ G as a {label → label} substitution dict."""
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

        Every σ ∈ G is applied to u.  The Atom constructor handles label sorting
        and absorbs permutation parity into the sign for ANTISYMMETRIC operations,
        so the result is always a list of canonical (sorted-labels) atoms.

        Returns a list with no duplicates in the same order as first encountered
        while iterating over G.  The first element is always u itself (from the
        identity permutation).
        """
        seen: set[Atom] = set()
        result: list[Atom] = []
        for perm_map in self._all_perm_maps():
            a = self._apply(perm_map, u)
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

        Every σ ∈ G is applied simultaneously to both atoms.  The Atom
        constructor handles label sorting and absorbs permutation parity into
        the sign for ANTISYMMETRIC operations, so the result is always a
        canonical (sorted-labels) pair.

        Returns a list with no duplicates in the same order as they are first
        encountered while iterating over G.  The first element is always
        (u, v) itself (from the identity permutation).
        """
        seen: set[tuple[Atom, Atom]] = set()
        result: list[tuple[Atom, Atom]] = []
        for perm_map in self._all_perm_maps():
            pair = (self._apply(perm_map, u), self._apply(perm_map, v))
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
