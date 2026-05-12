"""
symatom.group — TheGroup: the discrete symmetry group acting on atom-pairs.

TheGroup encapsulates the group G that acts simultaneously on both atoms of a
pair.  In the common case G = S_{n_1} × ... × S_{n_m} (the direct product of
symmetric groups, one per species), but the class is designed for extension to
additional symmetries such as spatial parity (ParityOp), BeamExchange, or other
discrete operations that may be incorporated in future without breaking the
interface.

Key objects
-----------
SignCorrelationType
    Enum describing how the G-orbit of (u, v) relates the sign degrees of
    freedom of u and v.  Needed for compressed polynomial encodings and for
    determining which sign-variant orbits are genuinely distinct.

TheGroup
    Frozen dataclass.  Constructed from a Context.  Provides:
      order()                     — |G|
      orbit(u, v)                 — the G-orbit of the atom-pair
      stabiliser_size(u, v)       — |Stab_G(u, v)|
      orbit_size(u, v)            — |G| / |Stab_G(u, v)|
      in_orbit(candidate, rep)    — membership test
      sign_correlation_type(u, v) — how signs are coupled inside the orbit
      companion_orbits(u, v)      — all distinct orbits related by sign flips
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from enum import Enum
from itertools import permutations as _perms, product as _iproduct
from .atoms import Atom, ArgumentSymmetry
from .context import Context


# ---------------------------------------------------------------------------
# SignCorrelationType
# ---------------------------------------------------------------------------

class SignCorrelationType(Enum):
    """
    Describes how the G-orbit of an atom-pair (u, v) relates the sign degrees
    of freedom of u and v.

    The label TYPE_XY encodes:
      X = 1  — u's sign does not flip independently within the orbit
      X = 2  — u's sign can flip independently of v's sign within the orbit
      Y = 1  — v's sign does not flip independently within the orbit
      Y = 2  — v's sign can flip independently of u's sign within the orbit

    "Independently" here means: there exists an orbit element where that
    atom's sign differs from the canonical representative but the *other*
    atom's sign is unchanged.

    Representative cases (for G = S_n acting on a single group):

    TYPE_11
        Neither sign flips independently.  This covers two structurally
        different situations:

        (a) Both operations are SYMMETRIC — signs are always +1, nothing
            can flip at all.

        (b) Both ANTISYMMETRIC with full label overlap (u and v share all
            their labels).  Every permutation flips both signs by the same
            amount, so only (u, v) and (-u, -v) appear in the orbit; the
            mixed-sign pairs (-u, v) and (u, -v) form a *separate* orbit.

    TYPE_12
        Only v's sign flips freely.  Typical: u SYMMETRIC, v ANTISYMMETRIC.

    TYPE_21
        Only u's sign flips freely.  Typical: u ANTISYMMETRIC, v SYMMETRIC.

    TYPE_22
        Both signs flip freely (are independent).  Typical: both
        ANTISYMMETRIC with disjoint label sets (zero overlap).

    Important: the type is NOT determined by ArgumentSymmetry alone.  For
    ANTISYMMETRIC × ANTISYMMETRIC pairs it depends on the overlap structure
    (i.e. on which concrete labels u and v share), and in principle on the
    full group action, ParityBehaviour, etc.  TheGroup.sign_correlation_type
    always computes it from the actual orbit.
    """
    TYPE_11 = "11"
    TYPE_12 = "12"
    TYPE_21 = "21"
    TYPE_22 = "22"


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
    A discrete group acting on atom-pairs.

    In the common case the group is G = S_{n_1} × ... × S_{n_m}: each factor
    S_{n_g} permutes the labels of the g-th VectorGroup, and the permutation is
    applied *simultaneously* to both atoms of a pair.  This is the physically
    natural symmetry for systems of indistinguishable particles.

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
    context : Context
        Defines the groups and their label sets.  TheGroup acts on atoms whose
        labels are drawn from context.all_labels.
    """
    context: Context

    # ------------------------------------------------------------------
    # Group structure
    # ------------------------------------------------------------------

    def order(self) -> int:
        """Return |G| = ∏_g n_g!  (total number of group elements)."""
        result = 1
        for g in self.context.groups:
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
        group_perm_lists = [list(_perms(g.labels)) for g in self.context.groups]
        for combo in _iproduct(*group_perm_lists):
            perm_map = {}
            for g, perm in zip(self.context.groups, combo):
                for orig, new in zip(g.labels, perm):
                    perm_map[orig] = new
            yield perm_map

    # ------------------------------------------------------------------
    # Orbit and stabiliser
    # ------------------------------------------------------------------

    def orbit(self, u: Atom, v: Atom) -> list[tuple[Atom, Atom]]:
        """
        Return all distinct (Atom, Atom) pairs in the G-orbit of (u, v).

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

    def stabiliser_size(self, u: Atom, v: Atom) -> int:
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
        for g in self.context.groups:
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

    def orbit_size(self, u: Atom, v: Atom) -> int:
        """Return |orbit(u, v)| = |G| / |Stab_G(u, v)|."""
        return self.order() // self.stabiliser_size(u, v)

    # ------------------------------------------------------------------
    # Membership
    # ------------------------------------------------------------------

    def in_orbit(
        self,
        candidate: tuple[Atom, Atom],
        representative: tuple[Atom, Atom],
    ) -> bool:
        """
        Return True iff candidate ∈ G · representative.

        Equivalent to asking whether candidate and representative are in the
        same G-orbit.  Currently O(|G|) brute-force; algebraic shortcut
        is planned for Step 3.
        """
        orb_set = set(self.orbit(representative[0], representative[1]))
        return candidate in orb_set

    # ------------------------------------------------------------------
    # Sign structure
    # ------------------------------------------------------------------

    def sign_correlation_type(self, u: Atom, v: Atom) -> SignCorrelationType:
        """
        Determine how the G-orbit of (u, v) relates the sign degrees of freedom.

        A sign "flips freely" for one atom iff there exists an orbit element
        where that atom's sign differs from the canonical representative while
        the other atom's sign is *unchanged*.  This is the empirical definition
        — it is computed from the actual orbit, not assumed from ArgumentSymmetry.

        For SYMMETRIC operations the sign is always +1 and can never flip, so
        those atoms always contribute a "1" to the type label.

        See SignCorrelationType for a full description of each type.
        """
        orb = self.orbit(u, v)

        # u_flips_freely: ∃ (a, b) in orbit where a.sign ≠ u.sign AND b.sign == v.sign
        u_flips_freely = any(
            a.sign != u.sign and b.sign == v.sign
            for a, b in orb
        )
        # v_flips_freely: ∃ (a, b) in orbit where b.sign ≠ v.sign AND a.sign == u.sign
        v_flips_freely = any(
            b.sign != v.sign and a.sign == u.sign
            for a, b in orb
        )

        if u_flips_freely and v_flips_freely:
            return SignCorrelationType.TYPE_22
        elif u_flips_freely:
            return SignCorrelationType.TYPE_21
        elif v_flips_freely:
            return SignCorrelationType.TYPE_12
        else:
            return SignCorrelationType.TYPE_11

    # ------------------------------------------------------------------
    # Companion orbits
    # ------------------------------------------------------------------

    def companion_orbits(
        self, u: Atom, v: Atom
    ) -> list[list[tuple[Atom, Atom]]]:
        """
        Return all distinct G-orbits that are sign-related to orbit(u, v).

        The sign-related variants of (u, v) are:
          (u,   v)   — always included
          (-u,  v)   — only if u.operation is ANTISYMMETRIC
          (u,  -v)   — only if v.operation is ANTISYMMETRIC
          (-u, -v)   — only if both are ANTISYMMETRIC

        Each variant's orbit is computed and compared to those already seen.
        G-orbits partition the space, so two orbits are identical iff they
        share any element.

        Returns a list of unique orbits (each a list of (Atom, Atom) pairs),
        with the primary orbit(u, v) as the first entry.
        """
        antisym_u = u.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        antisym_v = v.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

        # Build the sign variants to consider
        sign_variants: list[tuple[Atom, Atom]] = [(u, v)]
        if antisym_u:
            neg_u = Atom(u.operation, u.labels, -u.sign)
            sign_variants.append((neg_u, v))
        if antisym_v:
            neg_v = Atom(v.operation, v.labels, -v.sign)
            sign_variants.append((u, neg_v))
        if antisym_u and antisym_v:
            neg_u = Atom(u.operation, u.labels, -u.sign)
            neg_v = Atom(v.operation, v.labels, -v.sign)
            sign_variants.append((neg_u, neg_v))

        # Collect unique orbits (orbits are either identical or disjoint)
        seen_orbit_sets: list[frozenset[tuple[Atom, Atom]]] = []
        result: list[list[tuple[Atom, Atom]]] = []
        for su, sv in sign_variants:
            orb = self.orbit(su, sv)
            orb_key = frozenset(orb)
            if orb_key not in seen_orbit_sets:
                seen_orbit_sets.append(orb_key)
                result.append(orb)
        return result
