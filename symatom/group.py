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

    Primary methods (use these):
      order()                       — |G|
      orbit(u, v)                   — the G-orbit of the atom-pair
      stabiliser_size(u, v)         — |Stab_G(u, v)|  (algebraic, O(n))
      orbit_size(u, v)              — |G| / |Stab_G(u, v)|  (algebraic, O(n))
      in_orbit(candidate, rep)      — membership test
      sign_correlation_type(u, v)   — how signs are coupled  (algebraic, O(n))

    Brute-force reference methods (*_brute suffix):
      orbit_brute(u, v)             — same as orbit(), O(∏ n_g!)
      in_orbit_brute(candidate, rep)— same as in_orbit(), O(∏ n_g!)
      sign_correlation_type_brute(u, v) — same as sign_correlation_type(), O(∏ n_g!)

    All *_brute methods are permanent O(∏ n_g!) reference implementations.
    They must not be removed; they serve as ground-truth for validating faster
    replacements.  Cross-validation tests in test_group.py verify that every
    primary method agrees with its *_brute counterpart.

Step history
------------
Step 2: DirectOrbitEnumerator rewritten to use TheGroup.orbit() rather than the
    old OrbitUnion approach.
Step 3: sign_correlation_type() made algebraic O(n); *_brute variants added for
    all four orbit-related methods; in_orbit/orbit primary
    methods still delegate to their *_brute counterparts pending Step 4.
Step 4 (planned): replace orbit(), in_orbit() with O(orbit_size)
    direct combinatorial algorithms.
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

    The achievable set is the subgroup of Z_2 × Z_2 consisting of all
    (sign_u, sign_v) pairs that appear in the orbit.  There are exactly five
    subgroups of Z_2 × Z_2, corresponding to the five types below.

    The label TYPE_XY encodes:
      X = 1  — u's sign does not flip independently within the orbit
      X = 2  — u's sign can flip independently of v's sign within the orbit
      Y = 1  — v's sign does not flip independently within the orbit
      Y = 2  — v's sign can flip independently of u's sign within the orbit

    "Independently" means: ∃ orbit element where that atom's sign differs from
    the canonical representative but the *other* atom's sign is unchanged.

    TYPE_11
        Achievable = {(+1,+1)}.
        Neither sign can change at all.  Typical: both operations SYMMETRIC,
        or an ANTISYMMETRIC cross-group atom (all labels from different groups,
        so inter-group sorting never changes under G).

        Orbit structure: {z_k}  (n elements)
        Embedding: embed z_pos directly.

    TYPE_NEG
        Achievable = {(+1,+1), (-1,-1)}.
        Neither sign flips independently; only the correlated (simultaneous)
        negation of both signs is achievable.  Typical: both ANTISYMMETRIC
        with full per-group overlap in every group where both have ≥2 labels,
        e.g. AA pair where u and v share all of their labels within each group.

        Orbit structure: {z_k, -z_k}  (2n elements, closed under z → -z)
        Invariant: {z_k²}  (n complex values, since (-z_k)² = z_k²)
        Embedding: embed z_pos² directly.

        Note: the achievable set {(+1,+1), (-1,-1)} is the unique non-trivial
        subgroup of Z_2 × Z_2 that contains the correlated element but no
        independent flip.  It is distinct from TYPE_11 (no flip at all) and
        requires a different polynomial embedding.

    TYPE_12
        Achievable = {(+1,+1), (+1,-1)}.
        Only v's sign flips freely.  Typical: u SYMMETRIC, v ANTISYMMETRIC
        with ≥2 labels from the same group.

        Orbit structure: {z_k, conj(z_k)}  (2n elements, conjugate-closed)
        Embedding: embed z_full={z_pos, conj(z_pos)}, extract real polynomial
        coefficients.

    TYPE_21
        Achievable = {(+1,+1), (-1,+1)}.
        Only u's sign flips freely.  Typical: u ANTISYMMETRIC, v SYMMETRIC.

        Orbit structure: {z_k, -conj(z_k)}  (2n elements, anti-conj-closed)
        Embedding: embed z_full={z_pos, -conj(z_pos)}, extract even-real/
        odd-imaginary polynomial coefficients.

    TYPE_22
        Achievable = {(+1,+1), (-1,+1), (+1,-1), (-1,-1)}.
        Both signs flip freely (independent).  Typical: both ANTISYMMETRIC
        with different label sets in at least one group.

        Orbit structure: {z_k, conj(z_k), -z_k, -conj(z_k)}  (4n elements)
        Invariant: {z_k², conj(z_k²)}  (2n values, conjugate-closed)
        Embedding: embed z_sq_full={z_pos², conj(z_pos²)}, extract real
        polynomial coefficients.

    CRITICAL: the type is NOT determined by ArgumentSymmetry alone.
    SA pairs are NOT always TYPE_12 — a cross-group ANTISYMMETRIC atom
    (≤1 label per group) has a permanently +1 sign, giving TYPE_11 instead.
    Similarly AS pairs can be TYPE_11, and AA pairs can be TYPE_11, TYPE_NEG,
    or TYPE_22 depending on the per-group label overlap structure.
    Always compute the type via TheGroup.sign_correlation_type(), never by
    inspecting op.argument_symmetry directly.
    """
    TYPE_11  = "11"
    TYPE_NEG = "NEG"
    TYPE_12  = "12"
    TYPE_21  = "21"
    TYPE_22  = "22"


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

        Currently delegates to orbit_brute().  A direct O(orbit_size)
        combinatorial implementation is planned for Step 4.
        """
        return self.orbit_brute(u, v)

    def orbit_brute(self, u: Atom, v: Atom) -> list[tuple[Atom, Atom]]:
        """
        Brute-force O(∏ n_g!) orbit enumeration.  Permanent reference.

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

        Currently delegates to in_orbit_brute().  An algebraic O(n)
        implementation is planned for Step 4.
        """
        return self.in_orbit_brute(candidate, representative)

    def in_orbit_brute(
        self,
        candidate: tuple[Atom, Atom],
        representative: tuple[Atom, Atom],
    ) -> bool:
        """
        Brute-force O(∏ n_g!) membership test.  Permanent reference.

        Equivalent to asking whether candidate and representative are in the
        same G-orbit.
        """
        orb_set = set(self.orbit_brute(representative[0], representative[1]))
        return candidate in orb_set

    # ------------------------------------------------------------------
    # Sign structure
    # ------------------------------------------------------------------

    def sign_correlation_type(self, u: Atom, v: Atom) -> SignCorrelationType:
        """
        Algebraic O(n) determination of how the G-orbit couples signs.

        Algorithm — per-group achievable-set product
        --------------------------------------------
        The overall sign of an ANTISYMMETRIC atom is the product of per-group
        parities: sign = ∏_g parity_g, where parity_g is the parity of
        sorting the atom's labels from group g into ascending order.  Because
        labels from different groups never interleave (each group's labels sort
        as a contiguous block), different groups contribute independently.

        For each group g we compute the set of achievable (parity_u_g, parity_v_g)
        pairs under all σ_g ∈ S_{n_g}:

        • parity_u_g can only be −1 if u is ANTISYMMETRIC AND ku_g ≥ 2.
        • parity_v_g can only be −1 if v is ANTISYMMETRIC AND kv_g ≥ 2.
        • If neither can flip: achievable_g = {(+1, +1)}.
        • If only u can flip:  achievable_g = {(+1,+1), (−1,+1)}.
        • If only v can flip:  achievable_g = {(+1,+1), (+1,−1)}.
        • If both can flip AND u_g == v_g (full overlap in g):
              achievable_g = {(+1,+1), (−1,−1)}  — signs coupled.
        • If both can flip AND u_g ≠ v_g (partial/zero overlap in g):
              achievable_g = {(+1,+1), (+1,−1), (−1,+1), (−1,−1)} — independent.

        The overall achievable set is the group-product of per-group sets
        (one factor per group, independent choices across groups).

        Sign-correlation type is then read from the overall set:
          (−1,+1) ∈ achievable ↔ u_flips_freely  → bit X of the type label
          (+1,−1) ∈ achievable ↔ v_flips_freely  → bit Y of the type label

        Cross-validation: sign_correlation_type_brute() computes the same
        value by inspecting the actual orbit; both must always agree.
        """
        antisym_u = u.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        antisym_v = v.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

        u_set = set(u.labels)
        v_set = set(v.labels)

        # achievable: subset of {(+1,+1),(+1,-1),(-1,+1),(-1,-1)}
        # Represented as a Python set; start with the identity element.
        achievable: set[tuple[int, int]] = {(1, 1)}

        for g in self.context.groups:
            g_set = set(g.labels)
            u_g = u_set & g_set
            v_g = v_set & g_set
            ku_g = len(u_g)
            kv_g = len(v_g)

            pu_can_flip = antisym_u and ku_g >= 2
            pv_can_flip = antisym_v and kv_g >= 2

            if not pu_can_flip and not pv_can_flip:
                g_ach: set[tuple[int, int]] = {(1, 1)}
            elif pu_can_flip and not pv_can_flip:
                g_ach = {(1, 1), (-1, 1)}
            elif not pu_can_flip and pv_can_flip:
                g_ach = {(1, 1), (1, -1)}
            elif u_g == v_g:
                # Both can flip but full overlap in g: parities locked together.
                g_ach = {(1, 1), (-1, -1)}
            else:
                # Both can flip, partial/zero overlap: fully independent.
                g_ach = {(1, 1), (1, -1), (-1, 1), (-1, -1)}

            # Extend achievable by multiplying with this group's achievable set.
            achievable = {
                (pu * qu, pv * qv)
                for (pu, pv) in achievable
                for (qu, qv) in g_ach
            }

        u_flips = (-1, 1) in achievable
        v_flips = (1, -1) in achievable
        both_flip = (-1, -1) in achievable

        if u_flips and v_flips:
            return SignCorrelationType.TYPE_22
        elif u_flips:
            return SignCorrelationType.TYPE_21
        elif v_flips:
            return SignCorrelationType.TYPE_12
        elif both_flip:
            return SignCorrelationType.TYPE_NEG
        else:
            return SignCorrelationType.TYPE_11

    def sign_correlation_type_brute(self, u: Atom, v: Atom) -> SignCorrelationType:
        """
        Brute-force O(∏ n_g!) sign-correlation computation.  Permanent reference.

        A sign "flips freely" for one atom iff there exists an orbit element
        where that atom's sign differs from the canonical representative while
        the other atom's sign is *unchanged*.  This is the empirical definition
        — it is computed from the actual orbit, not assumed from ArgumentSymmetry.

        For SYMMETRIC operations the sign is always +1 and can never flip, so
        those atoms always contribute a "1" to the type label.

        See SignCorrelationType for a full description of each type.
        """
        orb = self.orbit_brute(u, v)

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
        # both_flip: ∃ (a, b) in orbit where both signs differ from the canonical
        both_flip = any(
            a.sign != u.sign and b.sign != v.sign
            for a, b in orb
        )

        if u_flips_freely and v_flips_freely:
            return SignCorrelationType.TYPE_22
        elif u_flips_freely:
            return SignCorrelationType.TYPE_21
        elif v_flips_freely:
            return SignCorrelationType.TYPE_12
        elif both_flip:
            return SignCorrelationType.TYPE_NEG
        else:
            return SignCorrelationType.TYPE_11

