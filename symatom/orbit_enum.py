from __future__ import annotations
from abc import ABC, abstractmethod
from symatom.atoms import Atom


class OrbitEnumerator(ABC):
    """
    Strategy interface for enumerating the G-orbit of an atom-pair given its
    PairFlavour.  Two implementations are provided:

    - BruteForceOrbitEnumerator  (correct by inspection; always available)
    - DirectOrbitEnumerator      (fast combinatorial construction)

    Pass the desired enumerator to Plan; symcoder calls
    plan.orbit_enumerator.orbit_elements(pf, context) without caring which
    implementation is active.
    """

    @abstractmethod
    def orbit_elements(self, pf, context) -> list:
        """
        Return all (Atom, Atom) pairs in the G-orbit of pf's canonical
        representative, given a context.
        """


class BruteForceOrbitEnumerator(OrbitEnumerator):
    """
    Delegates to pf.orbit_elements(context).

    pf.orbit_elements applies every element of G = S_{n_1} × ... × S_{n_m}
    simultaneously to the canonical pair representative, collecting unique
    images.  This is correct by direct inspection of the group action and
    serves as the permanent reference implementation against which
    DirectOrbitEnumerator is validated.

    Both enumerators must return the same set for every PairFlavour.
    """

    def orbit_elements(self, pf, context) -> list:
        return pf.orbit_elements(context)


class DirectOrbitEnumerator(OrbitEnumerator):
    """
    Returns the true G-orbit of a PairFlavour's canonical representative.

    This implementation delegates to TheGroup(context).orbit(u_canon, v_canon),
    where (u_canon, v_canon) is constructed from pf in the same way as
    pf.orbit_elements.  The result is identical to BruteForceOrbitEnumerator
    for all PairFlavours — the parametrised cross-comparison tests in
    test_orbit_enum.py verify this.

    Relationship to BruteForceOrbitEnumerator
    -----------------------------------------
    BruteForceOrbitEnumerator is the permanent reference: it delegates to
    pf.orbit_elements(context) and must never be changed.

    DirectOrbitEnumerator is the production path: it currently uses
    TheGroup.orbit (also brute-force O(∏ n_g!)), but is designed to be
    progressively replaced by algebraic shortcut methods in Step 3 of the
    architecture redesign — without touching BruteForce.

    Previous design (OrbitUnion, now removed)
    ------------------------------------------
    The earlier implementation iterated over label assignments from the
    PairFlavour combinatorial structure and emitted all sign combinations
    per assignment (±u × ±v).  This produced the *union* of all G-orbits of
    a given PairFlavour type — an OrbitUnion, not a true G-orbit.  For
    ANTISYMMETRIC pairs with non-zero overlap the OrbitUnion was larger than
    the true orbit, causing incorrect orbit sizes and breaking the invariant
    that Direct == BruteForce.  The OrbitUnion approach has been replaced by
    TheGroup-based enumeration in Step 2.
    """

    def orbit_elements(self, pf, context) -> list:
        from .group import TheGroup

        # Build the canonical representative in the same way pf.orbit_elements does
        group_sizes = tuple(g.size for g in context.groups)
        if pf.count(group_sizes) == 0:
            return []

        u_labels, v_labels = [], []
        for g, ku, kv, s in zip(
            context.groups,
            pf.flavour_u.counts,
            pf.flavour_v.counts,
            pf.overlap,
        ):
            u_labels.extend(g.labels[:ku])
            v_labels.extend(g.labels[:s])
            v_labels.extend(g.labels[ku:ku + kv - s])

        u_canon = Atom(pf.op_u, tuple(u_labels), sign=+1)
        v_canon = Atom(pf.op_v, tuple(v_labels), sign=+1)

        return TheGroup(context).orbit(u_canon, v_canon)
