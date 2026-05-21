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
    Reference enumerator — delegates to pf.orbit_elements(context).

    pf.orbit_elements() in turn delegates to context.the_group.orbit_brute_pair(),
    which is the single authoritative brute-force orbit implementation inside
    TheGroup.  All orbit logic lives in TheGroup; this class is a thin wrapper
    that makes the reference path available through the OrbitEnumerator interface.

    Role: permanent reference implementation, validated against DirectOrbitEnumerator
    in test_orbit_enum.py.  Both enumerators must return the same set for every
    PairFlavour.
    """

    def orbit_elements(self, pf, context) -> list:
        return pf.orbit_elements(context)


class DirectOrbitEnumerator(OrbitEnumerator):
    """
    Production enumerator — delegates to context.the_group.orbit_pair().

    Calls pf.canonical_pair(context) to obtain the canonical representative
    atom-pair (the single definition of "canonical" lives on PairFlavour), then
    asks TheGroup to enumerate the orbit.

    Relationship to BruteForceOrbitEnumerator
    -----------------------------------------
    Both enumerators delegate orbit computation to TheGroup — BruteForce via
    TheGroup.orbit_brute_pair() (always brute-force), Direct via
    TheGroup.orbit_pair() (currently also brute-force, but designed to be
    replaced by an algebraic O(orbit_size) implementation in Step 4 without
    touching BruteForce).  The parametrised cross-comparison tests in
    test_orbit_enum.py verify that both return the same set for every
    PairFlavour.

    Previous design (OrbitUnion, now removed)
    ------------------------------------------
    The earlier implementation iterated over label assignments and emitted all
    ±u × ±v sign combinations, producing the *union* of all G-orbits of a
    given PairFlavour type rather than a single true G-orbit.  For ANTISYMMETRIC
    pairs with non-zero overlap the OrbitUnion was larger than the true orbit.
    Replaced by TheGroup-based enumeration in Step 2.
    """

    def orbit_elements(self, pf, context) -> list:
        type_sizes = tuple(g.size for g in context.types)
        if pf.count(type_sizes) == 0:
            return []
        u_canon, v_canon = pf.canonical_pair(context)
        return context.the_group.orbit_pair(u_canon, v_canon)
