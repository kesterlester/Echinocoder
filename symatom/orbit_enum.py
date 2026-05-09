from __future__ import annotations
from abc import ABC, abstractmethod


class OrbitEnumerator(ABC):
    """
    Strategy interface for enumerating the G-orbit of an atom-pair given its
    PairFlavour.  Two implementations are provided:

    - BruteForceOrbitEnumerator  (correct by inspection; always available)
    - DirectOrbitEnumerator      (fast combinatorial construction; not yet implemented)

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
    Cross-product + filter implementation.  Delegates to pf.orbit_elements(),
    which enumerates FlavouredOperator.atoms() for each side and keeps only
    pairs whose pair_flavour_of matches pf.

    Correct by inspection: every atom in FlavouredOperator.atoms() is a
    distinct concrete atom with the right flavour; pair_flavour_of is the
    canonical orbit-type test.  Use this for testing and as the reference
    against which DirectOrbitEnumerator is validated.
    """

    def orbit_elements(self, pf, context) -> list:
        return pf.orbit_elements(context)


class DirectOrbitEnumerator(OrbitEnumerator):
    """
    Combinatorial construction from the PairFlavour overlap structure.

    For large n (target: n ~ 30 per species) the cross-product approach in
    BruteForceOrbitEnumerator is O(atoms_u × atoms_v) ≫ O(orbit_size).
    This implementation will instead iterate directly over the combinatorial
    structure encoded in pf.overlap — choosing shared labels, then u-only
    labels, then v-only labels — producing exactly orbit_size pairs without
    any filtering step.

    NOT YET IMPLEMENTED.  All tests that exercise this path are marked
    pytest.mark.skip until the implementation is complete.  Once implemented,
    the parametrised BruteForce vs Direct tests in test_orbit_enum.py serve
    as the correctness cross-check.
    """

    def orbit_elements(self, pf, context) -> list:
        raise NotImplementedError(
            "DirectOrbitEnumerator is not yet implemented.  "
            "Use BruteForceOrbitEnumerator for now."
        )
