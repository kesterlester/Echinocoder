"""
symcoder.encoders.sort_encoder
===============================
SortEncoder: simplest possible atom-orbit embedding — evaluate every atom in the
orbit, then sort the results.

Accepts all three OrbitSpec forms.  Has a low default priority (0.5) so that
more sophisticated encoders are preferred when available.  Useful as a fallback
and as a concrete reference implementation to validate the registry infrastructure.
"""
from __future__ import annotations

import numpy as np

from symatom.context import Plan

from ..eval import evaluate
from ._base import (
    AtomOrbitEncoder,
    OrbitSpec,
    OrbitSpecForm,
    EncodingCapability,
    EncodingResult,
)

# Priority lower than PolyEncoder so the polynomial embedding is preferred
# when both are capable.
_PRIORITY = 0.5


class SortEncoder(AtomOrbitEncoder):
    """
    Embeds an atom orbit as a sorted vector of real evaluations.

    output_dim equals the orbit size (one real per atom).  The embedding is
    permutation-invariant because sorting is permutation-invariant.

    For ANTISYMMETRIC operations the orbit contains {a_k, -a_k} pairs, so the
    sorted output is antisymmetric about 0 and only half the entries carry
    independent information.  This encoder does not exploit that redundancy;
    PolyEncoder provides a more compressed alternative for those orbits.

    Implementation guide (stubs below raise NotImplementedError)
    ------------------------------------------------------------

    assess()
    ~~~~~~~~
    Determine orbit_size from the spec form:

        return EncodingCapability(
            can_encode=True,
            output_dim=orbit_size,
            method_name="sort",
            priority=_PRIORITY,
            metadata={"orbit_size": orbit_size},
        )

    encode()
    ~~~~~~~~
    Obtain the list of Atoms for the orbit:

    REPRESENTATIVE_ATOM:
        atom = spec.payload
        # Construct the FlavouredOperator as above, then iterate its atoms.
        fo = FlavouredOperator(atom.operation, flavour, plan.context)
        atoms = list(fo.atoms_one_per_sign())

    FLAVOURED_OPERATOR:
        fo = spec.payload
        atoms = list(fo.atoms_one_per_sign())

    EXPLICIT_ORBIT:
        atoms = spec.payload  # already a list[Atom]

    Evaluate each atom:
        from symcoder.eval import evaluate
        # event maps label strings to numpy arrays, e.g. {"p1": array([...]), ...}
        vals = [evaluate(a, event) for a in atoms]

    Sort and return:
        values = np.sort(np.array(vals, dtype=np.float64))
        return EncodingResult(
            values=values,
            metadata={"method": "sort", "orbit_size": len(atoms)},
        )
    """

    def assess(self, spec: OrbitSpec, plan: Plan) -> EncodingCapability:
        """
        A sort encoder can simply sort a set PROVIDED that the set itself is
        invariant under the group action contained within the plan.
        For example:

           * if the plan's group G permutes {a,b,c} then it's     fine to sort encode S={a.b, a.c, b.c} since G.S = {S}.
           * if the plan's group G permutes {a,b,c} then it's NOT fine to sort encode S={a.b, a.c} since G.S != {S}.
           * if the plan's group G permutes {a,b} then it's     fine to sort encode S={a.b} since G.S = {S}.
           * if the plan's group G permutes {a,b} then it's NOT fine to sort encode S={eps2(a.b)} since G.S != {S}.

        So, most naive implementation takes every element of the group, uses it to 
        modify a list of atoms, and then see if that list is invariant.
        """

        if False and spec.form == OrbitSpecForm.EXPLICIT_ORBIT:
            raise NotImplementedError("We should also check that the orbit really is an orbit, but checking that (1) it is a set, and (2) its set of elements does not change when any group elelment acts on it.")
            the_orbit = spec.payload
            assert len(the_orbit)>=1 # TODO: can remove in future.
            assert len(the_orbit) == plan.context.the_group.orbit_size(the_orbit[0]) # TODO: can remove in future.
            return EncodingCapability( # TODO: remove large amount of code repetition (here and ~20 lines further on).
                can_encode=True,
                output_dim=len(the_orbit),
                method_name="JustEvalAllAtomsInOrbitAndSort", # TODO: name should not be hardcoded in multiple places
                priority=_PRIORITY,
                metadata={"the_orbit": the_orbit},
            )
        # Set canonical representative u:
        if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
            fo = spec.payload
            u = fo.canonical_representative()
        elif spec.form == OrbitSpecForm.REPRESENTATIVE_ATOM:
            u = spec.payload
        else:
            u = None
        if u is not None:
            the_orbit = plan.context.the_group.orbit(u) # What we need to eval and sort encode.
            the_orbit_size = len(the_orbit)
            assert the_orbit_size == plan.context.the_group.orbit_size(u) # Just a check. TODO: remove if never seems to fail.

            return EncodingCapability(
                can_encode=True,
                output_dim=the_orbit_size,
                method_name="JustEvalAllAtomsInOrbitAndSort", # TODO: name should not be hardcoded in multiple places
                priority=_PRIORITY,
                metadata={"the_orbit": the_orbit},
            )
        print(f"sort_encoder did not recognise type of OrbitSpec supplied {spec.form} with payload {spec.payload}.")
        return EncodingCapability(
                can_encode=False,
            )

    def encode(self, capability: EncodingCapability, event: dict, plan: Plan) -> EncodingResult:
        atoms = capability.metadata["the_orbit"]
        print(f"Sort encoder is evaluating {atoms} on Event {event}")
        vals = [evaluate(a, event) for a in atoms]
        print(f"COW COW COW actual size was {len(vals)}.")
        assert len(atoms) == len(vals) # TODO: Can probably delete
        values = np.sort(np.array(vals, dtype=np.float64))
        assert len(atoms) == len(values) # TODO: Can probably delete
        return EncodingResult(
            values=values,
            metadata={"method": "sort", "orbit_size": len(atoms)},
        )
