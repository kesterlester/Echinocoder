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

    REPRESENTATIVE_ATOM:
        atom = spec.payload  # type: Atom
        # Build a FlavouredOperator from the representative atom to count the orbit.
        # The flavour records how many of the atom's labels come from each group.
        from symatom.rep import FlavouredOperator, Flavour
        flavour = Flavour(tuple(
            sum(1 for lbl in atom.labels if lbl in set(g.labels))
            for g in plan.context.types
        ))
        fo = FlavouredOperator(atom.operation, flavour, plan.context)
        orbit_size = fo.count()

    FLAVOURED_OPERATOR:
        fo = spec.payload  # type: FlavouredOperator
        orbit_size = fo.count()

    EXPLICIT_ORBIT:
        orbit_size = len(spec.payload)

    Then return:
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
        atoms = list(fo.atoms())

    FLAVOURED_OPERATOR:
        fo = spec.payload
        atoms = list(fo.atoms())

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
        # A sort encoder can simply sort a set PROVIDED that the set itself is invariant under the group action contained within the plan.
        # e.g.:
        #    * if the plan's group G permutes {a,b,c} then it's     fine to sort encode S={a.b, a.c, b.c} since G.S = {S}.
        #    * if the plan's group G permutes {a,b,c} then it's NOT fine to sort encode S={a.b, a.c} since G.S != {S}.
        #    * if the plan's group G permutes {a,b} then it's     fine to sort encode S={a.b} since G.S = {S}.
        #    * if the plan's group G permutes {a,b} then it's NOT fine to sort encode S={eps2(a.b)} since G.S != {S}.
        # So, most naive implementation takes every element of the group, uses it to modify a list of atoms, and then see if that list is invariant:

        if spec.form == OrbitSpecForm.EXPLICIT_ORBIT:
            orbit = spec.payload
            return EncodingCapability(
                can_encode=True,
                output_dim=len(orbit),
                method_name="sortMOOMOOMOO",
                priority=_PRIORITY,
                metadata={"orbit_size_moo": 100},
            )
        if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
            fo = spec.payload
            return EncodingCapability(
                can_encode=True,
                output_dim=fo.count(),
                method_name="sortMOOMOOMOO",
                priority=_PRIORITY,
                metadata={"orbit_size_moo": 100},
            )
        print(f"sort_encoder did not recognise type of OrbitSpec supplied {spec.form}")
        return EncodingCapability(
                can_encode=False,
            )

    def encode(self, spec: OrbitSpec, event: dict, plan: Plan) -> EncodingResult:
        if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
            fo = spec.payload
            atoms = list(fo.atoms())
            assert len(atoms) == fo.count()
        else:
            print(f"sort_encoder did not recognise type of OrbitSpec supplied {spec.form}")
            raise NotImplementedError(
                "SortEncoder.encode() is a stub — implement atom evaluation and sorting "
                "(see class docstring)"
            )
        vals = [evaluate(a, event) for a in atoms]
        values = np.sort(np.array(vals, dtype=np.float64))
        return EncodingResult(
            values=values,
            metadata={"method": "sort", "orbit_size": len(atoms)},
        )
