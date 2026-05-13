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
        # TODO: implement orbit-size logic (see class docstring for all three forms)
        raise NotImplementedError(
            "SortEncoder.assess() is a stub — implement the orbit-size calculation "
            "for all three OrbitSpecForm variants (see class docstring)"
        )

    def encode(self, spec: OrbitSpec, event: dict, plan: Plan) -> EncodingResult:
        # TODO: implement atom enumeration + evaluate + sort (see class docstring)
        raise NotImplementedError(
            "SortEncoder.encode() is a stub — implement atom evaluation and sorting "
            "(see class docstring)"
        )
