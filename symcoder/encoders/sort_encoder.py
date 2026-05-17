"""
symcoder.encoders.sort_encoder
===============================
SortEncoderFactory / SortEncoder: simplest possible atom-orbit embedding —
evaluate every atom in the orbit, then sort the results.

SortEncoderFactory handles all three OrbitSpec forms.  It has a low default
priority (0.5) so that more sophisticated encoders are preferred when available.
Useful as a fallback and as a concrete reference implementation.
"""
from __future__ import annotations

import numpy as np

from symatom.context import Plan

from ..eval import evaluate
from ._base import (
    AtomOrbitEncoder,
    AtomOrbitEncoderFactory,
    OrbitSpec,
    OrbitSpecForm,
    EncodingResult,
)

_PRIORITY    = 0.5
_METHOD_NAME = "sort"


# ---------------------------------------------------------------------------
# SortEncoder — the ready-to-use encoder produced by SortEncoderFactory
# ---------------------------------------------------------------------------

class SortEncoder(AtomOrbitEncoder):
    """
    Embeds an atom orbit as a sorted vector of real evaluations.

    Constructed by SortEncoderFactory.assess() with the orbit pre-computed.
    encode(event) requires only the event data.

    For ANTISYMMETRIC operations the orbit contains {a_k, -a_k} pairs, so the
    sorted output is antisymmetric about 0 and only half the entries carry
    independent information.  This encoder does not exploit that redundancy;
    PolyEncoderFactory provides a more compressed alternative for those orbits.
    """

    def __init__(self, orbit: list) -> None:
        self._orbit = list(orbit)

    @property
    def output_dim(self) -> int:
        return len(self._orbit)

    @property
    def priority(self) -> float:
        return _PRIORITY

    @property
    def method_name(self) -> str:
        return _METHOD_NAME

    def encode(self, event: dict) -> EncodingResult:
        print(f"Sort encoder is evaluating {self._orbit} on Event {event}")
        vals = [evaluate(a, event) for a in self._orbit]
        values = np.sort(np.array(vals, dtype=np.float64))
        return EncodingResult(
            values=values,
            metadata={"method": _METHOD_NAME, "orbit_size": len(self._orbit)},
        )


# ---------------------------------------------------------------------------
# SortEncoderFactory — creates SortEncoder instances
# ---------------------------------------------------------------------------

class SortEncoderFactory(AtomOrbitEncoderFactory):
    """
    Factory that produces SortEncoder instances for any well-formed OrbitSpec.

    assess() resolves the orbit from the spec (using TheGroup.orbit() for
    REPRESENTATIVE_ATOM and FLAVOURED_OPERATOR forms, or taking the payload
    directly for EXPLICIT_ORBIT), then returns a single SortEncoder wrapping
    that orbit.

    Returns [] for unrecognised spec forms.
    """

    def assess(self, spec: OrbitSpec, plan: Plan) -> list[AtomOrbitEncoder]:
        """
        A sort encoder can simply sort a set PROVIDED that the set itself is
        invariant under the group action contained within the plan.
        For example:

           * if the plan's group G permutes {a,b,c} then it's     fine to sort encode S={a.b, a.c, b.c} since G.S = {S}.
           * if the plan's group G permutes {a,b,c} then it's NOT fine to sort encode S={a.b, a.c} since G.S != {S}.
           * if the plan's group G permutes {a,b} then it's     fine to sort encode S={a.b} since G.S = {S}.
           * if the plan's group G permutes {a,b} then it's NOT fine to sort encode S={eps2(a.b)} since G.S != {S}.

        So, the most naive implementation takes every element of the group, uses it to
        modify a list of atoms, and then sees if that list is invariant.
        """

        # TODO: For EXPLICIT_ORBIT we should verify the payload really is a G-orbit:
        #   (1) it is a set, and (2) its set of elements does not change when any
        #   group element acts on it.  Disabled for now (if False) while the check
        #   is not yet implemented.
        if False and spec.form == OrbitSpecForm.EXPLICIT_ORBIT:  # noqa: SIM210
            raise NotImplementedError(
                "We should also check that the orbit really is an orbit, by "
                "checking that (1) it is a set, and (2) its set of elements does "
                "not change when any group element acts on it."
            )
            orbit = spec.payload
            assert len(orbit) >= 1                                          # TODO: can remove in future
            assert len(orbit) == plan.context.the_group.orbit_size(orbit[0])  # TODO: can remove in future
            return [SortEncoder(orbit)]  # TODO: remove large code repetition (here and ~15 lines further on)

        if spec.form == OrbitSpecForm.EXPLICIT_ORBIT:
            orbit = spec.payload
            return [SortEncoder(orbit)]

        if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
            u = spec.payload.canonical_representative()
        elif spec.form == OrbitSpecForm.REPRESENTATIVE_ATOM:
            u = spec.payload
        else:
            return []

        orbit = plan.context.the_group.orbit(u)
        assert len(orbit) == plan.context.the_group.orbit_size(u)  # TODO: remove if never seems to fail
        return [SortEncoder(orbit)]
