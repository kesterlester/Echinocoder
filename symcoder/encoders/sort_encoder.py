"""
symcoder.encoders.sort_encoder
===============================
Two sort-based atom-orbit encoders:

SortEncoder / SortEncoderFactory
    Simplest possible embedding: evaluate every atom in the orbit, sort the
    results.  Output length == orbit size.  Handles all orbit types; used as
    the universal fallback.  Priority 0.5.

HalfSortEncoder / HalfSortEncoderFactory
    Compressed embedding for orbits that can be exactly partitioned into
    (atom, -atom) pairs.  For such an orbit {+a,-a, +b,-b, ...} it suffices
    to store sort([|eval(a)|, |eval(b)|, ...]), halving the output length.
    Eligibility is determined solely by physically attempting to partition the
    orbit into negative pairs — NOT by inspecting ArgumentSymmetry flags.
    Output length == orbit_size / 2.  Priority 1.0 (preferred over SortEncoder
    when applicable).
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
from symatom.atoms import are_negatives
from symatom.context import Plan

from ..eval import evaluate
from ..decoded_types import AnnotatedMultisetOfReals
from ._base import (
    AtomOrbitEncoder,
    AtomOrbitEncoderFactory,
    OrbitSpec,
    OrbitSpecForm,
    EncodingResult,
)

_SORT_PRIORITY      = 0.5
_HALF_SORT_PRIORITY = 1.0


# ---------------------------------------------------------------------------
# Shared helper: resolve any OrbitSpec form to a concrete list of atoms
# ---------------------------------------------------------------------------

def _orbit_from_spec(spec: OrbitSpec, plan: Plan) -> list | None:
    """
    Return the orbit atom list for any supported spec form, or None if the
    spec form is unrecognised.

    Sort-based encoders (SortEncoder, HalfSortEncoder) are only valid when the
    atom list is a genuine G-orbit — i.e. closed under the group action in the
    plan — because sorting only produces a group-invariant embedding if applying
    any group element merely permutes the list (leaving the sorted result unchanged).

    Examples of why this matters:

      * G permutes {a,b,c}: sort-encoding {a·b, a·c, b·c} is fine — G maps this
        set to itself.
      * G permutes {a,b,c}: sort-encoding {a·b, a·c} is NOT fine — G.{a·b,a·c}
        contains b·c, which is outside the list.
      * G permutes {a,b}: sort-encoding {eps2(a,b)} is NOT fine — G maps it to
        {eps2(b,a)} = {-eps2(a,b)}, a different set.

    For FLAVOURED_OPERATOR and REPRESENTATIVE_ATOM specs the orbit is computed
    by the group directly, so validity is guaranteed.  For EXPLICIT_ORBIT the
    caller is trusted to have supplied a genuine orbit; no runtime check is made
    (implementing the check would require verifying G-closure over the full atom
    list, which is non-trivial and was left as a TODO).
    """
    if spec.form == OrbitSpecForm.EXPLICIT_ORBIT:
        return spec.payload
    if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
        u = spec.payload.canonical_representative()
        return plan.context.the_group.orbit(u)
    if spec.form == OrbitSpecForm.REPRESENTATIVE_ATOM:
        return plan.context.the_group.orbit(spec.payload)
    return None


# ---------------------------------------------------------------------------
# Shared helper: attempt to partition an orbit into (atom, -atom) pairs
# ---------------------------------------------------------------------------

def _try_neg_pair_partition(orbit: list) -> list | None:
    """
    Try to partition orbit into (atom, -atom) pairs.

    Returns a list of one representative per pair (the sign=+1 atom of each
    pair), or None if the orbit cannot be so partitioned.

    The check is purely structural — atoms are grouped by (operation, labels)
    and each group must contain exactly two atoms that satisfy are_negatives().
    No inference from ArgumentSymmetry flags is made.
    """
    if len(orbit) % 2 != 0:
        return None

    groups: dict = defaultdict(list)
    for atom in orbit:
        groups[(atom.operation, atom.labels)].append(atom)

    representatives = []
    for _key, pair in groups.items():
        if len(pair) != 2:
            return None
        a, b = pair
        if not are_negatives(a, b):
            return None
        representatives.append(a if a.sign == +1 else b)

    return representatives


# ---------------------------------------------------------------------------
# SortEncoder
# ---------------------------------------------------------------------------

class SortEncoder(AtomOrbitEncoder):
    """
    Embeds an atom orbit as a sorted vector of real evaluations.

    Constructed by SortEncoderFactory.assess() with the orbit pre-computed.
    encode(event) requires only the event data.

    For orbits whose atoms come in (atom, -atom) pairs, HalfSortEncoderFactory
    provides a more compressed alternative and will be preferred when both
    factories are registered.
    """

    def __init__(self, orbit: list) -> None:
        self._orbit = list(orbit)

    @property
    def output_dim(self) -> int:
        return len(self._orbit)

    @property
    def priority(self) -> float:
        return _SORT_PRIORITY

    @property
    def method_name(self) -> str:
        return "sort"

    def encode(self, event: dict) -> EncodingResult:
        vals = [evaluate(a, event) for a in self._orbit]
        values = np.sort(np.array(vals, dtype=np.float64))
        return EncodingResult(
            values=values,
            metadata={"method": "sort", "orbit_size": len(self._orbit)},
        )

    def decode(self, values: np.ndarray) -> AnnotatedMultisetOfReals:
        """Decode a sorted orbit slice back to an annotated multiset of reals.

        The stored values ARE the orbit evaluations (just sorted); no numerical
        inversion is needed.  The full orbit is returned as-is.
        """
        return AnnotatedMultisetOfReals(
            values=list(values),
            atoms=list(self._orbit),
        )


# ---------------------------------------------------------------------------
# SortEncoderFactory
# ---------------------------------------------------------------------------

class SortEncoderFactory(AtomOrbitEncoderFactory):
    """
    Factory that produces SortEncoder instances for any well-formed OrbitSpec.

    Handles all three OrbitSpec forms.  Returns [] only for unrecognised forms.
    Use as the universal fallback in an OrbitEncoderFactory.
    """

    def assess(self, spec: OrbitSpec, plan: Plan) -> list[AtomOrbitEncoder]:
        orbit = _orbit_from_spec(spec, plan)
        if orbit is None:
            return []
        return [SortEncoder(orbit)]


# ---------------------------------------------------------------------------
# HalfSortEncoder
# ---------------------------------------------------------------------------

class HalfSortEncoder(AtomOrbitEncoder):
    """
    Compressed orbit encoder for orbits that partition into (atom, -atom) pairs.

    Stores one representative atom per pair (the sign=+1 member).  encode()
    evaluates each representative, takes the absolute value, and sorts the
    results.  Output length is half the orbit size.

    Constructed by HalfSortEncoderFactory.assess(), which verifies the pairing
    condition physically before returning this encoder.
    """

    def __init__(self, orbit: list, representatives: list) -> None:
        self._orbit = list(orbit) # Full orbit.
        self._representatives = list(representatives) # Half of the orbit, whose elements are evaluated and modded.

    @property
    def output_dim(self) -> int:
        return len(self._representatives)

    @property
    def priority(self) -> float:
        return _HALF_SORT_PRIORITY

    @property
    def method_name(self) -> str:
        return "half_sort"

    @property
    def notional_output_dim(self) -> int:
        return 2 * len(self._representatives)

    def encode(self, event: dict) -> EncodingResult:
        vals = [abs(evaluate(a, event)) for a in self._representatives]
        values = np.sort(np.array(vals, dtype=np.float64))
        return EncodingResult(
            values=values,
            metadata={"method": "half_sort", "orbit_size": 2 * len(self._representatives)},
        )

    def decode(self, values: np.ndarray) -> AnnotatedMultisetOfReals:
        """Decode a half-sorted orbit slice back to the full annotated multiset.

        The stored values are |eval(+atom)| for each representative.  The full
        orbit also contains the negatives, so we expand: [v1,…,vn] → [v1,…,vn,
        -v1,…,-vn].  Atom list similarly includes both signs.
        """
        from symatom.atoms import Atom
        pos_vals = list(values)
        neg_vals = [-v for v in pos_vals]
        neg_atoms = [Atom(r.operation, r.labels, sign=-1) for r in self._representatives]
        return AnnotatedMultisetOfReals(
            values=pos_vals + neg_vals,
            atoms=list(self._representatives) + neg_atoms,
        )


# ---------------------------------------------------------------------------
# HalfSortEncoderFactory
# ---------------------------------------------------------------------------

class HalfSortEncoderFactory(AtomOrbitEncoderFactory):
    """
    Factory that produces HalfSortEncoder instances for orbits that can be
    partitioned into (atom, -atom) pairs.

    assess() physically attempts the partition using _try_neg_pair_partition().
    Returns [] — deferring to a fallback factory — if the orbit cannot be so
    partitioned or if the spec form is unrecognised.
    """

    def assess(self, spec: OrbitSpec, plan: Plan) -> list[AtomOrbitEncoder]:
        orbit = _orbit_from_spec(spec, plan)
        if orbit is None:
            return []
        reps = _try_neg_pair_partition(orbit)
        if reps is None:
            return []
        return [HalfSortEncoder(orbit, reps)]
