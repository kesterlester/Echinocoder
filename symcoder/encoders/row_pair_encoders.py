"""
symcoder.encoders.row_pair_encoders
=====================================
Concrete row-pair encoders and their factories for Phase 2 encoding.

Each factory's assess() physically partitions the pair orbit to decide whether
it can handle the pair.  No factory calls _sign_correlation_type_from_pf; each
is self-contained.

Factories and what they handle
--------------------------------
SelfPairEncoderFactory    self-pairs (any type)    → NullPairEncoder    output_dim=0
Type11PairEncoderFactory  any non-self pair        → Type11PairEncoder  output_dim=2n  (universal fallback)
Type12PairEncoderFactory  TYPE_12 pairs            → Type12PairEncoder  output_dim=2n
Type21PairEncoderFactory  TYPE_21 pairs            → Type21PairEncoder  output_dim=2n
Type22PairEncoderFactory  TYPE_22 pairs            → Type22PairEncoder  output_dim=2n
NegPairEncoderFactory     TYPE_NEG pairs           → NegPairEncoder     output_dim=2n

Here n is the number of base atom-pairs (one per sign-equivalence class in the
orbit).  All five non-null encoders produce the same output_dim=2n.  The
overlap-block manager picks among offered encoders by minimum output_dim; for
ties it takes the first in the list (favouring the specialist registered first).

Type11 as universal fallback
-----------------------------
Type11PairEncoderFactory accepts any non-self-pair without inspection.
Specialists (Type12/21/22/Neg) each have a smaller achievable set and therefore
fewer orbit elements — but because output_dim=2n is the SAME for all types, the
min(output_dim) selection prefers specialists (smaller n means smaller 2n).
For a pure TYPE_11 orbit only specialists fail their assessments, leaving
Type11 as the sole offer.

Physical partition
------------------
For TYPE_21/12/22/NEG factories, assess() groups orbit elements by
(base_u, base_v) and checks the group structure.  The representative stored in
each encoder is the unique element with u.sign == +1 per group (or u.sign==+1,
v.sign==+1 for TYPE_22).  encode() evaluates those representatives and
reconstructs the full invariant polynomial without consulting the orbit again.

Self-pair null encoding
-----------------------
Phase 1 stores every individual row orbit.  A self-pair (u=v) z_k is fully
determined by Phase 1, so SelfPairEncoderFactory returns NullPairEncoder
(output_dim=0) for self-pairs and [] for all other pairs.
"""
from __future__ import annotations

import dataclasses
from abc import abstractmethod
from collections import defaultdict

import numpy as np

from symatom.atoms import are_negatives
from symcoder.sign_correlation import _complex_to_reals
from symcoder.pairs import _is_self_pair
from symcoder.describe import SegmentInfo, _assoc_example
from symcoder.eval import evaluate
from ._embedder import zip_embed
from .pair_base import PairOrbitEncoder, PairOrbitEncoderFactory, PairOrbitSpec, EncodingResult

_PRIORITY = 0.5


# ---------------------------------------------------------------------------
# Orbit partition helpers
# ---------------------------------------------------------------------------

def _pair_orbit_atoms(pf, plan) -> list:
    """Return all (u_atom, v_atom) pairs in the G-orbit of pf."""
    return list(plan.orbit_enumerator.orbit_elements(pf, plan.context))


def _base_atom(atom):
    """Return the positive-sign representative of atom (same op, same labels, sign=+1)."""
    from symatom.atoms import Atom
    return Atom(atom.operation, atom.labels, sign=+1)


def _group_by_base(pairs: list) -> dict:
    """
    Group (u_atom, v_atom) pairs by (base_u, base_v).

    base_x = Atom(x.operation, x.labels, sign=+1), i.e. the canonical
    sign=+1 form of x.  The returned dict maps each (base_u, base_v) key to
    the list of orbit elements that share that base.
    """
    groups: dict = defaultdict(list)
    for u, v in pairs:
        key = (_base_atom(u), _base_atom(v))
        groups[key].append((u, v))
    return dict(groups)


def _try_type21_partition(pairs: list):
    """
    Return positive-u representatives if the orbit is TYPE_21; else None.

    TYPE_21 criterion: every (base_u, base_v) group has exactly 2 elements
    where are_negatives(u1, u2) is True and v1 == v2 (v is unaffected).
    Representative: the element with u.sign == +1.
    """
    groups = _group_by_base(pairs)
    reps = []
    for group in groups.values():
        if len(group) != 2:
            return None
        (u1, v1), (u2, v2) = group
        if not (are_negatives(u1, u2) and v1 == v2):
            return None
        rep = (u1, v1) if u1.sign == +1 else (u2, v2)
        reps.append(rep)
    return reps


def _try_type12_partition(pairs: list):
    """
    Return positive-v representatives if the orbit is TYPE_12; else None.

    TYPE_12 is TYPE_21 with u/v roles swapped: every group has 2 elements
    where u1 == u2 and are_negatives(v1, v2).
    Implemented as a thin wrapper around _try_type21_partition.
    Representative: the element with v.sign == +1.
    """
    swapped = [(v, u) for u, v in pairs]
    swapped_reps = _try_type21_partition(swapped)
    if swapped_reps is None:
        return None
    return [(u, v) for v, u in swapped_reps]


def _try_type22_partition(pairs: list):
    """
    Return (+u, +v) representatives if the orbit is TYPE_22; else None.

    TYPE_22 criterion: every (base_u, base_v) group has exactly 4 elements
    covering all four sign combinations {(+,+),(+,-),(-,+),(-,-)}.
    Representative: the element with u.sign==+1 and v.sign==+1.
    """
    groups = _group_by_base(pairs)
    reps = []
    for group in groups.values():
        if len(group) != 4:
            return None
        signs = frozenset((u.sign, v.sign) for u, v in group)
        if signs != frozenset({(1, 1), (1, -1), (-1, 1), (-1, -1)}):
            return None
        rep = next((u, v) for u, v in group if u.sign == 1 and v.sign == 1)
        reps.append(rep)
    return reps


def _try_neg_partition(pairs: list):
    """
    Return sign=+1-u representatives if the orbit is TYPE_NEG; else None.

    TYPE_NEG criterion: every (base_u, base_v) group has exactly 2 elements
    where are_negatives(u1, u2) AND are_negatives(v1, v2).

    This covers both NEG forms:
      Form 1 — {(+u,+v),(-u,-v)}: representative is (+u,+v)
      Form 2 — {(+u,-v),(-u,+v)}: representative is (+u,-v)

    In both forms the orbit of z = eval(u)+i*eval(v) is {z, -z}, so z² is
    invariant regardless of which form.  Picking the rep with u.sign==+1
    works for both.
    """
    groups = _group_by_base(pairs)
    reps = []
    for group in groups.values():
        if len(group) != 2:
            return None
        (u1, v1), (u2, v2) = group
        if not (are_negatives(u1, u2) and are_negatives(v1, v2)):
            return None
        rep = (u1, v1) if u1.sign == +1 else (u2, v2)
        reps.append(rep)
    return reps


# ---------------------------------------------------------------------------
# _BasePairEncoder — shared scaffolding for all row-pair encoders
# ---------------------------------------------------------------------------

class _BasePairEncoder(PairOrbitEncoder):
    """
    Base class for concrete row-pair encoders.

    Constructed with a list of (u_atom, v_atom) representatives — the subset
    of orbit elements that suffice for invariant construction.  encode() calls
    evaluate() on these directly; no orbit re-enumeration at encode time.

    _NOTIONAL_FACTOR is the compression ratio |achievable_set|:
    TYPE_11=1, TYPE_12=TYPE_21=NEG=2, TYPE_22=4.  Used in describe() to report
    the uncompressed notional_length.
    """
    _METHOD_NAME:     str = "?"
    _NOTIONAL_FACTOR: int = 1

    def __init__(self, reps: list, pf, plan) -> None:
        self._reps       = reps
        self._pf         = pf
        self._plan       = plan
        self._n          = len(reps)

    @property
    def output_dim(self) -> int:
        return 2 * self._n

    @property
    def priority(self) -> float:
        return _PRIORITY

    @property
    def method_name(self) -> str:
        return self._METHOD_NAME

    def describe(self) -> list[SegmentInfo]:
        types = self._plan.context.types
        return [SegmentInfo(
            kind            = "ASSOC",
            start           = 0,   # placeholder; overlap-block encoder adjusts
            length          = self.output_dim,
            op_u            = self._pf.op_u.name,
            flavour_u       = tuple(self._pf.flavour_u.counts),
            op_v            = self._pf.op_v.name,
            flavour_v       = tuple(self._pf.flavour_v.counts),
            overlap         = tuple(self._pf.overlap),
            symmetry_class  = self._METHOD_NAME,
            notional_length = self._NOTIONAL_FACTOR * self.output_dim,
            method_name     = self.method_name,
            example         = _assoc_example(
                self._pf.op_u.name, self._pf.flavour_u.counts,
                self._pf.op_v.name, self._pf.flavour_v.counts,
                self._pf.overlap, types,
            ),
        )]

    @abstractmethod
    def encode(self, event: dict) -> EncodingResult: ...


# ---------------------------------------------------------------------------
# NullPairEncoder / SelfPairEncoderFactory
# ---------------------------------------------------------------------------

class NullPairEncoder(_BasePairEncoder):
    """
    Encodes a self-pair with zero outputs.

    A self-pair (u=v) adds no information beyond what Phase 1 already stores.
    This encoder signals that fact by contributing no values.
    """
    _METHOD_NAME = "null_self"

    def __init__(self, pf, plan) -> None:
        type_sizes = tuple(g.size for g in plan.context.types)
        self._pf   = pf
        self._plan = plan
        self._reps = []
        self._n    = pf.count(type_sizes)

    @property
    def output_dim(self) -> int:
        return 0

    def encode(self, event: dict) -> EncodingResult:
        return EncodingResult(values=np.array([], dtype=np.float64))

    def describe(self) -> list[SegmentInfo]:
        types = self._plan.context.types
        return [SegmentInfo(
            kind            = "NULL_SELF",
            start           = 0,
            length          = 0,
            op_u            = self._pf.op_u.name,
            flavour_u       = tuple(self._pf.flavour_u.counts),
            op_v            = self._pf.op_v.name,
            flavour_v       = tuple(self._pf.flavour_v.counts),
            overlap         = tuple(self._pf.overlap),
            symmetry_class  = self._METHOD_NAME,
            notional_length = 2 * self._n,
            method_name     = self.method_name,
            example         = _assoc_example(
                self._pf.op_u.name, self._pf.flavour_u.counts,
                self._pf.op_v.name, self._pf.flavour_v.counts,
                self._pf.overlap, types,
            ),
        )]


class SelfPairEncoderFactory(PairOrbitEncoderFactory):
    """Returns NullPairEncoder for self-pairs; [] for all other pairs."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return [NullPairEncoder(spec.pf, plan)]
        return []


# ---------------------------------------------------------------------------
# Type11PairEncoder / Type11PairEncoderFactory  (TYPE_11, universal fallback)
# ---------------------------------------------------------------------------

class Type11PairEncoder(_BasePairEncoder):
    """
    TYPE_11 encoding: no sign structure to exploit.

    Stores the full orbit (all elements have sign (+1,+1) for a true TYPE_11
    orbit, or all elements for non-TYPE_11 when used as fallback).
    Embeds all n z-values directly: n complex → 2n reals.
    """
    _METHOD_NAME     = "11"
    _NOTIONAL_FACTOR = 1

    def encode(self, event: dict) -> EncodingResult:
        z = np.array(
            [complex(evaluate(u, event), evaluate(v, event)) for u, v in self._reps],
            dtype=complex,
        )
        coeffs, _, _ = zip_embed(z)
        return EncodingResult(
            values=_complex_to_reals(coeffs),
            metadata={"method": self._METHOD_NAME},
        )


class Type11PairEncoderFactory(PairOrbitEncoderFactory):
    """
    Returns Type11PairEncoder for any non-self-pair; [] for self-pairs.

    Universal fallback: no orbit structure is assumed.  All orbit elements are
    stored as representatives.  Specialists always win the min(output_dim)
    selection when applicable (fewer representatives → smaller output_dim).
    """

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        reps = _pair_orbit_atoms(spec.pf, plan)
        return [Type11PairEncoder(reps, spec.pf, plan)]


# ---------------------------------------------------------------------------
# Type21PairEncoder / Type21PairEncoderFactory  (TYPE_21)
# ---------------------------------------------------------------------------

class Type21PairEncoder(_BasePairEncoder):
    """
    TYPE_21 encoding: u-sign flips freely, v-sign is fixed.

    For each representative z = eval(+u)+i*eval(v), the full orbit element is
    {z, -conj(z)} since eval(-u) = -eval(u).  The polynomial of {z, -conj(z)}
    has even-degree coefficients real and odd-degree purely imaginary; only n
    independent reals are needed → 2n reals.
    """
    _METHOD_NAME     = "21"
    _NOTIONAL_FACTOR = 2

    def encode(self, event: dict) -> EncodingResult:
        z_reps = np.array(
            [complex(evaluate(u, event), evaluate(v, event)) for u, v in self._reps],
            dtype=complex,
        )
        z_full = np.concatenate([z_reps, -np.conj(z_reps)])
        coeffs, _, _ = zip_embed(z_full)
        return EncodingResult(
            values=_complex_to_reals(coeffs.imag[0::2] + 1j * coeffs.real[1::2]),
            metadata={"method": self._METHOD_NAME},
        )


class Type21PairEncoderFactory(PairOrbitEncoderFactory):
    """Returns Type21PairEncoder if the orbit partitions as TYPE_21; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        reps = _try_type21_partition(_pair_orbit_atoms(spec.pf, plan))
        if reps is None:
            return []
        return [Type21PairEncoder(reps, spec.pf, plan)]


# ---------------------------------------------------------------------------
# Type12PairEncoder / Type12PairEncoderFactory  (TYPE_12)
# ---------------------------------------------------------------------------

class Type12PairEncoder(_BasePairEncoder):
    """
    TYPE_12 encoding: v-sign flips freely, u-sign is fixed.

    For each representative z = eval(u)+i*eval(+v), the full orbit is
    {z, conj(z)} since eval(-v) = -eval(v).  The polynomial of {z, conj(z)}
    has real coefficients; n independent reals → 2n reals.
    """
    _METHOD_NAME     = "12"
    _NOTIONAL_FACTOR = 2

    def encode(self, event: dict) -> EncodingResult:
        z_reps = np.array(
            [complex(evaluate(u, event), evaluate(v, event)) for u, v in self._reps],
            dtype=complex,
        )
        z_full = np.concatenate([z_reps, np.conj(z_reps)])
        coeffs, _, _ = zip_embed(z_full)
        r = coeffs.real
        return EncodingResult(
            values=_complex_to_reals(r[0::2] + 1j * r[1::2]),
            metadata={"method": self._METHOD_NAME},
        )


class Type12PairEncoderFactory(PairOrbitEncoderFactory):
    """
    Returns Type12PairEncoder if the orbit partitions as TYPE_12; [] otherwise.

    Thin wrapper: TYPE_12 partition = TYPE_21 partition with u/v roles swapped.
    """

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        reps = _try_type12_partition(_pair_orbit_atoms(spec.pf, plan))
        if reps is None:
            return []
        return [Type12PairEncoder(reps, spec.pf, plan)]


# ---------------------------------------------------------------------------
# Type22PairEncoder / Type22PairEncoderFactory  (TYPE_22)
# ---------------------------------------------------------------------------

class Type22PairEncoder(_BasePairEncoder):
    """
    TYPE_22 encoding: both signs flip freely (all four combinations present).

    Invariant: z² for representative z = eval(+u)+i*eval(+v).
    The orbit of z² is {z², conj(z²)}; polynomial has real coefficients →
    n independent reals → 2n reals.
    """
    _METHOD_NAME     = "22"
    _NOTIONAL_FACTOR = 4

    def encode(self, event: dict) -> EncodingResult:
        z_reps = np.array(
            [complex(evaluate(u, event), evaluate(v, event)) for u, v in self._reps],
            dtype=complex,
        )
        w = z_reps ** 2
        w_full = np.concatenate([w, np.conj(w)])
        coeffs, _, _ = zip_embed(w_full)
        r = coeffs.real
        return EncodingResult(
            values=_complex_to_reals(r[0::2] + 1j * r[1::2]),
            metadata={"method": self._METHOD_NAME},
        )


class Type22PairEncoderFactory(PairOrbitEncoderFactory):
    """Returns Type22PairEncoder if the orbit partitions as TYPE_22; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        reps = _try_type22_partition(_pair_orbit_atoms(spec.pf, plan))
        if reps is None:
            return []
        return [Type22PairEncoder(reps, spec.pf, plan)]


# ---------------------------------------------------------------------------
# NegPairEncoder / NegPairEncoderFactory  (TYPE_NEG)
# ---------------------------------------------------------------------------

class NegPairEncoder(_BasePairEncoder):
    """
    TYPE_NEG encoding: correlated sign flip — both u and v negate together.

    Two orbit forms arise:
      Form 1 — {(+u,+v),(-u,-v)}: z = eval(+u)+i*eval(+v),  orbit = {z,-z}
      Form 2 — {(+u,-v),(-u,+v)}: z = eval(+u)+i*eval(-v),  orbit = {z,-z}

    In both cases the orbit of z is {z,-z}, so z² is invariant:
    (-z)² = z².  Embeds z² directly: n complex → 2n reals.

    Note: forming {z², conj(z²)} (as in TYPE_22) would be wrong — conj(z²) is
    not in the orbit under TYPE_NEG, and adding it loses Im(z²) = 2·Re(z)·Im(z).
    """
    _METHOD_NAME     = "neg"
    _NOTIONAL_FACTOR = 2

    def encode(self, event: dict) -> EncodingResult:
        z_reps = np.array(
            [complex(evaluate(u, event), evaluate(v, event)) for u, v in self._reps],
            dtype=complex,
        )
        w = z_reps ** 2
        coeffs, _, _ = zip_embed(w)
        return EncodingResult(
            values=_complex_to_reals(coeffs),
            metadata={"method": self._METHOD_NAME},
        )


class NegPairEncoderFactory(PairOrbitEncoderFactory):
    """Returns NegPairEncoder if the orbit partitions as TYPE_NEG; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        reps = _try_neg_partition(_pair_orbit_atoms(spec.pf, plan))
        if reps is None:
            return []
        return [NegPairEncoder(reps, spec.pf, plan)]


# ---------------------------------------------------------------------------
# Convenience: the full set of standard row-pair factories
# ---------------------------------------------------------------------------

def standard_row_pair_factories() -> list[PairOrbitEncoderFactory]:
    """
    Return one instance of each of the six standard row-pair factories in the
    recommended registration order.

    Specialists are listed before Type11PairEncoderFactory so that when a
    specialist and the fallback both offer encoders with the same output_dim,
    the specialist is taken (first in list wins ties in min() scan).

    Suitable for passing directly to OverlapBlockEncoderFactory:
        OverlapBlockEncoderFactory(standard_row_pair_factories())
    """
    return [
        SelfPairEncoderFactory(),
        NegPairEncoderFactory(),
        Type22PairEncoderFactory(),
        Type21PairEncoderFactory(),
        Type12PairEncoderFactory(),
        Type11PairEncoderFactory(),
    ]
