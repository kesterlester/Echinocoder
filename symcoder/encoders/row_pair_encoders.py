"""
symcoder.encoders.row_pair_encoders
=====================================
Concrete row-pair encoders and their factories for Phase 2 encoding.

Each factory's assess() checks the PairFlavour's sign correlation type (via
_sign_correlation_type_from_pf) and returns [] for pairs it cannot handle.
An overlap-block manager queries all factories for each pair in its block and
selects among the offered encoders.

Factories and what they handle
--------------------------------
SelfPairEncoderFactory   self-pairs (any type)    → NullPairEncoder   output_dim=0
SSPairEncoderFactory     TYPE_11 pairs            → SSPairEncoder      output_dim=2n
SAPairEncoderFactory     TYPE_12 pairs            → SAPairEncoder      output_dim=2n
ASPairEncoderFactory     TYPE_21 pairs            → ASPairEncoder      output_dim=2n
AAPairEncoderFactory     TYPE_22 pairs            → AAPairEncoder      output_dim=2n
NEGPairEncoderFactory    TYPE_NEG pairs           → NEGPairEncoder     output_dim=2n

Here n = pf.count(type_sizes) — the number of base atom-pairs in the orbit.

All five non-null encoders produce the same output_dim = 2n.  The overlap-block
manager picks among offered encoders, typically by minimum output_dim; for ties
it takes the first in the list (favouring the more specialised factory that was
registered earlier).

Self-pair null encoding
-----------------------
Phase 1 stores every individual row orbit.  A self-pair (u=v) evaluation z_k
is fully determined by Phase 1 alone, so encoding it adds no information.
SelfPairEncoderFactory therefore returns NullPairEncoder (output_dim=0) for
self-pairs and [] for all other pairs.  Because output_dim=0 is the minimum
possible, the overlap-block manager's min(output_dim) selection naturally picks
this encoder whenever it is offered.
"""
from __future__ import annotations

import dataclasses
from abc import abstractmethod

import numpy as np

from symcoder.sign_correlation import (
    SignCorrelationType,
    _sign_correlation_type_from_pf,
    _complex_to_reals,
)
from symcoder.pairs import eval_pair_orbit_positive, _is_self_pair
from symcoder.describe import SegmentInfo, _assoc_example, _symmetry_class
from ._embedder import zip_embed
from .pair_base import PairOrbitEncoder, PairOrbitEncoderFactory, PairOrbitSpec, EncodingResult

_PRIORITY = 0.5


# ---------------------------------------------------------------------------
# _BasePairEncoder — shared scaffolding for all row-pair encoders
# ---------------------------------------------------------------------------

class _BasePairEncoder(PairOrbitEncoder):
    """
    Base class for concrete row-pair encoders.

    Subclasses set _METHOD_NAME and implement encode().  The describe() method
    generates an ASSOC SegmentInfo with start=0 (a placeholder; the
    overlap-block encoder adjusts it to the correct absolute position).
    NullPairEncoder overrides both describe() and output_dim.

    _NOTIONAL_FACTOR is the ratio of the SS-equivalent size to the actual output
    size, i.e. |achievable_set|.  SS=1, SA=AS=NEG=2, AA=4.  Used to populate
    notional_length in describe() so the full column shows the uncompressed size.
    """
    _METHOD_NAME:     str = "?"
    _NOTIONAL_FACTOR: int = 1

    def __init__(self, pf, plan) -> None:
        self._pf         = pf
        self._plan       = plan
        self._type_sizes = tuple(g.size for g in plan.context.types)
        self._n          = pf.count(self._type_sizes)

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
            kind           = "ASSOC",
            start          = 0,   # placeholder; overlap-block encoder adjusts
            length         = self.output_dim,
            op_u           = self._pf.op_u.name,
            flavour_u      = tuple(self._pf.flavour_u.counts),
            op_v           = self._pf.op_v.name,
            flavour_v      = tuple(self._pf.flavour_v.counts),
            overlap        = tuple(self._pf.overlap),
            symmetry_class = _symmetry_class(self._pf),
            notional_length= self._NOTIONAL_FACTOR * self.output_dim,
            method_name    = self.method_name,
            example        = _assoc_example(
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

    A self-pair (u=v) adds no information beyond what Phase 1 already stores
    (Phase 1 encodes every individual row orbit, from which any (u,u) pair is
    recoverable).  This encoder signals that fact by contributing no values.
    """
    _METHOD_NAME = "null_self"

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
            symmetry_class  = _symmetry_class(self._pf),
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
# SSPairEncoder / SSPairEncoderFactory  (TYPE_11)
# ---------------------------------------------------------------------------

class SSPairEncoder(_BasePairEncoder):
    """
    TYPE_11 encoding: achievable set = {(+1,+1)}, no sign changes.

    Embeds z_pos directly: n roots → n complex polynomial coefficients → 2n reals.
    """
    _METHOD_NAME = "ss"

    def encode(self, event: dict) -> EncodingResult:
        z_pos = np.array(eval_pair_orbit_positive(self._pf, self._plan, event), dtype=complex)
        assert len(z_pos) == self._n
        coeffs, _, _ = zip_embed(z_pos)
        return EncodingResult(
            values=_complex_to_reals(coeffs),
            metadata={"method": self._METHOD_NAME},
        )


class SSPairEncoderFactory(PairOrbitEncoderFactory):
    """Returns SSPairEncoder for TYPE_11 pairs; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        if _sign_correlation_type_from_pf(spec.pf) == SignCorrelationType.TYPE_11:
            return [SSPairEncoder(spec.pf, plan)]
        return []


# ---------------------------------------------------------------------------
# SAPairEncoder / SAPairEncoderFactory  (TYPE_12: v-sign flips freely)
# ---------------------------------------------------------------------------

class SAPairEncoder(_BasePairEncoder):
    """
    TYPE_12 encoding: achievable set = {(+1,+1),(+1,-1)}, only v's sign flips.

    Orbit: {z_k, conj(z_k)}  (conjugate-closed, 2n elements).
    Polynomial has real coefficients; packs n independent reals as n complex → 2n reals.
    """
    _METHOD_NAME     = "sa"
    _NOTIONAL_FACTOR = 2

    def encode(self, event: dict) -> EncodingResult:
        z_pos = np.array(eval_pair_orbit_positive(self._pf, self._plan, event), dtype=complex)
        assert len(z_pos) == self._n
        z_full = np.concatenate([z_pos, np.conj(z_pos)])
        coeffs, _, _ = zip_embed(z_full)
        r = coeffs.real
        return EncodingResult(
            values=_complex_to_reals(r[0::2] + 1j * r[1::2]),
            metadata={"method": self._METHOD_NAME},
        )


class SAPairEncoderFactory(PairOrbitEncoderFactory):
    """Returns SAPairEncoder for TYPE_12 pairs; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        if _sign_correlation_type_from_pf(spec.pf) == SignCorrelationType.TYPE_12:
            return [SAPairEncoder(spec.pf, plan)]
        return []


# ---------------------------------------------------------------------------
# ASPairEncoder / ASPairEncoderFactory  (TYPE_21: u-sign flips freely)
# ---------------------------------------------------------------------------

class ASPairEncoder(_BasePairEncoder):
    """
    TYPE_21 encoding: achievable set = {(+1,+1),(-1,+1)}, only u's sign flips.

    Orbit: {z_k, -conj(z_k)}  (anti-conjugate-closed, 2n elements).
    Even-degree polynomial coefficients are real; odd-degree are purely imaginary.
    Packs as: coeffs.imag[0::2] + 1j*coeffs.real[1::2] → n complex → 2n reals.
    """
    _METHOD_NAME     = "as"
    _NOTIONAL_FACTOR = 2

    def encode(self, event: dict) -> EncodingResult:
        z_pos = np.array(eval_pair_orbit_positive(self._pf, self._plan, event), dtype=complex)
        assert len(z_pos) == self._n
        z_full = np.concatenate([z_pos, -np.conj(z_pos)])
        coeffs, _, _ = zip_embed(z_full)
        return EncodingResult(
            values=_complex_to_reals(coeffs.imag[0::2] + 1j * coeffs.real[1::2]),
            metadata={"method": self._METHOD_NAME},
        )


class ASPairEncoderFactory(PairOrbitEncoderFactory):
    """Returns ASPairEncoder for TYPE_21 pairs; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        if _sign_correlation_type_from_pf(spec.pf) == SignCorrelationType.TYPE_21:
            return [ASPairEncoder(spec.pf, plan)]
        return []


# ---------------------------------------------------------------------------
# AAPairEncoder / AAPairEncoderFactory  (TYPE_22: both signs flip freely)
# ---------------------------------------------------------------------------

class AAPairEncoder(_BasePairEncoder):
    """
    TYPE_22 encoding: achievable set = all four sign combinations.

    Invariant: {z_k², conj(z_k²)}  (conjugate-closed, 2n elements).
    Polynomial has real coefficients; packs n independent reals as n complex → 2n reals.
    """
    _METHOD_NAME     = "aa"
    _NOTIONAL_FACTOR = 4

    def encode(self, event: dict) -> EncodingResult:
        z_pos = np.array(eval_pair_orbit_positive(self._pf, self._plan, event), dtype=complex)
        assert len(z_pos) == self._n
        w = z_pos ** 2
        w_full = np.concatenate([w, np.conj(w)])
        coeffs, _, _ = zip_embed(w_full)
        r = coeffs.real
        return EncodingResult(
            values=_complex_to_reals(r[0::2] + 1j * r[1::2]),
            metadata={"method": self._METHOD_NAME},
        )


class AAPairEncoderFactory(PairOrbitEncoderFactory):
    """Returns AAPairEncoder for TYPE_22 pairs; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        if _sign_correlation_type_from_pf(spec.pf) == SignCorrelationType.TYPE_22:
            return [AAPairEncoder(spec.pf, plan)]
        return []


# ---------------------------------------------------------------------------
# NEGPairEncoder / NEGPairEncoderFactory  (TYPE_NEG: correlated negation only)
# ---------------------------------------------------------------------------

class NEGPairEncoder(_BasePairEncoder):
    """
    TYPE_NEG encoding: achievable set = {(+1,+1),(-1,-1)}, correlated negation.

    Invariant: {z_k²}  (n elements, since (-z_k)² = z_k²).
    Embeds z_pos² directly: n roots → n complex coefficients → 2n reals.

    Note: forming {z_k², conj(z_k²)} (as in TYPE_22) would be wrong here —
    conj(z_k²) is NOT in the orbit under TYPE_NEG, and adding it loses
    Im(z_k²) = 2·Re(z_k)·Im(z_k).
    """
    _METHOD_NAME     = "neg"
    _NOTIONAL_FACTOR = 2

    def encode(self, event: dict) -> EncodingResult:
        z_pos = np.array(eval_pair_orbit_positive(self._pf, self._plan, event), dtype=complex)
        assert len(z_pos) == self._n
        w = z_pos ** 2
        coeffs, _, _ = zip_embed(w)
        return EncodingResult(
            values=_complex_to_reals(coeffs),
            metadata={"method": self._METHOD_NAME},
        )


class NEGPairEncoderFactory(PairOrbitEncoderFactory):
    """Returns NEGPairEncoder for TYPE_NEG pairs; [] otherwise."""

    def assess(self, spec: PairOrbitSpec, plan) -> list[PairOrbitEncoder]:
        if _is_self_pair(spec.pf):
            return []
        if _sign_correlation_type_from_pf(spec.pf) == SignCorrelationType.TYPE_NEG:
            return [NEGPairEncoder(spec.pf, plan)]
        return []


# ---------------------------------------------------------------------------
# Convenience: the full set of standard row-pair factories
# ---------------------------------------------------------------------------

def standard_row_pair_factories() -> list[PairOrbitEncoderFactory]:
    """
    Return one instance of each of the six standard row-pair factories in the
    recommended registration order (most specialised first so that min(output_dim)
    selection ties go to the specialist).

    Suitable for passing directly to OverlapBlockEncoderFactory:
        OverlapBlockEncoderFactory(standard_row_pair_factories())
    """
    return [
        SelfPairEncoderFactory(),
        SSPairEncoderFactory(),
        SAPairEncoderFactory(),
        ASPairEncoderFactory(),
        AAPairEncoderFactory(),
        NEGPairEncoderFactory(),
    ]
