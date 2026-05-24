"""
symcoder.encoders.orbit_encoder
================================
OrbitEncoder / OrbitEncoderFactory: top-level Phase 1 encoder.

Mirrors Phase2Encoder / Phase2EncoderFactory: constructor-injected sub-factories
replace the mutable AtomOrbitEncoderRegistry pattern.

Usage
-----
    factory = OrbitEncoderFactory([SortEncoderFactory(), PolyEncoderFactory()])
    orbit_enc = factory.build(plan)
    values = orbit_enc.encode(event).values
    segs   = orbit_enc.describe(start_offset=0)
"""
from __future__ import annotations

import numpy as np
from symatom import repS

from ._base import AtomOrbitEncoderFactory, AtomOrbitEncoder, OrbitSpec, EncodingResult
from ..describe import SegmentInfo, Phase1Tree, _orbit_example


class OrbitEncoder:
    """
    Encodes all FlavouredOperator orbits for a plan.

    Built by OrbitEncoderFactory.build(plan).  Holds one chosen AtomOrbitEncoder
    per FlavouredOperator (selected as min output_dim across all offered encoders).
    """

    def __init__(self, selections: list, types) -> None:
        # selections: list of (FlavouredOperator, AtomOrbitEncoder)
        self._selections = selections
        self._types = types

    @property
    def output_dim(self) -> int:
        return sum(enc.output_dim for _fo, enc in self._selections)

    def encode(self, event: dict) -> EncodingResult:
        parts = [enc.encode(event).values for _fo, enc in self._selections]
        values = np.concatenate(parts) if parts else np.array([], dtype=np.float64)
        return EncodingResult(values=values)

    def describe(self, start_offset: int = 0) -> Phase1Tree:
        orbits = []
        cursor = start_offset
        for fo, enc in self._selections:
            length = enc.output_dim
            orbits.append(SegmentInfo(
                kind             = "ORBIT",
                start            = cursor,
                length           = length,
                op_u             = fo.operation,
                flavour_u        = tuple(fo.flavour.counts),
                method_name      = enc.method_name,
                notional_length  = enc.notional_output_dim,
                example          = _orbit_example(fo.operation.name, fo.flavour.counts, self._types),
            ))
            cursor += length
        return Phase1Tree(orbits=orbits)


class OrbitEncoderFactory:
    """
    Builds an OrbitEncoder for a plan from constructor-injected sub-factories.

    Sub-factories are tried in order; all offered encoders are collected and the
    one with the smallest output_dim is selected (matching the Phase 2 selection
    policy).  Raises RuntimeError if no factory can handle a given orbit.
    """

    def __init__(self, factories: list[AtomOrbitEncoderFactory]) -> None:
        self._factories = list(factories)

    def build(self, plan) -> OrbitEncoder:
        fo_list = repS(plan.context, plan.operations)
        selections = []
        for fo in fo_list:
            assert fo.count_of_atoms_one_per_sign() > 0
            spec = OrbitSpec.from_flavoured_operator(fo)
            offered = []
            for factory in self._factories:
                offered.extend(factory.assess(spec, plan))
            if not offered:
                raise RuntimeError(
                    f"No factory in OrbitEncoderFactory could handle orbit {fo!r}. "
                    f"({len(self._factories)} factory/factories registered.)"
                )
            selections.append((fo, min(offered, key=lambda enc: enc.output_dim)))
        return OrbitEncoder(selections, plan.context.types)
