"""
symcoder.encoders.phase2_encoder
==================================
Top-level Phase 2 encoding: iterate all overlap blocks, select a block encoder
for each, concatenate their outputs.

Phase2Encoder
    Ready-to-use encoder for all overlap blocks in one Plan.  Holds a list of
    (OverlapBlockSpec, OverlapBlockEncoder) pairs produced by Phase2EncoderFactory.
    encode(event) concatenates all block outputs.
    describe(start_offset) returns a flat list of SegmentInfo with absolute
    positions (Phase2Encoder adds start_offset so Phase 1 positions are respected).

Phase2EncoderFactory
    Takes a list of OverlapBlockEncoderFactory instances at construction (the
    "inner registry").  build(plan) enumerates all overlap blocks from the
    plan's repS / canonical_pair_flavours, asks each block factory to assess
    each block, and picks the first factory that returns a non-empty response.
"""
from __future__ import annotations

import dataclasses
from itertools import groupby
from typing import Any

import numpy as np

from symatom import repS
from symatom.rep import canonical_pair_flavours
from symcoder.describe import SegmentInfo, _block_key
from .pair_base import EncodingResult
from .overlap_block import OverlapBlockSpec, OverlapBlockEncoder, OverlapBlockEncoderFactory


class Phase2Encoder:
    """
    Ready-to-use encoder for all Phase 2 overlap blocks in one Plan.

    encode(event) calls each block encoder in order and concatenates outputs.
    describe(start_offset) returns a flat SegmentInfo list with absolute array
    positions: each block encoder's relative positions are shifted by the
    cumulative cursor (starting from start_offset, which accounts for Phase 1).
    """

    def __init__(self, block_encoders: list[tuple[OverlapBlockSpec, OverlapBlockEncoder]]) -> None:
        self._block_encoders = block_encoders

    @property
    def output_dim(self) -> int:
        return sum(enc.output_dim for _, enc in self._block_encoders)

    def encode(self, event: dict) -> EncodingResult:
        parts = []
        for _, enc in self._block_encoders:
            result = enc.encode(event)
            parts.append(result.values)
        values = np.concatenate(parts) if parts else np.array([], dtype=np.float64)
        return EncodingResult(values=values)

    def describe(self, start_offset: int = 0) -> list[SegmentInfo]:
        segs   = []
        cursor = start_offset
        for _, enc in self._block_encoders:
            block_segs = enc.describe()  # relative positions (start=0 for first)
            for seg in block_segs:
                segs.append(dataclasses.replace(seg, start=seg.start + cursor))
            cursor += enc.output_dim
        return segs


class Phase2EncoderFactory:
    """
    Factory that builds a Phase2Encoder for a given Plan.

    Constructor arguments
    ---------------------
    block_factories : list[OverlapBlockEncoderFactory]
        The overlap-block factories to consult for each block.  Factories are
        tried in list order; the first one that returns a non-empty assess()
        response is used for that block.

    Usage
    -----
        phase2_factory = Phase2EncoderFactory([
            OverlapBlockEncoderFactory(standard_row_pair_factories())
        ])
        phase2_enc = phase2_factory.build(plan)
        values = phase2_enc.encode(event)
        segs   = phase2_enc.describe(start_offset=phase1_output_dim)
    """

    def __init__(self, block_factories: list[OverlapBlockEncoderFactory]) -> None:
        self._block_factories = block_factories

    def build(self, plan: Any) -> Phase2Encoder:
        """
        Enumerate all overlap blocks from plan's pair structure, select a block
        encoder for each, and return a ready-to-use Phase2Encoder.
        """
        fo_list = repS(plan.context, plan.operations)
        pf_list = canonical_pair_flavours(fo_list, plan.context)

        block_encoders: list[tuple[OverlapBlockSpec, OverlapBlockEncoder]] = []

        for _bk, block_iter in groupby(pf_list, key=_block_key):
            block = list(block_iter)
            spec  = OverlapBlockSpec(block=block)

            chosen_enc: OverlapBlockEncoder | None = None
            for factory in self._block_factories:
                encoders = factory.assess(spec, plan)
                if encoders:
                    chosen_enc = encoders[0]
                    break

            if chosen_enc is not None:
                block_encoders.append((spec, chosen_enc))

        return Phase2Encoder(block_encoders)
