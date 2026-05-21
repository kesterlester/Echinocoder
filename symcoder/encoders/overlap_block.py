"""
symcoder.encoders.overlap_block
=================================
Overlap-block level of the Phase 2 encoding hierarchy.

An OVERLAP BLOCK is the set of all PairFlavours sharing the same
(op_u, flavour_u, op_v, flavour_v) — i.e. they differ only in their overlap
counts.

OverlapBlockSpec
    Wraps the list of PairFlavours that form one overlap block.

OverlapBlockEncoder
    A ready-to-use, bound encoder for one complete overlap block.  Holds the
    PairSelection list chosen during assess(), and delegates encode/describe
    to the selected row-pair encoders.

OverlapBlockEncoderFactory
    Receives a list of row-pair factories at construction time (no external
    registry object; the list is the registry).  assess() queries those
    factories for each pair in the block, applies a min(output_dim) selection,
    then optionally applies the complementarity drop (the non-null pair with
    the largest output_dim is excluded and recorded as NULL_COMP).

    Constructor keyword arg: use_complementarity_drop (default True).
    Set False for a debug/uncompressed variant.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Any

import numpy as np

from symcoder.describe import SegmentInfo, OverlapBlockNode, _assoc_example
from .pair_base import PairOrbitEncoder, PairOrbitEncoderFactory, PairOrbitSpec, EncodingResult


# ---------------------------------------------------------------------------
# OverlapBlockSpec
# ---------------------------------------------------------------------------

@dataclass
class OverlapBlockSpec:
    """
    Specification for one overlap block: the list of PairFlavours that share
    (op_u, flavour_u, op_v, flavour_v) and differ only in their overlap counts.
    """
    block: list  # list[PairFlavour]


# ---------------------------------------------------------------------------
# PairSelection  (internal record used by OverlapBlockEncoder)
# ---------------------------------------------------------------------------

@dataclass
class PairSelection:
    """
    Records the encoder chosen for one PairFlavour within an overlap block.

    is_comp_drop=True means this pair was excluded by the complementarity drop
    rule; its encoder is stored for metadata (output_dim, method_name) but
    encode() will never be called on it.
    """
    pf:           Any             # PairFlavour
    encoder:      PairOrbitEncoder
    is_comp_drop: bool = False


# ---------------------------------------------------------------------------
# OverlapBlockEncoder
# ---------------------------------------------------------------------------

class OverlapBlockEncoder:
    """
    Ready-to-use encoder for one overlap block.  Constructed by
    OverlapBlockEncoderFactory.assess() with all row-pair selections pre-made.

    encode(event) calls each selected (non-comp-dropped) row-pair encoder in
    order and concatenates their outputs.

    describe() returns a flat list of SegmentInfo, one per pair in the block,
    using:
      kind="NULL_SELF"  for null (self-pair) encoders
      kind="ASSOC"      for non-null, non-dropped encoders
      kind="NULL_COMP"  for the complementarity-dropped pair

    start values in describe() are relative to the start of this block (i.e.
    the first segment has start=0).  Phase2Encoder offsets them to absolute
    positions.
    """

    def __init__(self, selections: list[PairSelection], plan: Any) -> None:
        self._selections = selections
        self._plan       = plan

    @property
    def output_dim(self) -> int:
        return sum(
            sel.encoder.output_dim
            for sel in self._selections
            if not sel.is_comp_drop
        )

    def encode(self, event: dict) -> EncodingResult:
        parts = []
        for sel in self._selections:
            if sel.is_comp_drop:
                continue
            result = sel.encoder.encode(event)
            parts.append(result.values)
        values = np.concatenate(parts) if parts else np.array([], dtype=np.float64)
        return EncodingResult(values=values)

    def decode(self, values: np.ndarray, phase1_results: dict) -> list:
        """Decode this block's slice of the Phase 2 encoded array.

        Parameters
        ----------
        values : np.ndarray
            Sub-array produced by encode() for this block (non-comp-dropped data only).
        phase1_results : dict[(str, tuple), AnnotatedMultisetOfReals]
            Phase 1 decoded multisets keyed by (op_name, flavour_counts_tuple).
            Required only for NULL_SELF entries; ignored for ASSOC entries.

        Returns
        -------
        list[AnnotatedMultisetOfRealPairs]
            One entry per non-comp-dropped selection, in selection order.
            NULL_SELF entries are reconstructed from Phase 1; ASSOC entries are
            decoded by the corresponding row-pair encoder.
        """
        from symcoder.decoded_types import AnnotatedMultisetOfRealPairs

        results = []
        cursor = 0
        for sel in self._selections:
            if sel.is_comp_drop:
                continue
            enc = sel.encoder
            if enc.output_dim == 0:
                # NULL_SELF: (u=v) pairs are fully determined by Phase 1.
                # The pairs are (a, a) for each atom a in the shared orbit.
                key = (sel.pf.op_u.name, tuple(sel.pf.flavour_u.counts))
                phase1 = phase1_results[key]
                results.append(AnnotatedMultisetOfRealPairs(
                    pairs=[(v, v) for v in phase1.values],
                    atom_pairs=[(a, a) for a in phase1.atoms],
                ))
            else:
                chunk = values[cursor:cursor + enc.output_dim]
                results.append(enc.decode(chunk))
                cursor += enc.output_dim
        return results

    def describe(self) -> OverlapBlockNode:
        types  = self._plan.context.types
        segs   = []
        cursor = 0  # relative cursor within this block

        for sel in self._selections:
            pf = sel.pf
            if sel.is_comp_drop:
                segs.append(SegmentInfo(
                    kind            = "NULL_COMP",
                    start           = cursor,
                    length          = 0,
                    op_u            = pf.op_u.name,
                    flavour_u       = tuple(pf.flavour_u.counts),
                    op_v            = pf.op_v.name,
                    flavour_v       = tuple(pf.flavour_v.counts),
                    overlap         = tuple(pf.overlap),
                    symmetry_class  = sel.encoder.method_name,
                    notional_length = sel.encoder.output_dim,
                    method_name     = sel.encoder.method_name,
                    example         = _assoc_example(
                        pf.op_u.name, pf.flavour_u.counts,
                        pf.op_v.name, pf.flavour_v.counts,
                        pf.overlap, types,
                    ),
                ))
                # cursor does not advance — this pair contributes nothing
            else:
                enc_segs = sel.encoder.describe()
                for seg in enc_segs:
                    segs.append(dataclasses.replace(seg, start=seg.start + cursor))
                cursor += sel.encoder.output_dim

        pf0 = self._selections[0].pf
        return OverlapBlockNode(
            op_u      = pf0.op_u.name,
            flavour_u = tuple(pf0.flavour_u.counts),
            op_v      = pf0.op_v.name,
            flavour_v = tuple(pf0.flavour_v.counts),
            segments  = segs,
        )


# ---------------------------------------------------------------------------
# OverlapBlockEncoderFactory
# ---------------------------------------------------------------------------

class OverlapBlockEncoderFactory:
    """
    Factory for OverlapBlockEncoder instances.

    Constructor arguments
    ---------------------
    row_pair_factories : list[PairOrbitEncoderFactory]
        The row-pair factories to consult for each pair in a block.  These are
        the "inner registry" — pass them at construction rather than via a
        separate registry object.  Factories are queried in list order; all
        offered encoders are collected and the one with the smallest output_dim
        is selected per pair.

    use_complementarity_drop : bool (default True)
        When True (the "tight" mode), the non-null pair with the largest
        output_dim in the block is excluded (replaced by a NULL_COMP entry).
        When False (the "debug" mode), all non-null pairs are encoded, giving
        a redundant but fully inspectable output.

    Selection policy
    ----------------
    For each PairFlavour in the block:
      1. Query every row_pair_factory.assess(PairOrbitSpec(pf), plan).
      2. Collect all offered encoders.
      3. Select the one with the minimum output_dim.
         (For self-pairs, SelfPairEncoderFactory offers NullPairEncoder with
         output_dim=0, which always wins the minimum — this is how self-pair
         dropping is handled without special-case logic here.)

    Complementarity drop
    --------------------
    After selection, if use_complementarity_drop:
      Among selected encoders with output_dim > 0, find the one with the
      largest output_dim.  Mark it as is_comp_drop=True.  The overlap-block
      encoder will skip its encode() call and emit a NULL_COMP SegmentInfo.
    """

    def __init__(
        self,
        row_pair_factories: list[PairOrbitEncoderFactory],
        *,
        use_complementarity_drop: bool = True,
    ) -> None:
        self._row_pair_factories       = row_pair_factories
        self._use_complementarity_drop = use_complementarity_drop

    def assess(self, spec: OverlapBlockSpec, plan: Any) -> list[OverlapBlockEncoder]:
        """
        Produce a ready-to-use OverlapBlockEncoder for this block.

        Returns a one-element list (always succeeds as long as at least one
        row-pair factory handles each pair).  Returns [] only if no factory
        can handle any pair in the block (mis-configured factory list).
        """
        selections: list[PairSelection] = []

        for pf in spec.block:
            pair_spec = PairOrbitSpec.from_pair_flavour(pf)
            offered: list[PairOrbitEncoder] = []
            for factory in self._row_pair_factories:
                offered.extend(factory.assess(pair_spec, plan))

            if not offered:
                # No factory could handle this pair — return [] to signal failure.
                return []

            chosen = min(offered, key=lambda e: e.output_dim)
            selections.append(PairSelection(pf=pf, encoder=chosen))

        if self._use_complementarity_drop and selections:
            non_null = [
                (i, sel) for i, sel in enumerate(selections)
                if sel.encoder.output_dim > 0
            ]
            if non_null:
                comp_idx, _ = max(non_null, key=lambda x: x[1].encoder.output_dim)
                old = selections[comp_idx]
                selections[comp_idx] = PairSelection(
                    pf=old.pf, encoder=old.encoder, is_comp_drop=True
                )

        return [OverlapBlockEncoder(selections, plan)]
