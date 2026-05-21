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
from symatom.atoms import Atom


# ---------------------------------------------------------------------------
# Module-level helpers for the decode-planning pass
# ---------------------------------------------------------------------------

def _neg(atom: Atom) -> Atom:
    """Return a copy of *atom* with its sign negated."""
    return Atom(atom.operation, atom.labels, sign=-atom.sign)


def _plan_owned_atom_pairs(spec, selections, plan) -> tuple:
    """
    Decode-planning pass: assign every atom-pair in the block's full Cartesian
    product to exactly one selection.

    Algorithm
    ---------
    1. Build u_orbit = G-orbit of the u-side canonical seed atom.
       Build v_orbit = G-orbit of the v-side canonical seed atom.
       L = u_orbit × v_orbit  (full Cartesian product as a set).
    2. For each selection (in order):
       a. Obtain the base G-orbit via TheGroup.orbit_pair(u_canon, v_canon).
       b. Claim those pairs from L (discard + record in owned list).
       c. Try all 7 non-identity sign/swap transforms of (u_canon, v_canon):
            (neg_u, v), (u, neg_v), (neg_u, neg_v),
            (v, u), (neg_v, u), (v, neg_u), (neg_v, neg_u)
          For each transform that still has unclaimed pairs in L, claim them
          and record the transform tag.
    3. Assert L is empty (partition identity).

    The swap transforms are a no-op for cross-type blocks (op_u ≠ op_v) because
    the swapped pairs live in op_v × op_u space, which is disjoint from L.

    Returns
    -------
    tuple[list[list[tuple[Atom, Atom]]], list[list[str]]]
        (owned_lists, tag_lists) where:
          owned_lists[i] = all (Atom, Atom) pairs owned by selections[i].
          tag_lists[i]   = list of value-transform tags for partner orbits
                           claimed by selections[i] (does not include the base
                           orbit, which corresponds to the encoder's own pairs).
        Union of owned_lists = full Cartesian product; sets are disjoint.
    """
    context = plan.context
    pf0 = spec.block[0]

    # Build canonical seed atoms for the u-side and v-side of this block.
    # All PairFlavours in the block share the same (op_u, flavour_u, op_v,
    # flavour_v), so pf0 is representative.
    u_labels: list = []
    for g, ku in zip(context.types, pf0.flavour_u.counts):
        u_labels.extend(g.labels[:ku])
    u_seed = Atom(pf0.op_u, tuple(u_labels), sign=+1)

    v_labels: list = []
    for g, kv in zip(context.types, pf0.flavour_v.counts):
        v_labels.extend(g.labels[:kv])
    v_seed = Atom(pf0.op_v, tuple(v_labels), sign=+1)

    u_orbit = context.the_group.orbit(u_seed)
    v_orbit = context.the_group.orbit(v_seed)

    # L = full Cartesian product as a set of (Atom, Atom) pairs.
    L: set = {(u, v) for u in u_orbit for v in v_orbit}

    owned_lists: list = [[] for _ in selections]
    tag_lists:   list = [[] for _ in selections]

    for i, sel in enumerate(selections):
        u_can, v_can = sel.pf.canonical_pair(context)

        # --- Claim the base G-orbit ---
        base_orbit = context.the_group.orbit_pair(u_can, v_can)
        for pair in base_orbit:
            if pair in L:
                L.discard(pair)
                owned_lists[i].append(pair)

        # --- Try all non-identity sign/swap transforms ---
        # Sign-negation is only valid for ANTISYMMETRIC operations (sign=+1 is
        # the only valid sign for SYMMETRIC operations).
        from symatom import ArgumentSymmetry
        u_anti = (u_can.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)
        v_anti = (v_can.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC)

        # Each entry: (tag, (u_transform_atom, v_transform_atom))
        # The tag encodes the value-space transform to apply to base decoded pairs.
        named_transforms = []
        # Pure sign flips
        if u_anti:
            named_transforms.append(("neg_u",        (_neg(u_can), v_can        )))
        if v_anti:
            named_transforms.append(("neg_v",        (u_can,       _neg(v_can)  )))
        if u_anti and v_anti:
            named_transforms.append(("neg_both",     (_neg(u_can), _neg(v_can)  )))
        # Swap (u,v) → (v,u): only reaches L when op_u == op_v (same-type blocks)
        named_transforms.append(    ("swap",         (v_can,       u_can        )))
        if v_anti:
            named_transforms.append(("swap_neg_v",   (_neg(v_can), u_can        )))
        if u_anti:
            named_transforms.append(("swap_neg_u",   (v_can,       _neg(u_can)  )))
        if u_anti and v_anti:
            named_transforms.append(("swap_neg_both",(_neg(v_can), _neg(u_can)  )))

        for tag, (u_t, v_t) in named_transforms:
            partner_orbit = context.the_group.orbit_pair(u_t, v_t)
            new_pairs = [p for p in partner_orbit if p in L]
            if new_pairs:
                for pair in new_pairs:
                    L.discard(pair)
                    owned_lists[i].append(pair)
                tag_lists[i].append(tag)

    assert len(L) == 0, (
        f"_plan_owned_atom_pairs: {len(L)} atom-pair(s) in block "
        f"({pf0.op_u.name}×{pf0.op_v.name}) not claimed by any selection — "
        "partition identity violated"
    )
    return owned_lists, tag_lists


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

    owned_atom_pairs is populated by the decode-planning pass
    (_plan_owned_atom_pairs) inside OverlapBlockEncoderFactory.assess().  It
    lists every (Atom, Atom) pair in the block's full Cartesian product that
    belongs to this selection: the base G-orbit of pf.canonical_pair PLUS any
    NULL_SIGN_OR_SWAP partner orbits claimed by this selection.  The union
    over all selections equals the full Cartesian product (partition identity).
    Used only in decode(); encode() is unchanged.
    """
    pf:               Any             # PairFlavour
    encoder:          PairOrbitEncoder
    is_comp_drop:     bool = False
    owned_atom_pairs:          list = dataclasses.field(default_factory=list)
    owned_value_transform_tags: list = dataclasses.field(default_factory=list)


# ---------------------------------------------------------------------------
# Value-space transforms for partner-orbit subtraction in decode()
# ---------------------------------------------------------------------------

_VALUE_TRANSFORMS: dict = {
    # tag          → function (u_val, v_val) → (new_u, new_v)
    "neg_u":        lambda u, v: (-u,  v),
    "neg_v":        lambda u, v: ( u, -v),
    "neg_both":     lambda u, v: (-u, -v),
    "swap":         lambda u, v: ( v,  u),
    "swap_neg_v":   lambda u, v: (-v,  u),
    "swap_neg_u":   lambda u, v: ( v, -u),
    "swap_neg_both":lambda u, v: (-v, -u),
}


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
            One entry per selection (including the comp-dropped one), in selection
            order.  NULL_SELF entries are reconstructed from Phase 1; ASSOC entries
            are decoded by the corresponding row-pair encoder; the NULL_COMP entry
            is reconstructed by multiset complement (Z_full ∖ NULL_SELF ∖ ASSoC).
        """
        from collections import Counter
        from symcoder.decoded_types import AnnotatedMultisetOfRealPairs

        results = []
        comp_drop_idx = None   # position of the NULL_COMP slot in results
        cursor = 0

        for sel in self._selections:
            if sel.is_comp_drop:
                # Record the position; fill with None as placeholder.
                comp_drop_idx = len(results)
                results.append(None)
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
                u_key = (sel.pf.op_u.name, tuple(sel.pf.flavour_u.counts))
                v_key = (sel.pf.op_v.name, tuple(sel.pf.flavour_v.counts))
                u_ph1 = phase1_results.get(u_key)
                v_ph1 = phase1_results.get(v_key)
                chunk = values[cursor:cursor + enc.output_dim]
                results.append(enc.decode(chunk, u_phase1=u_ph1, v_phase1=v_ph1))
                cursor += enc.output_dim

        if comp_drop_idx is not None:
            # Reconstruct the NULL_COMP association by multiset complement.
            #
            # All selections in a block share (op_u, flavour_u, op_v, flavour_v),
            # so we can use any selection's pf to look up U and V.
            pf0  = self._selections[0].pf
            u_key = (pf0.op_u.name, tuple(pf0.flavour_u.counts))
            v_key = (pf0.op_v.name, tuple(pf0.flavour_v.counts))
            U = phase1_results[u_key].values
            V = phase1_results[v_key].values

            # Build full Cartesian product as a Counter.
            # After polishing, all decoded values are bit-for-bit copies of Phase 1
            # values, so exact Counter subtraction is always correct.
            assert not any(u != u for u in U), "NaN in Phase 1 U values"
            assert not any(v != v for v in V), "NaN in Phase 1 V values"
            Z_full: Counter = Counter()
            for u in U:
                for v in V:
                    Z_full[(u, v)] += 1

            # Subtract every active selection's owned value pairs from Z_full.
            #
            # Each selection's encoder returns the BASE G-orbit value pairs in
            # results[i].pairs.  Partner orbits (NULL_SIGN_OR_SWAP variants)
            # are owned by the same selection and their value pairs are obtained
            # by applying the recorded value-space transforms to the base pairs.
            #
            # _VALUE_TRANSFORMS maps tag → (u,v) → (u',v').  The tags were
            # computed at assess() time by _plan_owned_atom_pairs and stored in
            # sel.owned_value_transform_tags.  Only transforms that actually
            # claimed new atom pairs get a tag, so no double-subtraction occurs.
            for i, sel in enumerate(self._selections):
                if sel.is_comp_drop:
                    continue
                base_pairs = results[i].pairs
                # Subtract base orbit pairs.
                for pair in base_pairs:
                    Z_full[pair] -= 1
                # Subtract partner orbit pairs via value-space transforms.
                for tag in sel.owned_value_transform_tags:
                    fn = _VALUE_TRANSFORMS[tag]
                    for u_val, v_val in base_pairs:
                        Z_full[fn(u_val, v_val)] -= 1

            assert all(cnt >= 0 for cnt in Z_full.values()), (
                "NULL_COMP reconstruction: Counter went negative — "
                "partition identity violated (bug in encoder or polishing)"
            )

            null_comp_pairs = [
                pair for pair, cnt in Z_full.items() for _ in range(cnt)
            ]

            results[comp_drop_idx] = AnnotatedMultisetOfRealPairs(
                pairs=null_comp_pairs,
                atom_pairs=list(self._selections[comp_drop_idx].owned_atom_pairs),
            )

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

        # Decode-planning pass: assign owned atom-pairs to each selection.
        # This runs AFTER the complementarity drop so that is_comp_drop flags
        # are already set (the planning pass itself does not use is_comp_drop;
        # all selections receive their owned_atom_pairs regardless of drop status).
        owned_lists, tag_lists = _plan_owned_atom_pairs(spec, selections, plan)
        selections = [
            dataclasses.replace(sel, owned_atom_pairs=owned,
                                owned_value_transform_tags=tags)
            for sel, owned, tags in zip(selections, owned_lists, tag_lists)
        ]

        return [OverlapBlockEncoder(selections, plan)]
