from __future__ import annotations
from dataclasses import dataclass
from itertools import groupby
from symatom.atoms import ArgumentSymmetry
from symatom.rep import canonical_pair_flavours
from symatom import repL


def _symmetry_class(pf) -> str:
    u = "A" if pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC else "S"
    v = "A" if pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC else "S"
    return f"{u}{v}"


def _block_key(pf):
    """Key that identifies an OVERLAP BLOCK: fixed (op_u, flavour_u, op_v, flavour_v)."""
    u = (pf.op_u.name, pf.op_u.rank, pf.op_u.parity,
         pf.op_u.argument_symmetry.value, pf.flavour_u.counts)
    v = (pf.op_v.name, pf.op_v.rank, pf.op_v.parity,
         pf.op_v.argument_symmetry.value, pf.flavour_v.counts)
    return (u, v)


@dataclass(frozen=True)
class SegmentInfo:
    """
    Metadata for one entry in the describe_encoding() output.

    kind
    ----
    "ORBIT"  Phase 1: sort-encoded single-operator orbit.
             Fields op_v, flavour_v, overlap, symmetry_class are None.
             sign_compressed=True: ANTISYMMETRIC operation — the orbit contains
             {+u, -u} pairs; only |eval(sign=+1 atom)| is stored per combination
             (n values, n = base label combinations).
             sign_compressed=False: SYMMETRIC/UNSTRUCTURED — eval(u) for every
             atom stored (fo.count() values).

    "ASSOC"  Phase 2: compressed pair encoding of one association.
             One PairFlavour within an OVERLAP BLOCK, encoded via _embed_compressed.
             All fields are set.

    "NULL"   Phase 2: dropped association — deducible from orbit headers + other
             associations in the same OVERLAP BLOCK by complementation in the
             association table.  length == 0; all fields are set.
             These are the largest-count association per OVERLAP BLOCK.

    Fields
    ------
    start           : first index (inclusive) in the output array (real units)
    length          : number of real values in this segment (0 for NULL).
                      ORBIT: fo.count() // 2 if sign_compressed else fo.count().
                      ASSOC: two reals per association pair (2 * pf.count()),
                             since each complex polynomial coefficient = (re, im).
    op_u            : name of operation u (or the single operation for ORBIT)
    flavour_u       : counts tuple — labels from each group going into op_u
    op_v            : name of operation v (None for ORBIT)
    flavour_v       : counts tuple for op_v (None for ORBIT)
    overlap         : shared-label counts per group (None for ORBIT)
    symmetry_class  : "SS", "SA", "AS", or "AA" (None for ORBIT)
    sign_compressed : True if ANTISYMMETRIC ORBIT (5c compression applied),
                      False if SYMMETRIC/UNSTRUCTURED ORBIT, None for ASSOC/NULL.

    Human-readable display
    ----------------------
    str(seg)  — single-line summary.

    Machine-readable export
    -----------------------
    seg.to_dict()  — plain dict; json.dumps-able for downstream tools.

    Notes
    -----
    NULL entries appear at their logical position in the OVERLAP BLOCK sequence
    but contribute zero coefficients; start == stop for a NULL entry.
    Consecutive ASSOC/NULL entries within the same block are ordered by
    canonical_pair_flavours' overlap sort (lex ascending).
    """
    kind:            str
    start:           int
    length:          int
    op_u:            str
    flavour_u:       tuple
    op_v:            str   | None = None
    flavour_v:       tuple | None = None
    overlap:         tuple | None = None
    symmetry_class:  str   | None = None
    sign_compressed: bool  | None = None

    @property
    def stop(self) -> int:
        return self.start + self.length

    def __str__(self) -> str:
        idx = f"[{self.start}:{self.stop}]"
        fl_u = ",".join(str(c) for c in self.flavour_u)

        if self.kind == "ORBIT":
            sc = "  [sign_compressed]" if self.sign_compressed else ""
            return f"{idx}  {self.op_u}  ORBIT  u=({fl_u})  len={self.length}{sc}"

        fl_v = ",".join(str(c) for c in self.flavour_v)
        ov   = ",".join(str(c) for c in self.overlap)
        base = (
            f"{idx}"
            f"  {self.op_u}×{self.op_v}"
            f"  {self.symmetry_class}"
            f"  u=({fl_u})  v=({fl_v})  shared=({ov})"
            f"  len={self.length}"
        )
        if self.kind == "NULL":
            base += "  NULL_ENCODING(deducible_from_uv_overlap_block)"
        return base

    def to_dict(self) -> dict:
        d = {
            "kind":      self.kind,
            "start":     self.start,
            "stop":      self.stop,
            "length":    self.length,
            "op_u":      self.op_u,
            "flavour_u": list(self.flavour_u),
        }
        if self.kind == "ORBIT":
            d["sign_compressed"] = self.sign_compressed
        else:
            d.update({
                "op_v":           self.op_v,
                "flavour_v":      list(self.flavour_v),
                "overlap":        list(self.overlap),
                "symmetry_class": self.symmetry_class,
            })
        return d


def describe_encoding(plan) -> list[SegmentInfo]:
    """
    Return a list of SegmentInfo objects describing the full structure of the
    vector produced by encode(plan, event).

    The list mirrors encode()'s two-phase structure exactly:

      Phase 1 — ORBIT segments (one per distinct FlavouredOperator in repL):
        kind="ORBIT", length=ceil(fo.count()/2).  Appear first.

      Phase 2 — ASSOC and NULL segments (one per PairFlavour):
        Within each OVERLAP BLOCK (fixed op_u/flavour_u/op_v/flavour_v) the
        largest-count PairFlavour is marked kind="NULL" (length=0, cursor does
        not advance); all others are kind="ASSOC".

    NULL segments present as real encoding entries with zero length — their
    position in the block sequence is preserved so a decoder can reconstruct
    which association was dropped.

    Usage examples
    --------------
    # Human-readable table:
    for seg in describe_encoding(plan):
        print(seg)

    # Machine-readable JSON:
    import json
    json.dumps([s.to_dict() for s in describe_encoding(plan)])

    # Total encoded length (should match len(encode(plan, event))):
    sum(s.length for s in describe_encoding(plan))

    This function is pure — it does not evaluate any event data.
    """
    fo_list = repL(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    group_sizes = tuple(g.size for g in plan.context.groups)

    segments = []
    cursor = 0

    # Phase 1: one ORBIT segment per distinct FlavouredOperator.
    seen_fo_keys: set = set()
    for fo in fo_list:
        fo_key = (fo.operation.name, fo.operation.rank, fo.operation.parity,
                  fo.operation.argument_symmetry.value, fo.flavour.counts)
        if fo_key in seen_fo_keys:
            continue
        seen_fo_keys.add(fo_key)
        if fo.count() == 0:
            continue
        antisym = fo.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        n = fo.count() // 2 if antisym else fo.count()
        segments.append(SegmentInfo(
            kind             = "ORBIT",
            start            = cursor,
            length           = n,
            op_u             = fo.operation.name,
            flavour_u        = tuple(fo.flavour.counts),
            sign_compressed  = antisym,
        ))
        cursor += n

    # Phase 2: one ASSOC or NULL segment per PairFlavour, grouped into blocks.
    for _bkey, block_iter in groupby(pf_list, key=_block_key):
        block = list(block_iter)
        max_idx = max(range(len(block)), key=lambda i: block[i].count(group_sizes))
        for i, pf in enumerate(block):
            sc = _symmetry_class(pf)
            if i == max_idx:
                segments.append(SegmentInfo(
                    kind           = "NULL",
                    start          = cursor,   # cursor does not advance
                    length         = 0,
                    op_u           = pf.op_u.name,
                    flavour_u      = tuple(pf.flavour_u.counts),
                    op_v           = pf.op_v.name,
                    flavour_v      = tuple(pf.flavour_v.counts),
                    overlap        = tuple(pf.overlap),
                    symmetry_class = sc,
                ))
            else:
                length = 2 * pf.count(group_sizes)   # n complex coeffs = 2n reals
                segments.append(SegmentInfo(
                    kind           = "ASSOC",
                    start          = cursor,
                    length         = length,
                    op_u           = pf.op_u.name,
                    flavour_u      = tuple(pf.flavour_u.counts),
                    op_v           = pf.op_v.name,
                    flavour_v      = tuple(pf.flavour_v.counts),
                    overlap        = tuple(pf.overlap),
                    symmetry_class = sc,
                ))
                cursor += length

    return segments
