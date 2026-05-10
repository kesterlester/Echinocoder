from __future__ import annotations
from dataclasses import dataclass
from symatom import ArgumentSymmetry
from symatom.rep import canonical_pair_flavours
from symatom import repL


def _symmetry_class(pf) -> str:
    u = "A" if pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC else "S"
    v = "A" if pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC else "S"
    return f"{u}{v}"


@dataclass(frozen=True)
class SegmentInfo:
    """
    Metadata for one contiguous block of complex coefficients in the encode()
    output vector.

    Fields
    ------
    start          : first index (inclusive) in the output array
    length         : number of complex coefficients in this segment
    op_u, op_v     : names of the two operations forming this PairFlavour
    flavour_u      : counts tuple — how many labels from each group go into op_u
    flavour_v      : counts tuple — how many labels from each group go into op_v
    overlap        : how many labels are shared between op_u and op_v per group
    symmetry_class : two-character code — "SS", "SA", "AS", or "AA"
                     (S = SYMMETRIC, A = ANTISYMMETRIC for op_u then op_v)

    Human-readable display
    ----------------------
    str(seg)  — single-line summary, suitable for printing a table.

    Machine-readable export
    -----------------------
    seg.to_dict()  — plain dict; json.dumps-able for downstream tools.

    Notes
    -----
    Segment lengths are computed analytically from pf.count(group_sizes) and
    do not require evaluating the encoding on any event.  A "dry-run with
    all-zero inputs" would give the same lengths but is unnecessary here.
    """
    start:          int
    length:         int
    op_u:           str
    op_v:           str
    flavour_u:      tuple
    flavour_v:      tuple
    overlap:        tuple
    symmetry_class: str

    @property
    def stop(self) -> int:
        return self.start + self.length

    def __str__(self) -> str:
        fl_u = ",".join(str(c) for c in self.flavour_u)
        fl_v = ",".join(str(c) for c in self.flavour_v)
        ov   = ",".join(str(c) for c in self.overlap)
        return (
            f"[{self.start}:{self.stop}]"
            f"  {self.op_u}×{self.op_v}"
            f"  {self.symmetry_class}"
            f"  u=({fl_u})  v=({fl_v})  shared=({ov})"
            f"  len={self.length}"
        )

    def to_dict(self) -> dict:
        return {
            "start":          self.start,
            "stop":           self.stop,
            "length":         self.length,
            "op_u":           self.op_u,
            "op_v":           self.op_v,
            "flavour_u":      list(self.flavour_u),
            "flavour_v":      list(self.flavour_v),
            "overlap":        list(self.overlap),
            "symmetry_class": self.symmetry_class,
        }


def describe_encoding(plan) -> list[SegmentInfo]:
    """
    Return a list of SegmentInfo objects describing the structure of the vector
    produced by encode(plan, event).

    The i-th SegmentInfo covers output[seg.start : seg.stop] and describes
    which PairFlavour generated that block, how many coefficients it contains,
    and what symmetry compression was applied.

    Usage examples
    --------------
    # Human-readable table:
    for seg in describe_encoding(plan):
        print(seg)

    # Machine-readable JSON:
    import json
    json.dumps([s.to_dict() for s in describe_encoding(plan)])

    The segments are listed in the same order as the coefficients in the
    encode() output, so output[seg.start : seg.stop] always corresponds to
    the PairFlavour described by that SegmentInfo.

    This function is pure — it does not evaluate any event data.
    """
    fo_list = repL(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    group_sizes = tuple(g.size for g in plan.context.groups)

    segments = []
    cursor = 0
    for pf in pf_list:
        length = pf.count(group_sizes)
        segments.append(SegmentInfo(
            start          = cursor,
            length         = length,
            op_u           = pf.op_u.name,
            op_v           = pf.op_v.name,
            flavour_u      = tuple(pf.flavour_u.counts),
            flavour_v      = tuple(pf.flavour_v.counts),
            overlap        = tuple(pf.overlap),
            symmetry_class = _symmetry_class(pf),
        ))
        cursor += length

    return segments
