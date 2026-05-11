from __future__ import annotations
from dataclasses import dataclass
from itertools import groupby
from symatom.atoms import ArgumentSymmetry
from symatom.rep import canonical_pair_flavours
from symatom import repL
from .pairs import _is_self_pair


def _orbit_example(op_name: str, flavour_counts: tuple, groups) -> str:
    """Canonical label example for an ORBIT, e.g. 'dot(a,p)'."""
    labels = []
    for group, k in zip(groups, flavour_counts):
        labels.extend(list(group.labels)[:k])
    return f"{op_name}({','.join(str(l) for l in labels)})"


def _assoc_example(op_u: str, flavour_u: tuple, op_v: str, flavour_v: tuple,
                   overlap: tuple, groups) -> str:
    """Canonical label example for an ASSOC/NULL, e.g. 'dot(p,q)×dot(a,p)'."""
    u_labels, v_labels = [], []
    for group, ku, kv, s in zip(groups, flavour_u, flavour_v, overlap):
        glabels = list(group.labels)
        u_g = glabels[:ku]
        shared = u_g[:s]
        remaining = [l for l in glabels if l not in u_g]
        v_g = shared + remaining[:kv - s]
        u_labels.extend(u_g)
        v_labels.extend(v_g)
    u_str = f"{op_u}({','.join(str(l) for l in u_labels)})"
    v_str = f"{op_v}({','.join(str(l) for l in v_labels)})"
    return f"{u_str}×{v_str}"


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
    "ORBIT"     Phase 1: sort-encoded single-operator orbit.
                Fields op_v, flavour_v, overlap, symmetry_class are None.
                sign_compressed=True: ANTISYMMETRIC operation — the orbit contains
                {+u, -u} pairs; only |eval(sign=+1 atom)| is stored per combination
                (n values, n = base label combinations).
                sign_compressed=False: SYMMETRIC/UNSTRUCTURED — eval(u) for every
                atom stored (fo.count() values).

    "ASSOC"     Phase 2: compressed pair encoding of one association.
                One PairFlavour within an OVERLAP BLOCK, encoded via _embed_compressed.
                All fields are set.

    "NULL_SELF" Phase 2: self-pairing association dropped because its z-values
                are entirely determined by Phase 1.  Condition: op_u==op_v,
                flavour_u==flavour_v, overlap==flavour_u (maximum overlap), so
                z_k = (1+i)*a_k where a_k comes from the Phase 1 orbit.
                length == 0; all fields are set.  Dropped before complementation.

    "NULL_COMP" Phase 2: largest remaining (non-self-pair) association per OVERLAP
                BLOCK, dropped because it is deducible by multiset complementation
                from the stored associations and Phase 1 orbit headers.
                length == 0; all fields are set.

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
    example         : decorative canonical label string, e.g. 'dot(a,p)' for
                      ORBIT or 'dot(p,q)×dot(a,p)' for ASSOC/NULL.  Set by
                      describe_encoding() from the plan's group labels; None if
                      not computed.  Not needed for decoding.

    Human-readable display
    ----------------------
    str(seg)  — fixed 9-token format (suitable for pipe through `column -t`):
        [start:stop]  kind  op_u  op_v  variant  u=(...)  v=(...)  shared=(...)  len=N
    where `variant` is the symmetry class (SS/SA/AS/AA) for pair rows, or
    SC (sign-compressed) / "." for ORBIT rows.  Unused fields show ".".
    An optional "  |  example" tail follows when example is set.

    Machine-readable export
    -----------------------
    seg.to_dict()  — plain dict; json.dumps-able for downstream tools.

    Notes
    -----
    NULL_SELF and NULL_COMP entries appear at their logical position in the
    OVERLAP BLOCK sequence but contribute zero coefficients; start == stop.
    Within each block: NULL_SELF entries (if any) appear before NULL_COMP (if any),
    because the drop order is self-pair first, complementation last.
    Consecutive entries within the same block are ordered by
    canonical_pair_flavours' overlap sort (lex ascending).
    """
    kind:             str
    start:            int
    length:           int
    op_u:             str
    flavour_u:        tuple
    op_v:             str   | None = None
    flavour_v:        tuple | None = None
    overlap:          tuple | None = None
    symmetry_class:   str   | None = None
    sign_compressed:  bool  | None = None
    notional_length:  int   | None = None  # length before drop; None → same as length
    example:          str   | None = None

    @property
    def stop(self) -> int:
        return self.start + self.length

    def __str__(self) -> str:
        idx     = f"[{self.start}:{self.stop}]"
        fl_u    = ",".join(str(c) for c in self.flavour_u)
        op_v    = self.op_v    if self.op_v    is not None else "."
        fl_v    = (",".join(str(c) for c in self.flavour_v)
                   if self.flavour_v is not None else ".")
        ov      = (",".join(str(c) for c in self.overlap)
                   if self.overlap   is not None else ".")
        # Single "variant" column: symmetry class (SS/SA/AS/AA) for pair rows,
        # SC (sign-compressed) or "." for ORBIT rows.
        if self.kind == "ORBIT":
            variant = "SC" if self.sign_compressed else "."
        else:
            variant = self.symmetry_class
        full    = self.notional_length if self.notional_length is not None else self.length
        ex      = f"  |  {self.example}" if self.example is not None else ""
        return (
            f"{idx}  {self.kind}  {self.op_u}  {op_v}  {variant}"
            f"  u=({fl_u})  v=({fl_v})  shared=({ov})"
            f"  len={self.length}  full={full}{ex}"
        )

    def to_dict(self) -> dict:
        full = self.notional_length if self.notional_length is not None else self.length
        d = {
            "kind":             self.kind,
            "start":            self.start,
            "stop":             self.stop,
            "length":           self.length,
            "notional_length":  full,
            "op_u":             self.op_u,
            "flavour_u":        list(self.flavour_u),
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
    groups = plan.context.groups

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
            notional_length  = n,
            op_u             = fo.operation.name,
            flavour_u        = tuple(fo.flavour.counts),
            sign_compressed  = antisym,
            example          = _orbit_example(fo.operation.name, fo.flavour.counts, groups),
        ))
        cursor += n

    # Phase 2: one ASSOC, NULL_SELF, or NULL_COMP segment per PairFlavour.
    # Drop order: self-pairs first (NULL_SELF, Phase-1 redundant), then the
    # largest remaining non-self-pair per block (NULL_COMP, complementation).
    for _bkey, block_iter in groupby(pf_list, key=_block_key):
        block = list(block_iter)
        non_self_idx = [i for i, pf in enumerate(block) if not _is_self_pair(pf)]
        comp_drop = (max(non_self_idx, key=lambda i: block[i].count(group_sizes))
                     if non_self_idx else None)
        for i, pf in enumerate(block):
            sc = _symmetry_class(pf)
            ex = _assoc_example(
                pf.op_u.name, pf.flavour_u.counts,
                pf.op_v.name, pf.flavour_v.counts,
                pf.overlap, groups,
            )
            notional = 2 * pf.count(group_sizes)
            common = dict(
                op_u            = pf.op_u.name,
                flavour_u       = tuple(pf.flavour_u.counts),
                op_v            = pf.op_v.name,
                flavour_v       = tuple(pf.flavour_v.counts),
                overlap         = tuple(pf.overlap),
                symmetry_class  = sc,
                notional_length = notional,
                example         = ex,
            )
            if _is_self_pair(pf):
                segments.append(SegmentInfo(
                    kind   = "NULL_SELF",
                    start  = cursor,
                    length = 0,
                    **common,
                ))
            elif i == comp_drop:
                segments.append(SegmentInfo(
                    kind   = "NULL_COMP",
                    start  = cursor,
                    length = 0,
                    **common,
                ))
            else:
                segments.append(SegmentInfo(
                    kind   = "ASSOC",
                    start  = cursor,
                    length = notional,
                    **common,
                ))
                cursor += notional

    return segments
