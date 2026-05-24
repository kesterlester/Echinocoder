from __future__ import annotations
from dataclasses import dataclass
from itertools import groupby
from symatom import Operation, repS
from symatom.rep import canonical_pair_flavours
from .pairs import _is_self_pair


def _orbit_example(op_name: str, flavour_counts: tuple, types) -> str:
    """Canonical label example for an ORBIT, e.g. 'dot(a,p)'."""
    labels = []
    for vt, k in zip(types, flavour_counts):
        labels.extend(list(vt.labels)[:k])
    return f"{op_name}({','.join(str(l) for l in labels)})"


def _assoc_example(op_u: str, flavour_u: tuple, op_v: str, flavour_v: tuple,
                   overlap: tuple, types) -> str:
    """Canonical label example for an ASSOC/NULL, e.g. 'dot(p,q)×dot(a,p)'."""
    u_labels, v_labels = [], []
    for vt, ku, kv, s in zip(types, flavour_u, flavour_v, overlap):
        glabels = list(vt.labels)
        u_g = glabels[:ku]
        shared = u_g[:s]
        remaining = [l for l in glabels if l not in u_g]
        v_g = shared + remaining[:kv - s]
        u_labels.extend(u_g)
        v_labels.extend(v_g)
    u_str = f"{op_u}({','.join(str(l) for l in u_labels)})"
    v_str = f"{op_v}({','.join(str(l) for l in v_labels)})"
    return f"{u_str}×{v_str}"



def _block_key(pf):
    """Key that identifies an OVERLAP BLOCK: fixed (op_u, flavour_u, op_v, flavour_v).

    Uses the Operation objects directly — they are frozen dataclasses whose
    equality and hash include eval_fn identity, so two operations with the same
    name but different eval_fn produce distinct keys.
    """
    return (pf.op_u, pf.flavour_u.counts, pf.op_v, pf.flavour_v.counts)


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
                atom stored (fo.count_of_atoms_one_per_sign() values).

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
                      ORBIT: fo.count_of_atoms_one_per_sign() // 2 if sign_compressed else fo.count_of_atoms_one_per_sign().
                      ASSOC: two reals per association pair (2 * pf.count()),
                             since each complex polynomial coefficient = (re, im).
    op_u            : operation u (or the single operation for ORBIT) — full Operation object
    flavour_u       : counts tuple — labels from each group going into op_u
    op_v            : operation v (None for ORBIT) — full Operation object
    flavour_v       : counts tuple for op_v (None for ORBIT)
    overlap         : shared-label counts per group (None for ORBIT)
    symmetry_class  : orbit type label — "11", "12", "21", "22", "neg",
                      or "null_self" (None for ORBIT)
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
    op_u:             Operation
    flavour_u:        tuple
    op_v:             Operation | None = None
    flavour_v:        tuple | None = None
    overlap:          tuple | None = None
    symmetry_class:   str   | None = None
    sign_compressed:  bool  | None = None
    notional_length:  int   | None = None  # length before drop; None → same as length
    method_name:      str   | None = None  # encoder method that produced this segment
    example:          str   | None = None

    @property
    def stop(self) -> int:
        return self.start + self.length

    def __str__(self) -> str:
        idx     = f"[{self.start}:{self.stop}]"
        fl_u    = ",".join(str(c) for c in self.flavour_u)
        op_v    = self.op_v.name if self.op_v is not None else "."
        fl_v    = (",".join(str(c) for c in self.flavour_v)
                   if self.flavour_v is not None else ".")
        ov      = (",".join(str(c) for c in self.overlap)
                   if self.overlap   is not None else ".")
        # Single "variant" column: symmetry class (SS/SA/AS/AA) for pair rows,
        # SC (sign-compressed) or "." for ORBIT rows.
        if self.kind == "ORBIT":
            if self.method_name is not None:
                variant = self.method_name
            else:
                variant = "SC" if self.sign_compressed else "."
        else:
            variant = self.symmetry_class
        full    = self.notional_length if self.notional_length is not None else self.length
        ex      = f"  |  {self.example}" if self.example is not None else ""
        return (
            f"{idx}  {self.kind}  {self.op_u.name}  {op_v}  {variant}"
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
            "op_u":             self.op_u.name,
            "flavour_u":        list(self.flavour_u),
        }
        if self.kind == "ORBIT":
            d["sign_compressed"] = self.sign_compressed
            d["method_name"]     = self.method_name
        else:
            d.update({
                "op_v":           self.op_v.name,
                "flavour_v":      list(self.flavour_v),
                "overlap":        list(self.overlap),
                "symmetry_class": self.symmetry_class,
            })
        return d


# ---------------------------------------------------------------------------
# Tree node types for hierarchical describe output
# ---------------------------------------------------------------------------

@dataclass
class Phase1Tree:
    """Flat container for all Phase 1 (orbit) segments."""
    orbits: list  # list[SegmentInfo]

    def flat(self) -> list:
        return list(self.orbits)

    def __iter__(self):
        return iter(self.orbits)

    def __len__(self):
        return len(self.orbits)

    def __getitem__(self, idx):
        return self.orbits[idx]


@dataclass
class OverlapBlockNode:
    """One overlap block in the Phase 2 tree: fixed (op_u, flavour_u, op_v, flavour_v)."""
    op_u:      Operation
    flavour_u: tuple
    op_v:      Operation
    flavour_v: tuple
    segments:  list  # list[SegmentInfo]

    def flat(self) -> list:
        return list(self.segments)

    def __str__(self) -> str:
        fl_u = ",".join(str(c) for c in self.flavour_u)
        fl_v = ",".join(str(c) for c in self.flavour_v)
        return f"OverlapBlock  {self.op_u.name} × {self.op_v.name}  u=({fl_u})  v=({fl_v})"


@dataclass
class Phase2Tree:
    """Container for all Phase 2 overlap block nodes."""
    blocks: list  # list[OverlapBlockNode]

    def flat(self) -> list:
        return [seg for block in self.blocks for seg in block.flat()]


@dataclass
class EncodingTree:
    """
    Top-level tree mirroring the full encoder hierarchy.

    phase1  — flat collection of orbit segments (Phase 1)
    phase2  — nested collection of overlap-block nodes (Phase 2)

    List-like protocol (flat iteration) is supported so that existing code
    using `for s in tree`, `tree[i]`, `len(tree)` continues to work unchanged.
    Use tree.phase1 / tree.phase2 for structured access.
    """
    phase1: Phase1Tree
    phase2: Phase2Tree

    def flat(self) -> list:
        return self.phase1.flat() + self.phase2.flat()

    def __iter__(self):
        return iter(self.flat())

    def __len__(self):
        return len(self.flat())

    def __getitem__(self, idx):
        return self.flat()[idx]

    def __eq__(self, other):
        if isinstance(other, EncodingTree):
            return self.phase1 == other.phase1 and self.phase2 == other.phase2
        if isinstance(other, list):
            return self.flat() == other
        return NotImplemented


# describe_encoding() has moved to encode.py so that it shares the same
# encoding loop as encode(), deriving all metadata from the registry's
# assess() responses rather than from a separate algebraic computation.
