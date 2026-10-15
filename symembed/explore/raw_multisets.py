"""Raw multisets (relabelling only): orbit construction in symembed versus embedding the n points directly."""
import sys, time
from dataclasses import dataclass
import numpy as np
sys.path.insert(0, ".")                  # run from the repository root
from symembed import EmbeddingOptions, GroupSpec
from symembed.build import _build_from_table
from symembed.bursar import bursar_size
from symembed.slots import iter_group_elements


@dataclass(frozen=True)
class RawSlot:
    kind: str
    positions: tuple
    def name(self, labels=()): return f"x_{self.positions[0]}{self.positions[1]}"
    def latex(self, labels=()): return self.name()
    def L_degree(self, k=1): return 1


def raw(n, d, emb):
    spec = GroupSpec((tuple(f"o{i}" for i in range(n)),))
    slots = tuple(RawSlot(f"f{f}", (i, f)) for f in range(d) for i in range(n))      # orbit order: by feature
    where = {(s.positions): t for t, s in enumerate(slots)}
    els = tuple(iter_group_elements(spec))
    index = [[where[(g.perm[s.positions[0]], s.positions[1])] for s in slots] for g in els]
    sign = [[1] * len(slots) for _ in els]
    _, meta = _build_from_table(spec, d, EmbeddingOptions(embedder=emb, reductions="none", max_outputs=10**9),
                                slots, els, index, sign, False, spec.species_index_blocks, vector_w=False)
    return meta.D


print(f"{'n':>2} {'d':>2} | {'orbit: graded':>14} {'nested':>8} {'cascade':>8} | {'direct Bursar':>13} {'Simplex 2':>9}")
for n, d in [(2, 2), (2, 4), (3, 2), (3, 4), (4, 4), (5, 4)]:
    row = []
    for emb in ("bursar_graded", "bursar_nested", "bursar_cascade"):
        t = time.time()
        row.append(raw(n, d, emb))
    print(f"{n:>2} {d:>2} | {row[0]:>14} {row[1]:>8} {row[2]:>8} | {bursar_size(n, d):>13} {2 * n * d + 1 - d:>9}", flush=True)
