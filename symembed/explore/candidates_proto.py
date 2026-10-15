"""Scratch prototype of candidates(problem): enumerate methods, build them, report size and degree statistics."""
import sys, time, itertools
from dataclasses import dataclass
import numpy as np
import torch
sys.path.insert(0, ".")                  # run from the repository root
from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding
from symembed.beam import BeamSpec, build_beam_embedding
from symembed.bursar import bursar_output_labels, bursar_size
from symembed.build import _build_from_table
from symembed.slots import GroupElement, iter_group_elements


def degrees(meta) -> list[int]:
    """Degree in w of each output (upper bound: the largest degree of any monomial that can occur)."""
    inner = getattr(meta, "inner", None)
    m = inner if (inner is not None and type(meta).__name__ == "BeamMeta") else meta
    name = type(m).__name__
    if name == "EmbeddingMeta":
        return [c for c, _ in m.base_output_labels] if not m.relations else \
               [m.base_output_labels[b][0] for b in m.kept_base_indices]
    if name == "GradedMeta":
        sel = m.selection()
        n_inv = len(m.packing.invariants)
        full = [1] * n_inv + [c for c, _, _, _ in sel]
        return [full[i] for i in m.kept]
    if name == "NestedMeta":
        inner_deg = max(degrees(m.inner))
        return [c * inner_deg for c, _ in bursar_output_labels(len(m.reps), m.inner.D)]
    if name == "CascadeMeta":
        p = m.plan
        deg = [1] * len(p.base_re)
        for st in p.stages:
            dmax = max([deg[i] for i in st.packed], default=0)
            new = [deg[i] for i in st.passed]
            labels = bursar_output_labels(st.p, st.W) if st.W else []
            new += [labels[q][0] * dmax for q in st.select]
            deg = new
        out = [1] * len(p.invariants)
        if p.s is None:
            out += deg + deg
        else:
            out += [d for d, f in zip(deg, p.final_imaginary) if not f] + [2 * max(deg, default=0)] * max(0, 2 * len(deg) - 1)
        if len(m.reps) == 1:
            return out
        return [c * max(out) for c, _ in bursar_output_labels(len(m.reps), p.count)]
    if name == "BeamMeta":        # projective baseline: plain Bursar over perms of (even, -q^2): degree 2c at most
        return [2 * c for c, _ in bursar_output_labels(meta.N, meta._width)]
    raise TypeError(name)


def report(label, build):
    t = time.time()
    try:
        _, meta = build()
    except Exception as e:
        print(f"  {label:44s} FAILED: {str(e)[:70]}")
        return
    d = degrees(meta)
    print(f"  {label:44s} D={meta.D:6d}  max deg {max(d):3d}  mean deg {np.mean(d):5.1f}  build {time.time()-t:5.1f}s",
          flush=True)


print("LHC: 3 jets (masses vary) + 2 leptons (common mass), axial + beam exchange (beam frame)")
spec = BeamSpec((("j1", "j2", "j3"), ("l1", "l2")), masses=(None, 0.0))
for emb in ("bursar_graded", "bursar_nested", "bursar_cascade"):
    report(f"beam frame, swap in group, {emb}", lambda: build_beam_embedding(spec, EmbeddingOptions(embedder=emb, reductions="safe" if emb == "bursar_graded" else "none")))
report("beam frame, swap first, bursar_cascade", lambda: build_beam_embedding(spec, EmbeddingOptions(embedder="bursar_cascade"), swap="projective"))
report("beam frame, plain RP baseline", lambda: build_beam_embedding(spec, swap="projective"))
report("beam frame, swap in group, plain bursar", lambda: build_beam_embedding(spec, EmbeddingOptions(max_outputs=10**8)))

print("\nSame particles, no beams: Lorentz, S(3) x S(2) on 5 vectors (k=4)")
g = GroupSpec((("j1", "j2", "j3"), ("l1", "l2")))
for emb in ("bursar_graded", "bursar_nested", "bursar_cascade"):
    report(f"dots + dets, {emb}", lambda: build_invariant_embedding(4, g, EmbeddingOptions(embedder=emb, reductions="safe" if emb == "bursar_graded" else "none")))
