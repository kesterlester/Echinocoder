"""
Outputs of definite parity (option ``parity_split``).

Parity p reflects every vector, which fixes every dot product d and negates every determinant e. An output o(L) that
mixes both behaviours is replaced by its even and odd parts,

    o_even(L) = (o(L) + o(pL)) / 2,        o_odd(L) = (o(L) - o(pL)) / 2,

which together determine o (their sum), so injectivity is unaffected. An output is kept as it is when it is provably
pure: an exact counting argument shows that every monomial that can occur in it has the same number of e factors
modulo 2. When parity is in G every output is already even, and nothing is split.
"""
from __future__ import annotations

import json

import numpy as np
import torch
from torch import nn

from .bursar import bursar_output_labels


def possible_e_parities(position_kinds: list, N: int) -> list[set]:
    """
    For each Bursarial coefficient E_{c,j} (c = 1..N, j = 0..c(W-1), in output order): the parities of the number of
    e factors in the monomials that can occur. ``position_kinds[j]`` is the set of kinds ("d", "e") involved in the
    coordinate at power x^j.
    """
    W = len(position_kinds)
    L = N * (W - 1) + 1
    steps = [{int(kind == "e") for kind in kinds} for kinds in position_kinds]
    reach = np.zeros((L, 2), dtype=bool)
    for j, ps in enumerate(steps):
        for p in ps:
            reach[j, p] = True
    out = []
    for c in range(1, N + 1):
        if c > 1:
            new = np.zeros_like(reach)
            for j, ps in enumerate(steps):
                for p in ps:
                    new[j:, :] |= np.roll(reach[:L - j, :], p, axis=1) if p else reach[:L - j, :]
            reach = new
        for j in range(c * (W - 1) + 1):
            out.append({int(p) for p in np.nonzero(reach[j])[0]})
    return out


def graded_parities(packing, slots, N: int) -> list[set]:
    """
    Possible e-parities of each graded coefficient E_{c,j} (in the order of packing.surviving(N)), by a refinement
    of the counting argument. E_{c,j} is a sum over patterns J (multisets of c occupied positions summing to j) of
    pieces F_J. Parity p acts on the coordinate at position a either as a sign (``S`` mode: only d's, sign +1, or
    only e's, sign -1), or, when s is used, as a sign times complex conjugation (``C`` mode: phi_a(pv) =
    +-conj(phi_a(v))). If every coordinate of J has S mode, F_J has the product of their signs as its parity. If
    every coordinate of J has C mode, so does F_J: since phi_a(s v) = conj(phi_a(v)) for every packed coordinate,
    replacing w by pw conjugates every entry, which is the same as relabelling g -> s g. Patterns mixing the modes
    (or using a coordinate with neither) are undetermined, and so is any E_{c,j} that has one.
    """
    import numpy as np
    from .cyclotomic import Cyclotomic
    F = Cyclotomic(packing.L)
    flip = np.array([1 if slot.kind == "d" else -1 for slot in slots], dtype=object)
    modes = {}
    for a, (U, weight) in enumerate(zip(packing.exact, packing.weights)):
        if weight is None:
            continue
        kinds = {slots[t].kind for t in range(len(slots)) if np.any(U[t])}
        s_mode = 0 if kinds == {"d"} else (1 if kinds == {"e"} else None)
        c_mode = None
        if packing.real:
            UP, Uc = U * flip[:, None], U @ F.conj_matrix.T
            c_mode = 0 if np.array_equal(UP, Uc) else (1 if np.array_equal(UP, -Uc) else None)
        modes[a] = (s_mode, c_mode)
    add = lambda x, y: None if x is None or y is None else (x + y) % 2
    reach = {0: {(0, 0)}}
    by_c = {}
    for c in range(1, N + 1):
        new: dict[int, set] = {}
        for j, states in reach.items():
            for a, (ms, mc) in modes.items():
                bucket = new.setdefault(j + a, set())
                for ss, sc in states:
                    bucket.add((add(ss, ms), add(sc, mc)))
        reach = new
        by_c[c] = reach
    out = []
    for c, j in packing.surviving(N):
        parities = set()
        for ss, sc in by_c[c].get(j, set()):
            if ss is None and sc is None:
                parities = {0, 1}
                break
            if ss is None or sc is None or ss == sc:
                parities.add(ss if sc is None else sc)
        out.append(parities)
    return out


def output_parities(meta) -> list[set]:
    """The possible e-parities of each output of a plain, involution or graded embedding (exact, combinatorial)."""
    from .graded import GradedMeta
    from .rootnorm import RootNormaliseMeta
    if isinstance(meta, RootNormaliseMeta):          # sign(f)|f|^(1/r) is odd in f: parities are unchanged
        return output_parities(meta.inner)
    if isinstance(meta, GradedMeta):
        p = meta.packing
        kinds = [{meta.slots[t].kind for t in range(meta.m) if abs(a[t]) > 1e-12 or abs(b[t]) > 1e-12}
                 for a, b in zip(p.plus, p.minus)]
        invariants = [{int(meta.slots[t].kind == "e") for t in range(meta.m) if abs(v[t]) > 1e-12}
                      for v in p.invariants]
        if p.W == 0:
            full = invariants
        else:
            by_cj = dict(zip(p.surviving(meta.N), graded_parities(p, meta.slots, meta.N)))
            full = invariants + [by_cj[(c, j)] for c, j, _, _ in meta.selection()]
        return [full[i] for i in meta.kept]
    if hasattr(meta, "base_output_labels") and hasattr(meta, "position_kinds"):
        by_label = dict(zip(meta.base_output_labels, possible_e_parities(meta.position_kinds(), meta.N)))
        return [by_label[label] for label in meta.output_labels]
    raise NotImplementedError(f"parity_split is not available for {type(meta).__name__}")


class ParitySplitEmbedding(nn.Module):
    """Wraps an embedding: provably pure outputs are kept; every other output is replaced by its even and odd parts."""

    def __init__(self, inner: nn.Module, parities: list, k: int):
        super().__init__()
        self.inner = inner
        D = len(parities)
        plan = []
        for i, ps in enumerate(parities):
            if len(ps) <= 1:
                plan.append(i)                         # pure (or identically zero): as it is
            else:
                plan.extend([D + i, 2 * D + i])        # even part, odd part
        self.register_buffer("plan", torch.tensor(plan, dtype=torch.int64))
        reflection = torch.ones(k, dtype=inner.g_sign.dtype)
        reflection[-1] = -1.0
        self.register_buffer("reflection", reflection)
        self.register_buffer("input_shape", inner.input_shape.clone())
        self.register_buffer("output_shape", torch.tensor([len(plan)], dtype=torch.int64))

    def forward(self, x) -> torch.Tensor:
        x = self.inner._prepare(x)
        o, o_reflected = self.inner(x), self.inner(x * self.reflection)
        stacked = torch.cat([o, (o + o_reflected) / 2, (o - o_reflected) / 2], dim=-1)
        return stacked[..., self.plan]


class ParitySplitMeta:
    """Metadata for a parity-split embedding: the inner embedding's metadata, plus which outputs were split."""

    def __init__(self, inner_meta, parities: list):
        self.inner, self.parities = inner_meta, parities
        self.relations = ()

    def __getattr__(self, name):                    # delegate everything else to the inner metadata
        return getattr(self.inner, name)

    @property
    def certification_level(self) -> int:
        return self.inner.certification_level

    def _inner_names(self) -> list[str]:
        names = getattr(self.inner, "output_names", None)       # also found through wrapper delegation
        if callable(names):
            return names()
        return [self.inner.base_name(b) for b in self.inner.kept_base_indices]     # EmbeddingMeta

    def output_names(self) -> list[str]:
        names, out = self._inner_names(), []
        for name, ps in zip(names, self.parities):
            if len(ps) == 0:
                out.append(f"{name} [zero]")
            elif len(ps) == 1:
                out.append(f"{name} [{'even' if 0 in ps else 'odd'}]")
            else:
                out.extend([f"even part of {name}", f"odd part of {name}"])
        return out

    @property
    def D(self) -> int:
        return len(self.output_names())

    @property
    def certification(self) -> str:
        return (self.inner.certification + " Parity split: each output not provably pure is replaced by its even and "
                "odd parts, (o(L) + o(pL))/2 and (o(L) - o(pL))/2, which together determine it.")

    def to_dict(self) -> dict:
        d = self.inner.to_dict()
        d["D"] = self.D
        d["certification"] = self.certification
        d["parity_split"] = {"split": sum(1 for ps in self.parities if len(ps) > 1),
                             "pure_even": sum(1 for ps in self.parities if ps == {0}),
                             "pure_odd": sum(1 for ps in self.parities if ps == {1})}
        d["outputs"] = [{"index": i, "name": nm} for i, nm in enumerate(self.output_names())]
        return d

    def to_json(self, **kwargs) -> str:
        return json.dumps(self.to_dict(), **kwargs)
