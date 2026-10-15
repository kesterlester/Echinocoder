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


def output_parities(meta) -> list[set]:
    """The possible e-parities of each output of a plain, involution or graded embedding (exact, combinatorial)."""
    from .graded import GradedMeta
    if isinstance(meta, GradedMeta):
        p = meta.packing
        kinds = [{meta.slots[t].kind for t in range(meta.m) if abs(a[t]) > 1e-12 or abs(b[t]) > 1e-12}
                 for a, b in zip(p.plus, p.minus)]
        invariants = [{int(meta.slots[t].kind == "e") for t in range(meta.m) if abs(v[t]) > 1e-12}
                      for v in p.invariants]
        if p.W == 0:
            full = invariants
        else:
            by_label = dict(zip(bursar_output_labels(meta.N, p.W), possible_e_parities(kinds, meta.N)))
            full = invariants + [by_label[(c, j)] for c, j, _, _ in meta.selection()]
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
        if hasattr(type(self.inner), "output_names"):
            return self.inner.output_names()
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
