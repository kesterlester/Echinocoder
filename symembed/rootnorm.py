"""
Root normalisation (option ``root_normalise``; design objective 3 in the notes).

An output f that is provably homogeneous of degree r in the components of L (an exact counting argument: every
monomial that can occur in it has the same degree) is replaced by sign(f)|f|^(1/r). That is a bijection of R, so
injectivity is unaffected, and the result scales linearly with L, which tames the dynamic range (useful in float32).
Outputs that are not provably homogeneous are left unchanged. The derivative of sign(f)|f|^(1/r) is infinite at f = 0;
there it is defined to be finite (|f| is clamped to the smallest positive normal number of the dtype).
"""
from __future__ import annotations

import json

import torch
from torch import nn


class _SignedRoot(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, exponent):
        ctx.save_for_backward(x, exponent)
        return torch.sign(x) * x.abs() ** exponent

    @staticmethod
    def backward(ctx, grad):
        x, exponent = ctx.saved_tensors
        tiny = torch.finfo(x.dtype).tiny
        return grad * exponent * x.abs().clamp_min(tiny) ** (exponent - 1), None


def output_degrees(meta) -> list[list[int]]:
    """Possible degrees in L of each output (exact, combinatorial), for plain, involution and graded embeddings."""
    from .graded import GradedMeta
    from .meta import possible_L_degrees
    if isinstance(meta, GradedMeta):
        return meta.possible_degrees()
    if hasattr(meta, "base_output_labels") and hasattr(meta, "position_kinds"):
        degrees = possible_L_degrees(meta.slots, meta.k, meta.N, meta.position_kinds())
        by_label = dict(zip(meta.base_output_labels, degrees))
        return [by_label[label] for label in meta.output_labels]
    raise NotImplementedError(f"root_normalise is not available for {type(meta).__name__}")


def numerically_zero(inner: nn.Module, degrees: list, n: int, k: int, samples: int = 16, seed: int = 0) -> list[bool]:
    """
    Outputs that vanish (to rounding) at random inputs, relative to the other outputs of the same degree. Such outputs
    are typically identically zero as functions of L (for example through the Gram relations), and their rounding
    noise would be amplified enormously by a root, so they are left unnormalised.
    """
    g = torch.Generator().manual_seed(seed)
    X = torch.randn(samples, n, k, dtype=torch.float64, generator=g)
    with torch.no_grad():
        Y = inner.double()(X).abs().max(dim=0).values
        inner.to(inner.g_sign.dtype)
    key = [tuple(d) for d in degrees]
    block = {}
    for i, kk in enumerate(key):
        block[kk] = max(block.get(kk, 0.0), Y[i].item())
    return [Y[i].item() <= 1e-9 * block[kk] for i, kk in enumerate(key)]


class RootNormaliseEmbedding(nn.Module):
    """Wraps an embedding, replacing each provably homogeneous output f of degree r by sign(f)|f|^(1/r)."""

    def __init__(self, inner: nn.Module, degrees: list, skip: list | None = None):
        super().__init__()
        self.inner = inner
        skip = skip or [False] * len(degrees)
        exponents = [1.0 / d[0] if len(d) == 1 and d[0] > 0 and not z else 1.0 for d, z in zip(degrees, skip)]
        self.register_buffer("exponent", torch.tensor(exponents, dtype=inner.g_sign.dtype))
        self.register_buffer("input_shape", inner.input_shape.clone())
        self.register_buffer("output_shape", inner.output_shape.clone())
        self.g_sign = inner.g_sign                         # dtype/device anchor, as for other modules
        self._prepare = inner._prepare

    def forward(self, x) -> torch.Tensor:
        return _SignedRoot.apply(self.inner(x), self.exponent)


class RootNormaliseMeta:
    """Metadata for a root-normalised embedding: the inner metadata, plus which outputs were normalised."""

    def __init__(self, inner_meta, degrees: list, skip: list | None = None):
        self.inner, self.degrees = inner_meta, degrees
        self.skip = skip or [False] * len(degrees)

    def __getattr__(self, name):
        return getattr(self.inner, name)

    @property
    def normalised(self) -> list[int]:
        return [i for i, (d, z) in enumerate(zip(self.degrees, self.skip)) if len(d) == 1 and d[0] > 0 and not z]

    @property
    def certification(self) -> str:
        return (self.inner.certification + f" Root normalisation: {len(self.normalised)} provably homogeneous "
                "outputs f of degree r replaced by sign(f)|f|^(1/r), a bijection of R, so injectivity is unaffected.")

    def to_dict(self) -> dict:
        d = self.inner.to_dict()
        d["certification"] = self.certification
        d["root_normalise"] = {"normalised": len(self.normalised), "unchanged": len(self.degrees) - len(self.normalised),
                               "numerically_zero_left_unnormalised": [i for i, z in enumerate(self.skip) if z],
                               "degrees": [dd[0] if len(dd) == 1 else None for dd in self.degrees]}
        return d

    def to_json(self, **kwargs) -> str:
        return json.dumps(self.to_dict(), **kwargs)
