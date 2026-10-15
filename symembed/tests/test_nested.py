"""Nested coset packing (embedder="bursar_nested")."""
import json

import pytest
import torch

from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding
from symembed.bursar import bursar_coefficients
from symembed.graded import GradedMeta
from symembed.nested import NestedMeta

NESTED = EmbeddingOptions(embedder="bursar_nested", dtype=torch.float64)
# A small case where nesting wins: parity in G with the determinants kept (reductions="none").
SMALL_WIN = (4, GroupSpec((("a", "b", "c"), ("p", "q")), parity=True))


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(97531)


def test_nested_used_when_it_wins():
    k, spec = SMALL_WIN
    _, meta = build_invariant_embedding(k, spec, NESTED)
    flat = build_invariant_embedding(k, spec, EmbeddingOptions(embedder="bursar_graded"))[1].D
    assert isinstance(meta, NestedMeta) and meta.D == 371 < flat == 460


def test_falls_back_to_graded_when_nesting_loses():
    _, meta = build_invariant_embedding(2, GroupSpec((("a", "b", "c"),)), NESTED)
    assert isinstance(meta, GradedMeta) and meta.D == 32


def test_outer_bursar_of_inner_outputs(gen):
    k, spec = SMALL_WIN
    embed, meta = build_invariant_embedding(k, spec, NESTED)
    X = torch.randn(2, spec.n, k, dtype=torch.float64, generator=gen)
    w = embed.w(X)
    inner = torch.stack([embed.inner.from_w(w[:, embed.g_index[t]] * embed.g_sign[t])
                         for t in range(embed.n_cosets)], dim=1)
    torch.testing.assert_close(embed(X), bursar_coefficients(inner))
    assert inner.shape[1] == len(meta.reps) == meta.N // len(meta.H)


def test_invariance_including_parity(gen):
    k, spec = SMALL_WIN
    embed, _ = build_invariant_embedding(k, spec, NESTED)
    X = torch.randn(spec.n, k, dtype=torch.float64, generator=gen)
    Y = embed(X)
    scale = Y.abs().max()
    torch.testing.assert_close(embed(X[[2, 0, 1, 4, 3]]), Y, rtol=1e-9, atol=1e-9 * scale)
    q, _ = torch.linalg.qr(torch.randn(k, k, dtype=torch.float64, generator=gen))    # det +1 or -1: parity in G
    torch.testing.assert_close(embed(X @ q.T), Y, rtol=1e-9, atol=1e-9 * scale)


def test_separation_smoke(gen):
    k, spec = SMALL_WIN
    embed, _ = build_invariant_embedding(k, spec, NESTED)
    X = torch.randn(20, spec.n, k, dtype=torch.float64, generator=gen)
    Y = embed(X)
    assert (torch.cdist(Y, Y) + torch.eye(20, dtype=torch.float64)).min() > 1e-8


def test_gradients(gen):
    # gradcheck's elementwise tolerances do not suit outputs of such high degree (Jacobian entries span ~1e46 down
    # to rounding noise), so compare the analytic Jacobian with central differences relative to its largest entry.
    k, spec = SMALL_WIN
    embed, _ = build_invariant_embedding(k, spec, NESTED)
    X = torch.randn(1, spec.n, k, dtype=torch.float64, generator=gen)
    J = torch.autograd.functional.jacobian(embed, X).reshape(-1, spec.n * k)
    eps = 1e-6
    numeric = torch.zeros_like(J)
    with torch.no_grad():
        for i in range(spec.n * k):
            d = torch.zeros(spec.n * k, dtype=torch.float64)
            d[i] = eps
            d = d.reshape(1, spec.n, k)
            numeric[:, i] = ((embed(X + d) - embed(X - d)) / (2 * eps)).reshape(-1)
    assert (J - numeric).abs().max() < 1e-6 * J.abs().max()


def test_meta():
    k, spec = SMALL_WIN
    _, meta = build_invariant_embedding(k, spec, NESTED)
    d = json.loads(meta.to_json())
    assert d["D"] == 371 and d["nested"]["order_H"] * len(d["nested"]["coset_representatives"]) == meta.N
    assert d["certification"].startswith("Level 0") and "Nested coset packing" in meta.latex()


def test_s5_count():
    _, meta = build_invariant_embedding(4, GroupSpec((tuple("abcde"),)), EmbeddingOptions(embedder="bursar_nested"))
    assert isinstance(meta, NestedMeta) and meta.D == 6775
