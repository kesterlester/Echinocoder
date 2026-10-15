"""Tests for the cascaded graded embedding (embedder="bursar_cascade")."""
import json

import pytest
import torch

from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding
from symembed.cascade import factorisations, plan_stage
from symembed.slots import iter_group_elements

S3S2 = GroupSpec((("a", "b", "c"), ("p", "q")))
CASCADE = EmbeddingOptions(embedder="bursar_cascade", dtype=torch.float64)


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(1928)


@pytest.mark.parametrize("k,spec,reductions,D", [
    (4, S3S2, "none", 63), (3, S3S2, "none", 95), (2, S3S2, "none", 109),
    (4, GroupSpec(S3S2.species, parity=True), "safe", 52), (3, GroupSpec((("a", "b", "c"),)), "none", 12),
    (4, GroupSpec((("a", "b", "c", "d"),)), "none", 153),          # nested: 3 cosets of a D_4
])
def test_counts_invariance_and_local_injectivity(k, spec, reductions, D, gen):
    embed, meta = build_invariant_embedding(k, spec, EmbeddingOptions(embedder="bursar_cascade", reductions=reductions,
                                                                      dtype=torch.float64))
    assert meta.D == D
    X = torch.randn(3, spec.n, k, dtype=torch.float64, generator=gen)
    Y = embed(X)
    for el in iter_group_elements(spec):
        Xg = X[:, list(el.perm)].clone()
        if el.parity:
            Xg[..., -1] = -Xg[..., -1]
        torch.testing.assert_close(embed(Xg), Y, rtol=1e-9, atol=1e-9 * Y.abs().max())
    Q, _ = torch.linalg.qr(torch.randn(k, k, dtype=torch.float64, generator=gen))
    if torch.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]                                          # a rotation, not a reflection
    torch.testing.assert_close(embed(X @ Q.T), Y, rtol=1e-9, atol=1e-9 * Y.abs().max())
    # Generic evidence of injectivity up to G (the proof is the certificate): full rank in w.
    w = torch.randn(1, meta.m, dtype=torch.float64, generator=gen)
    J = torch.autograd.functional.jacobian(embed.from_w, w).reshape(-1, meta.m)
    J = J / J.norm(dim=1, keepdim=True).clamp_min(1e-300)
    assert torch.linalg.matrix_rank(J, rtol=1e-9).item() == meta.m


def test_not_invariant_beyond_g(gen):
    embed, _ = build_invariant_embedding(4, S3S2, CASCADE)
    X = torch.randn(2, 5, 4, dtype=torch.float64, generator=gen)
    Y = embed(X)
    assert not torch.allclose(embed(X[:, [3, 1, 2, 0, 4]]), Y)      # a and p are different species
    Xr = X.clone()
    Xr[..., -1] = -Xr[..., -1]
    assert not torch.allclose(embed(Xr), Y)                          # parity is not in G


def test_stage_rule():
    # A stage passes Z_p-fixed coordinates on, and keeps coefficients with j = shift c (mod p).
    stage = plan_stage([0, 1, 2, 3, 4, 5, 3, 2], [False] * 8, 6, 2)
    assert stage.passed == [0, 2, 4, 7]
    assert all((j - stage.shift * 0) % 6 == (w + stage.shift) % 6 for j, w in zip(stage.position, stage.weight))
    assert sorted(factorisations(6)) == [[2, 3], [3, 2], [6]]


def test_meta():
    _, meta = build_invariant_embedding(4, S3S2, CASCADE)
    d = json.loads(meta.to_json())
    assert d["D"] == 63 and d["cascade"]["order_D"] == 12 and d["cascade"]["chain"] in ([2, 3], [3, 2])
    assert d["certification"].startswith("Level 0") and len(meta.output_names()) == 63
