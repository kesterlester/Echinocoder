"""Graded packing (embedder="bursar_graded"): properties of the packing, vanishing, invariance, interface."""
import json
from math import comb

import numpy as np
import pytest
import torch

from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding
from symembed.bursar import bursar_coefficients_complex, bursar_output_labels
from symembed.graded import build_packing, check_packing, matrices

GRADED = EmbeddingOptions(embedder="bursar_graded", dtype=torch.float64)
S2 = GroupSpec((("a", "b"),))
S3 = GroupSpec((("a", "b", "c"),))
S2S1 = GroupSpec((("a", "b"), ("p",)))
S3S2 = GroupSpec((("a", "b", "c"), ("p", "q")))
S4 = GroupSpec((("a", "b", "c", "d"),))


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(1357)


@pytest.mark.parametrize("k,spec,D", [(2, S2, 4), (2, S2S1, 11), (2, S3, 32), (3, S3S2, 203), (2, S4, 718),
                                      (3, S4, 469)])
def test_output_counts(k, spec, D):
    _, meta = build_invariant_embedding(k, spec, GRADED)
    assert meta.D == D


@pytest.mark.parametrize("k,spec", [(2, S3), (3, S3S2), (3, S4), (2, S2S1)])
def test_every_candidate_packing_has_the_defining_properties(k, spec):
    _, meta = build_invariant_embedding(k, spec)
    Ms = matrices(meta.index, meta.sign)
    I = np.eye(meta.m)
    involutions = [g for g in range(meta.N) if np.allclose(Ms[g] @ Ms[g], I)]
    checked = 0
    for h in range(0, meta.N, max(1, meta.N // 6)):
        for s in [g for g in involutions if np.allclose(Ms[g] @ Ms[h] @ Ms[g], Ms[h].T)][:3] + [None]:
            for drop in (False, True):
                check_packing(Ms, build_packing(Ms, h, s, drop))
                checked += 1
    assert checked > 10


@pytest.mark.parametrize("k,spec", [(2, S2), (2, S3), (3, S3S2), (2, S4)])
def test_omitted_coefficients_vanish_and_kept_ones_are_real(k, spec, gen):
    embed, meta = build_invariant_embedding(k, spec, GRADED)
    p = meta.packing
    X = torch.randn(spec.n, k, dtype=torch.float64, generator=gen)
    w = embed.w(X)
    gw = (w[embed.g_index] * embed.g_sign).unsqueeze(0)
    Z_re, Z_im = gw @ embed.graded_plus.T, gw @ embed.graded_minus.T
    scale = torch.sqrt(Z_re ** 2 + Z_im ** 2).max()
    re, im = bursar_coefficients_complex(Z_re / scale, Z_im / scale)
    surviving = set(p.surviving(meta.N))
    for i, (c, j) in enumerate(bursar_output_labels(meta.N, p.W)):
        natural = comb(meta.N, c) * p.W ** c
        if (c, j) not in surviving:
            assert abs(re[0, i]) + abs(im[0, i]) < 1e-12 * natural
        elif p.real:
            assert abs(im[0, i]) < 1e-12 * natural


@pytest.mark.parametrize("parity", [False, True])
def test_invariance(parity, gen):
    spec = GroupSpec((("a", "b", "c"), ("p", "q")), parity=parity)
    embed, _ = build_invariant_embedding(3, spec, GRADED)
    X = torch.randn(spec.n, 3, dtype=torch.float64, generator=gen)
    Y = embed(X)
    scale = Y.abs().max()
    torch.testing.assert_close(embed(X[[2, 0, 1, 4, 3]]), Y, rtol=1e-9, atol=1e-9 * scale)
    q, _ = torch.linalg.qr(torch.randn(3, 3, dtype=torch.float64, generator=gen))
    Y_turned = embed(X @ q.T)
    if parity or torch.linalg.det(q) > 0:
        torch.testing.assert_close(Y_turned, Y, rtol=1e-9, atol=1e-9 * scale)
    else:
        assert (Y_turned - Y).abs().max() > 1e-6 * scale


def test_separates_like_the_plain_embedding(gen):
    # Smoke test only: equivalent inputs agree, random inputs are distinguished.
    embed, _ = build_invariant_embedding(2, S3, GRADED)
    X = torch.randn(30, 3, 2, dtype=torch.float64, generator=gen)
    Y = embed(X)
    assert (embed(X[:, [2, 0, 1]]) - Y).abs().max() < 1e-9 * Y.abs().max()
    assert (torch.cdist(Y, Y) + torch.eye(30, dtype=torch.float64)).min() > 1e-8


def test_identity_h_reproduces_involution_packing():
    _, meta = build_invariant_embedding(2, S2S1)
    _, inv_meta = build_invariant_embedding(2, S2S1, EmbeddingOptions(embedder="bursar_involution"))
    Ms = matrices(meta.index, meta.sign)
    packing = build_packing(Ms, 0, inv_meta.involution, False)       # h = identity (element 0)
    check_packing(Ms, packing)
    assert packing.r == 1 and packing.output_count(meta.N) == inv_meta.D == 14


def test_trivial_group_outputs_w_itself():
    _, meta = build_invariant_embedding(2, GroupSpec((("a",), ("b",))), GRADED)
    assert meta.D == meta.m == 4


def test_parity_with_safe_reductions(gen):
    spec = GroupSpec((("a", "b", "c"),), parity=True)
    embed, meta = build_invariant_embedding(2, spec, EmbeddingOptions(embedder="bursar_graded", reductions="safe",
                                                                      dtype=torch.float64))
    assert meta.e_dropped and all(s.kind == "d" for s in meta.slots)
    X = torch.randn(spec.n, 2, dtype=torch.float64, generator=gen)
    X_reflected = X * torch.tensor([1.0, -1.0], dtype=torch.float64)
    torch.testing.assert_close(embed(X_reflected), embed(X))
    assert any("restriction of scalars" in r.get("reason", "") for r in meta.reductions_log())


def test_gradients(gen):
    embed, _ = build_invariant_embedding(2, S3, GRADED)
    X = torch.randn(2, 3, 2, dtype=torch.float64, generator=gen, requires_grad=True)
    assert torch.autograd.gradcheck(embed, (X,))


def test_meta():
    _, meta = build_invariant_embedding(2, S2, GRADED)
    d = json.loads(meta.to_json())
    assert d["D"] == 4 and d["certification_level"] == 0 and d["certification"].startswith("Level 0")
    assert d["graded"]["h"] == "(a b)" and d["graded"]["mode"] == "scalar"
    assert [o["name"] for o in d["outputs"]] == ["inv[0]", "inv[1]", "Re E_{2,0}", "Im E_{2,0}"]
    assert "Graded packing" in meta.latex()


# ----- exact structural reductions (unreachable powers; the c = 1 block) -----------------------------------------

def test_unreachable_powers_are_omitted():
    # With padding (unused positions), some powers that are multiples of r are not sums of c occupied positions.
    _, meta = build_invariant_embedding(3, S3S2, GRADED)
    p = meta.packing
    occupied = [j for j, wt in enumerate(p.weights) if wt is not None]
    assert len(occupied) < p.W
    reach = {0}
    surviving = set(p.surviving(meta.N))
    for c in range(1, meta.N + 1):
        reach = {a + b for a in reach for b in occupied}
        assert all(j in reach for cc, j in surviving if cc == c)


SAFE_GRADED = EmbeddingOptions(embedder="bursar_graded", reductions="safe", dtype=torch.float64)


@pytest.mark.parametrize("k,spec,D", [(3, S3S2, 201), (3, S4, 468), (2, S3, 32)])
def test_c1_reductions_counts(k, spec, D):
    options = EmbeddingOptions(embedder="bursar_graded", reductions="safe", linear_reduction_max_c=1,
                               dtype=torch.float64)
    _, meta = build_invariant_embedding(k, spec, options)
    assert meta.D == D
    assert all(meta.all_names()[rel.dropped].startswith(("Re E_{1,", "Im E_{1,")) for rel in meta.relations)


@pytest.mark.parametrize("k,spec,D", [(3, S3S2, 199), (3, S4, 468), (2, S4, 715)])
def test_exact_cyclotomic_reductions_counts(k, spec, D):
    # c >= 2 relations, found over Q(zeta_L) by restriction of scalars (the build also checks exactly, for every
    # c it reaches, that omitted coefficients vanish and kept ones are real).
    _, meta = build_invariant_embedding(k, spec, SAFE_GRADED)
    assert meta.D == D and meta.linear_reduction_max_c >= 2


def test_coefficients_in_cyclotomic_fields():
    from fractions import Fraction
    from symembed.cyclotomic import Cyclotomic
    from symembed.graded import coefficient_text, coefficient_value
    F = Cyclotomic(12)
    lam = tuple(Fraction(int(x)) for x in F.zeta_power(1) + F.zeta_power(-1))      # 2 cos(pi/6) = sqrt(3)
    assert coefficient_text(lam, 12) == "sqrt(3)"
    assert abs(coefficient_value(lam, 12) - 3 ** 0.5) < 1e-12


@pytest.mark.parametrize("k,spec", [(3, S3S2), (3, S4), (2, GroupSpec((("a", "b", "c"),), parity=True))])
def test_c1_reductions_replay(k, spec, gen):
    embed, meta = build_invariant_embedding(k, spec, SAFE_GRADED)
    full, _ = build_invariant_embedding(k, spec, EmbeddingOptions(embedder="bursar_graded", dtype=torch.float64,
                                                                   reductions="none")) if not spec.parity else (None, None)
    X = torch.randn(10, spec.n, k, dtype=torch.float64, generator=gen)
    Y = embed(X)
    rebuilt = meta.reconstruct_base(Y)
    assert rebuilt.shape[-1] == len(meta.all_names())
    if full is not None:
        reference = full(X)
        torch.testing.assert_close(rebuilt, reference, rtol=1e-10, atol=1e-10 * reference.abs().max())
    assert meta.certification_level == (1 if meta.relations else 0)
