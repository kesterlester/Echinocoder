"""Involution packing (embedder="bursar_involution"): packing algebra, real coefficients, invariance, reductions."""
import json

import numpy as np
import pytest
import sympy
import torch

from symembed import EmbeddingOptions, GroupSpec, InvariantEmbedding, build_invariant_embedding
from symembed.bursar import bursar_coefficients_complex
from symembed.involution import (choose_involution, is_involution, packing_from_matrix,
                                 packing_from_signed_involution, trace)

INV = EmbeddingOptions(embedder="bursar_involution", dtype=torch.float64)
INV_SAFE = EmbeddingOptions(embedder="bursar_involution", reductions="safe", dtype=torch.float64)
S2 = GroupSpec((("a", "b"),))
S3 = GroupSpec((("a", "b", "c"),))
S2S1 = GroupSpec((("a", "b"), ("p",)))
S3S2 = GroupSpec((("a", "b", "c"), ("p", "q")))


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(2468)


def matrix(row_index, row_sign):
    m = len(row_index)
    T = np.zeros((m, m), dtype=int)
    for s in range(m):
        T[s, row_index[s]] = row_sign[s]
    return T


def check_packing(packing, T):
    P, Q = np.array(packing.plus, dtype=float), np.array(packing.minus, dtype=float)
    # g0 acts as conjugation: real parts invariant, imaginary parts negated.
    assert np.array_equal(P @ T, P) and np.array_equal(Q @ T, -Q)
    # The packing is a bijection from R^m to the packed coordinates (zero padding rows aside).
    assert np.linalg.matrix_rank(np.vstack([P, Q])) == T.shape[0]


# ----- packing algebra -------------------------------------------------------------------------------------------

@pytest.mark.parametrize("k,spec", [(2, S2), (2, S3), (2, S2S1), (3, S3S2), (2, GroupSpec((("a", "b"),), True))])
def test_signed_packing_for_every_involution(k, spec):
    _, meta = build_invariant_embedding(k, spec)
    count = 0
    for row_index, row_sign in zip(meta.index, meta.sign):
        if is_involution(row_index, row_sign):
            packing = packing_from_signed_involution(row_index, row_sign)
            check_packing(packing, matrix(row_index, row_sign))
            assert packing.r == (meta.m + abs(trace(row_index, row_sign))) // 2
            count += 1
    assert count >= 2      # the identity and at least one transposition


def test_generic_packing_for_a_non_permutation_involution():
    rng = np.random.default_rng(0)
    for _ in range(5):
        while True:
            S = rng.integers(-3, 4, size=(5, 5))
            if round(abs(np.linalg.det(S))) == 1:
                break
        D = np.diag([1, 1, -1, 1, -1])
        T = np.rint(S @ D @ np.linalg.inv(S)).astype(int)
        assert np.array_equal(T @ T, np.eye(5, dtype=int))
        packing = packing_from_matrix(T)
        assert packing.r == 3
        check_packing(packing, T)


def test_generic_and_signed_packings_span_the_same_spaces():
    _, meta = build_invariant_embedding(3, S3S2)
    g = choose_involution(meta.index, meta.sign)
    T = matrix(meta.index[g], meta.sign[g])
    a, b = packing_from_signed_involution(meta.index[g], meta.sign[g]), packing_from_matrix(T)
    for x, y in ((a.plus, b.plus), (a.minus, b.minus)):
        x, y = np.array(x, dtype=float), np.array(y, dtype=float)
        assert np.linalg.matrix_rank(x) == np.linalg.matrix_rank(y) == np.linalg.matrix_rank(np.vstack([x, y]))


def test_choose_involution_minimises_trace():
    _, meta = build_invariant_embedding(2, S2)
    g = choose_involution(meta.index, meta.sign)
    assert meta.elements[g].cycle_notation(S2.labels) == "(a b)" and trace(meta.index[g], meta.sign[g]) == 0


# ----- the embedding ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("k,spec,D", [(2, S2, 5), (2, S2S1, 14), (2, S3, 90), (3, S3S2, 948)])
def test_output_counts(k, spec, D):
    _, meta = build_invariant_embedding(k, spec, INV)
    assert meta.D == D == len(meta.base_output_labels)


@pytest.mark.parametrize("k,spec", [(2, S2), (2, S2S1), (1, GroupSpec((("a", "b", "c"),)))])
def test_imaginary_parts_vanish_symbolically_and_numerically(k, spec, gen):
    embed, meta = build_invariant_embedding(k, spec, INV)
    expanded = meta.expanded_base_outputs()          # raises if any imaginary part is non-zero
    X = torch.randn(spec.n, k, dtype=torch.float64, generator=gen)
    w = embed.w(X)
    gw = (w[embed.g_index] * embed.g_sign).unsqueeze(0)
    re, im = bursar_coefficients_complex(gw @ embed.pack_plus.T, gw @ embed.pack_minus.T)
    assert im.abs().max() <= 1e-12 * re.abs().max()
    values = dict(zip(meta.slot_symbols(), w.tolist()))
    expected = torch.tensor([float(sympy.sympify(e).subs(values)) for e in expanded], dtype=torch.float64)
    torch.testing.assert_close(embed(X), expected, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("parity", [False, True])
def test_invariance(parity, gen):
    spec = GroupSpec((("a", "b", "c"), ("p", "q")), parity=parity)
    embed, _ = build_invariant_embedding(3, spec, INV)
    X = torch.randn(spec.n, 3, dtype=torch.float64, generator=gen)
    Y = embed(X)
    scale = Y.abs().max()
    torch.testing.assert_close(embed(X[[1, 2, 0, 4, 3]]), Y, rtol=1e-9, atol=1e-9 * scale)
    q, _ = torch.linalg.qr(torch.randn(3, 3, dtype=torch.float64, generator=gen))
    Y_turned = embed(X @ q.T)
    if parity or torch.linalg.det(q) > 0:
        torch.testing.assert_close(Y_turned, Y, rtol=1e-9, atol=1e-9 * scale)
    else:
        assert (Y_turned - Y).abs().max() > 1e-6 * scale


def test_separates_like_the_unpacked_embedding(gen):
    # Smoke test only: on random inputs and on group-related inputs, packed and unpacked outputs agree about
    # whether two inputs are equivalent.
    packed, _ = build_invariant_embedding(2, S2S1, INV)
    plain, _ = build_invariant_embedding(2, S2S1, EmbeddingOptions(dtype=torch.float64))
    X = torch.randn(30, 3, 2, dtype=torch.float64, generator=gen)
    Xg = X[:, [1, 0, 2]]
    for f in (packed, plain):
        assert (f(X) - f(Xg)).abs().max() < 1e-9 * f(X).abs().max()
        d = torch.cdist(f(X), f(X)) + torch.eye(30, dtype=torch.float64)
        assert d.min() > 1e-8


def test_gradients(gen):
    embed, _ = build_invariant_embedding(2, S2S1, INV)
    X = torch.randn(2, 3, 2, dtype=torch.float64, generator=gen, requires_grad=True)
    assert torch.autograd.gradcheck(embed, (X,))


def test_no_involution_is_an_error():
    with pytest.raises(ValueError, match="involution"):
        build_invariant_embedding(2, GroupSpec((("a",), ("b",))), INV)


# ----- with the safe reductions ----------------------------------------------------------------------------------

@pytest.mark.parametrize("k,spec", [(3, S3S2), (2, GroupSpec((("a", "b"), ("p",)), parity=True)), (2, S3)])
def test_safe_reductions_on_packed_outputs(k, spec, gen):
    embed, meta = build_invariant_embedding(k, spec, INV_SAFE)
    full = InvariantEmbedding(meta.n, k, meta.slots, meta.index, meta.sign, torch.float64, packing=meta.packing)
    X = torch.randn(10, spec.n, k, dtype=torch.float64, generator=gen)
    reference = full(X)
    torch.testing.assert_close(meta.reconstruct_base(embed(X)), reference, rtol=1e-10,
                               atol=1e-10 * reference.abs().max())


def test_relations_on_packed_outputs_hold_symbolically():
    _, meta = build_invariant_embedding(3, GroupSpec((("a", "b"), ("p",))), INV_SAFE)
    base = meta.expanded_base_outputs()
    for rel in meta.relations:
        rhs = sum((sympy.Rational(f.numerator, f.denominator) * base[t] for t, f in rel.terms), sympy.Integer(0))
        assert sympy.expand(base[rel.dropped] - rhs) == 0


def test_parity_with_omitted_determinants_uses_a_permutation_involution():
    _, meta = build_invariant_embedding(2, GroupSpec((("a", "b"), ("p",)), parity=True), INV_SAFE)
    assert meta.e_dropped and not meta.elements[meta.involution].parity


def test_meta_describes_packing():
    _, meta = build_invariant_embedding(2, S2, INV)
    d = json.loads(meta.to_json())
    assert d["packing"]["involution"] == "(a b)" and d["packing"]["trace"] == 0
    assert d["packing"]["coordinates"] == [{"re": "d_aa + d_bb", "im": "d_aa - d_bb"}, {"re": "d_ab", "im": "e_ab"}]
    assert "Involution packing (Level 0)" in d["certification"]
    assert "Involution packing with" in meta.latex(expanded=True)
