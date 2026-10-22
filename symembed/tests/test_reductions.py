"""Safe reductions (Level 0 omission of determinants with parity; Level 1 exact linear relations)."""
import json
from fractions import Fraction

import numpy as np
import pytest
import sympy
import torch

from symembed import EmbeddingOptions, GroupSpec, InvariantEmbedding, build_invariant_embedding
from symembed.reductions import _rational_reconstruct, exact_expansion, find_linear_relations

SAFE = EmbeddingOptions(reductions="safe", dtype=torch.float64)
EXAMPLE = GroupSpec((("a", "b"), ("p",)))
EXAMPLE_P = GroupSpec((("a", "b"), ("p",)), parity=True)
TWO_PAIRS = GroupSpec((("a", "b"), ("p", "q")))
TRIPLE_P = GroupSpec((("a", "b", "c"),), parity=True)


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(4321)


def unreduced_module(meta):
    """The Bursarial embedding on the same w and multiset, keeping every coefficient."""
    return InvariantEmbedding(meta.n, meta.k, meta.slots, meta.index, meta.sign, torch.float64)


# ----- the exact expansion ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("k,spec", [(2, EXAMPLE), (3, TWO_PAIRS), (2, GroupSpec((("a", "b", "c"),)))])
def test_exact_expansion_matches_sympy(k, spec):
    _, meta = build_invariant_embedding(k, spec, EmbeddingOptions(dtype=torch.float64))
    w = meta.slot_symbols()
    expanded = meta.expanded_base_outputs()
    c_max = 2
    position = 0
    for c, (monos, coeffs) in enumerate(exact_expansion(meta.index, meta.sign, meta.m, c_max), start=1):
        for j in range(coeffs.shape[0]):
            ours = sum(int(coeffs[j, t]) * sympy.Mul(*[w[s] for s in mono]) for t, mono in enumerate(monos))
            assert sympy.expand(ours - expanded[position]) == 0
            position += 1


# ----- exact linear algebra --------------------------------------------------------------------------------------

def test_find_linear_relations_synthetic():
    r0 = np.array([1, 2, 0, 3])
    r1 = np.array([0, 1, 5, -1])
    rows = np.array([r0, r1, 2 * r0, 0 * r0, r0 - r1, 3 * r0 + 7 * r1, [1, 0, 0, 0], 7 * r0 - 3 * r1])
    kept, relations = find_linear_relations(rows)
    assert kept == [0, 1, 6]
    rel = {r.dropped: dict(r.terms) for r in relations}
    assert rel == {2: {0: 2}, 3: {}, 4: {0: 1, 1: -1}, 5: {0: 3, 1: 7}, 7: {0: 7, 1: -3}}


def test_find_linear_relations_fractions():
    r0 = np.array([2, 4, 6, 1])
    r1 = np.array([3, 0, 9, 0])
    third = np.array([1, 0, 3, 0])                # = r1 / 3
    rows = np.array([r0, r1, third, 7 * third - 2 * r0])
    kept, relations = find_linear_relations(rows)
    assert kept == [0, 1]
    rel = {r.dropped: dict(r.terms) for r in relations}
    assert rel[2] == {1: Fraction(1, 3)}
    assert rel[3] == {0: -2, 1: Fraction(7, 3)}


def test_unreconstructable_relation_is_kept():
    # A dependency whose rational coefficient is too large to reconstruct modulo the prime must be kept: unsafe
    # drops are impossible, only missed compression.
    r0 = np.array([100003, 0, 200006])
    r1 = np.array([0, 1, 0])
    huge = np.array([100019, 0, 200038])          # = (100019/100003) * r0
    kept, relations = find_linear_relations(np.vstack([r0, r1, huge]))
    assert kept == [0, 1, 2] and relations == []


def test_rational_reconstruct():
    p = 2_147_483_629
    for frac in (Fraction(1), Fraction(-1), Fraction(1, 2), Fraction(-7, 3), Fraction(123, 4567)):
        residue = frac.numerator * pow(frac.denominator, p - 2, p) % p
        assert _rational_reconstruct(residue) == frac


# ----- the reduced embeddings ------------------------------------------------------------------------------------

def test_worked_example_relations():
    _, meta = build_invariant_embedding(2, EXAMPLE, SAFE)
    assert [meta.relation_summary(r) for r in meta.relations] == [
        "E_{1,1} = E_{1,0}", "E_{1,5} = E_{1,4}", "E_{1,6} = 0", "E_{1,8} = E_{1,7}"]
    assert meta.D == 26 - 4 and meta.certification_level == 1


def test_parity_omits_determinants():
    embed, meta = build_invariant_embedding(2, EXAMPLE_P, SAFE)
    assert meta.e_dropped
    assert all(s.kind == "d" for s in meta.slots) and meta.m == 6
    assert meta.N == 2 and not any(g.parity for g in meta.elements)
    assert embed.e_subsets.shape[0] == 0
    assert "Determinants omitted" in meta.certification


@pytest.mark.parametrize("k,spec", [(2, EXAMPLE), (2, EXAMPLE_P), (3, TWO_PAIRS), (3, TRIPLE_P),
                                    (2, GroupSpec((("a", "b", "c"), ("p", "q"))))])
def test_dropped_outputs_reconstruct_exactly(k, spec, gen):
    embed, meta = build_invariant_embedding(k, spec, SAFE)
    X = torch.randn(20, spec.n, k, dtype=torch.float64, generator=gen)
    full = unreduced_module(meta)(X)
    torch.testing.assert_close(meta.reconstruct_base(embed(X)), full, rtol=1e-12, atol=1e-12 * full.abs().max())
    torch.testing.assert_close(embed(X), full[:, meta.kept_base_indices])


@pytest.mark.parametrize("k,spec", [(2, EXAMPLE), (3, TWO_PAIRS), (3, TRIPLE_P)])
def test_relations_hold_symbolically(k, spec):
    _, meta = build_invariant_embedding(k, spec, SAFE)
    base = meta.expanded_base_outputs()
    for rel in meta.relations:
        rhs = sum((sympy.Rational(f.numerator, f.denominator) * base[t] for t, f in rel.terms), sympy.Integer(0))
        assert sympy.expand(base[rel.dropped] - rhs) == 0
        c = meta.base_output_labels[rel.dropped][0]
        assert all(meta.base_output_labels[t][0] == c and t in meta.kept_base_indices for t, _ in rel.terms)


@pytest.mark.parametrize("parity", [False, True])
def test_reduced_embedding_invariance(parity, gen):
    spec = GroupSpec((("a", "b", "c"), ("p", "q")), parity=parity)
    embed, _ = build_invariant_embedding(3, spec, SAFE)
    X = torch.randn(spec.n, 3, dtype=torch.float64, generator=gen)
    Y = embed(X)
    scale = Y.abs().max()
    torch.testing.assert_close(embed(X[[2, 0, 1, 4, 3]]), Y, rtol=1e-9, atol=1e-9 * scale)
    q, _ = torch.linalg.qr(torch.randn(3, 3, dtype=torch.float64, generator=gen))   # det +1 or -1
    Y_turned = embed(X @ q.T)
    if parity or torch.linalg.det(q) > 0:
        torch.testing.assert_close(Y_turned, Y, rtol=1e-9, atol=1e-9 * scale)


def test_c_max_zero_disables_linear_reduction():
    _, meta = build_invariant_embedding(2, EXAMPLE_P, EmbeddingOptions(reductions="safe", linear_reduction_max_c=0))
    assert meta.relations == () and meta.e_dropped and meta.certification_level == 0


def test_meta_records_reductions():
    _, meta = build_invariant_embedding(2, EXAMPLE, SAFE)
    d = json.loads(meta.to_json())
    assert d["D"] == 22 and d["base_D"] == 26 and d["certification_level"] == 1
    linear = [r for r in d["reductions"] if r["type"] == "linear"]
    assert [r["summary"] for r in linear][0] == "E_{1,1} = E_{1,0}"
    assert all(r["verified"].startswith("exact") for r in linear)
    assert [o["name"] for o in d["outputs"][:5]] == ["E_{1,0}", "E_{1,2}", "E_{1,3}", "E_{1,4}", "E_{1,7}"]
    text = meta.latex()
    assert "exact identity" in text and "E_{1,6}" in text
