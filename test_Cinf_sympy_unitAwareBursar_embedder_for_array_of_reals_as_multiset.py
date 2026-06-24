#!/usr/bin/env python
#
# Tests for the UnitAwareBursarialEncoder (writeup: DOCS/dimensional_multiset_embedding).
# Property-based rather than golden-value: we check the invariants the construction
# is supposed to have, on worked inputs.

import numpy as np
from Cinf_sympy_unitAwareBursar_embedder_for_array_of_reals_as_multiset import (
    UnitAwareBursarialEncoder, Embedder)

# (mass, length, length) over base order (M, L): k=3, B=2, two entries share unit L.
DIMS = [(1, 0), (0, 1), (0, 1)]


def _rng(seed=0):
    return np.random.default_rng(seed)


def _rand_data(rng, n, k):
    return rng.uniform(-5, 5, size=(n, k)).astype(np.float64)


def test_alias_is_class():
    assert Embedder is UnitAwareBursarialEncoder


def test_size_matches_embedding_length():
    enc = UnitAwareBursarialEncoder(DIMS)
    rng = _rng()
    for n in (2, 3, 4):
        emb, shape, _meta = enc.embed(_rand_data(rng, n, 3))
        assert shape == (n, 3)
        assert len(emb) == enc.size_from_n_k(n, 3) > 0


def test_size_minus_one_for_wrong_k():
    enc = UnitAwareBursarialEncoder(DIMS)              # promised k = 3
    assert enc.size_from_n_k(3, 2) == -1
    assert enc.size_from_n_k(3, 5) == -1
    assert enc.size_from_n_k(1, 2) == -1               # n==1 edge still guarded
    assert enc.size_from_n_k(3, 3) > 0


def test_permutation_invariance():
    enc = UnitAwareBursarialEncoder(DIMS)
    rng = _rng(1)
    data = _rand_data(rng, 4, 3)
    emb0, *_ = enc.embed(data)
    emb1, *_ = enc.embed(data[rng.permutation(4)])
    assert np.allclose(emb0, emb1, atol=1e-8)


def test_injectivity_necessary_checks():
    enc = UnitAwareBursarialEncoder(DIMS)
    rng = _rng(2)
    A = _rand_data(rng, 3, 3)
    B = _rand_data(rng, 3, 3)
    embA, *_ = enc.embed(A)
    embB, *_ = enc.embed(B)
    assert not np.allclose(embA, embB, atol=1e-6)                       # distinct -> distinct
    embA_relabelled, *_ = enc.embed(A[rng.permutation(3)])
    assert np.allclose(embA, embA_relabelled, atol=1e-8)               # relabel -> same
    # a repeated vector (a genuine multiset repeat) must still embed cleanly:
    Arep = np.array([[1., 2., 3.], [1., 2., 3.], [4., 5., 6.]])
    embrep, *_ = enc.embed(Arep)
    assert len(embrep) == enc.size_from_n_k(3, 3)


def test_dimensional_covariance():
    # Change base unit b by lambda: entry c's value scales by lambda^{d_{c,b}}, and
    # each coefficient must rescale by lambda^{(n-a) D_b - m_b}.
    enc = UnitAwareBursarialEncoder(DIMS)
    rng = _rng(3)
    lam, b, n = 2.0, 1, 3                              # b=1 is the length axis
    data = _rand_data(rng, n, 3)
    emb, *_ = enc.embed(data)

    data_scaled = data.copy()
    for c in range(3):
        data_scaled[:, c] *= lam ** DIMS[c][b]
    emb_scaled, *_ = enc.embed(data_scaled)

    support = enc._canonical_support(n)                # monomials (w, t, xi0, xi1)
    Db = enc._D[b]
    for idx, mono in enumerate(support):
        a, m_b = mono[0], mono[2 + b]
        if abs(emb[idx]) > 1e-9:
            predicted = lam ** ((n - a) * Db - m_b)
            assert np.isclose(emb_scaled[idx] / emb[idx], predicted, rtol=1e-6)


def test_negative_exponent_derived_units():
    # velocity (L/T) and energy (M L^2 / T^2) over base order (M, L, T).
    enc = UnitAwareBursarialEncoder([(0, 1, -1), (1, 2, -2)])
    rng = _rng(4)
    data = _rand_data(rng, 3, 2)
    emb, *_ = enc.embed(data)
    assert len(emb) == enc.size_from_n_k(3, 2)
    emb_perm, *_ = enc.embed(data[rng.permutation(3)])
    assert np.allclose(emb, emb_perm, atol=1e-8)


def test_uniform_units_match_unit_free_bursar_size():
    # Three lengths (rho = 0): cost should collapse to the dimensionless Bursar size.
    enc = UnitAwareBursarialEncoder([(1,), (1,), (1,)])
    for n in (2, 3, 4):
        expected = n + (3 - 1) * n * (n + 1) // 2      # n + (k-1) n(n+1)/2 with k=3
        assert enc.size_from_n_k(n, 3) == expected


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"[PASS] {name}")
    print("ALL TESTS PASSED")
