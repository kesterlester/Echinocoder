#!/usr/bin/env python
# Verifies that the fast embedders produce the same pre-hash intermediate state
# as the original embedders.  The final numerical embeddings differ because the
# two implementations use different hash functions (seeded numpy RNG vs
# sequential MD5 per dimension), but the inputs to the hash step must be
# identical for correctness of permutation-invariance and injectivity.

import numpy as np
import data_sources
import C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset as s1_orig
import C0HomDeg1_simplicialComplex_embedder_1_fast_for_array_of_reals_as_multiset as s1_fast
import C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset as s2_orig
import C0HomDeg1_simplicialComplex_embedder_2_fast_for_array_of_reals_as_multiset as s2_fast

embedder_s1_orig = s1_orig.Embedder()
embedder_s1_fast = s1_fast.Embedder()
embedder_s2_orig = s2_orig.Embedder()
embedder_s2_fast = s2_fast.Embedder()


def _check_pre_hash_state(data, orig_embedder, fast_embedder, label):
    sd_orig, cc_orig = orig_embedder._pre_hash_state(data)
    sd_fast, cc_fast = fast_embedder._pre_hash_state(data)
    np.testing.assert_array_almost_equal(sd_orig, sd_fast, decimal=10,
        err_msg=f"{label}: second_diffs mismatch for data shape {data.shape}")
    np.testing.assert_array_equal(cc_orig, cc_fast,
        err_msg=f"{label}: canonical_eji_counts mismatch for data shape {data.shape}")


_SHAPES = [(2, 3), (3, 3), (4, 3), (3, 2), (5, 4)]


def test_simplex1_pre_hash_state_matches():
    rng = np.random.default_rng(42)
    for n, k in _SHAPES:
        data = rng.standard_normal((n, k))
        _check_pre_hash_state(data, embedder_s1_orig, embedder_s1_fast, "S1")


def test_simplex2_pre_hash_state_matches():
    rng = np.random.default_rng(42)
    for n, k in _SHAPES:
        data = rng.standard_normal((n, k))
        _check_pre_hash_state(data, embedder_s2_orig, embedder_s2_fast, "S2")


def test_simplex1_pre_hash_state_matches_fixed():
    data = np.asarray([[4, 2, 3],
                       [-3, 5, 1],
                       [8, 9, 2],
                       [2, 7, 2]], dtype=float)
    _check_pre_hash_state(data, embedder_s1_orig, embedder_s1_fast, "S1 fixed")


def test_simplex2_pre_hash_state_matches_fixed():
    data = np.asarray([[4, 2, 3],
                       [-3, 5, 1],
                       [8, 9, 2],
                       [2, 7, 2]], dtype=float)
    _check_pre_hash_state(data, embedder_s2_orig, embedder_s2_fast, "S2 fixed")


if __name__ == "__main__":
    test_simplex1_pre_hash_state_matches_fixed()
    test_simplex2_pre_hash_state_matches_fixed()
    test_simplex1_pre_hash_state_matches()
    test_simplex2_pre_hash_state_matches()
    print("All pre-hash state comparison tests passed.")
