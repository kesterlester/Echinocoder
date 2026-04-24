"""
Fast vectorised re-implementation of Simplex-1 embedder.

Mathematically identical pipeline to
C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset.py
up to and including the canonical-form step.  The hash that maps each
canonical Eji_LinComb to a point in the unit hypercube is replaced with a
seeded numpy RNG (one MD5 seed per vertex instead of one MD5 call per
dimension per vertex), so numerical output values differ from the original
embedder — but permutation-invariance and injectivity are preserved.

Key differences vs the original:
  1. MSV indicator matrices are built with np.cumsum instead of Python sets.
  2. Eji_LinComb._eji_counts prefix sums use np.cumsum instead of iterative
     add() calls — reducing O((nk)^3) Python iterations to one O((nk)^2)
     numpy op.
  3. hash_to_point_in_unit_hypercube: one MD5 seed → np.random.default_rng
     instead of bigN sequential MD5 updates per vertex.
  4. Final weighted sum is a matrix multiply (second_diffs @ points).
"""

import hashlib
import numpy as np
from MultisetEmbedder import MultisetEmbedder
from typing import Any


class Embedder(MultisetEmbedder):

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        return MultisetEmbedder.embed_kOne_sorting(data), None

    def _pre_hash_state(self, data: np.ndarray):
        """Return (second_diffs, canonical_eji_counts) — the inputs to the hash step.

        second_diffs       : float64 array, shape (num_v,)
        canonical_eji_counts: uint16 array,  shape (num_v, n, k)

        These are identical to the corresponding values in the original
        implementation and can be used to verify correctness independently of
        the choice of hash function.
        """
        assert MultisetEmbedder.is_generic_data(data)
        n, k = data.shape
        nk = n * k
        num_v = nk - 1

        values = data.ravel().astype(np.float64)
        jj, ii = np.indices((n, k))
        perm1 = np.argsort(-values, kind='stable')
        sorted_values = values[perm1]
        sorted_js = jj.ravel()[perm1]
        sorted_is = ii.ravel()[perm1]

        first_diffs = sorted_values[:-1] - sorted_values[1:]

        single_inds = np.zeros((nk, n, k), dtype=np.uint16)
        single_inds[np.arange(nk), sorted_js, sorted_is] = 1
        msv_indicators = np.cumsum(single_inds, axis=0)[:num_v]

        perm2 = np.argsort(-first_diffs, kind='stable')
        sorted_deltas = first_diffs[perm2]
        sorted_msv = msv_indicators[perm2]

        next_d = np.empty(num_v, dtype=np.float64)
        next_d[:-1] = sorted_deltas[1:]
        next_d[-1] = 0.0
        scale = np.arange(1, num_v + 1, dtype=np.float64)
        second_diffs = scale * (sorted_deltas - next_d)

        eji_counts = np.cumsum(sorted_msv, axis=0)

        canonical = np.empty_like(eji_counts)
        for p in range(num_v):
            m = eji_counts[p]
            canonical[p] = m[np.lexsort(m.T[::-1])]

        return second_diffs, canonical

    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        assert MultisetEmbedder.is_generic_data(data)
        n, k = data.shape
        nk = n * k
        num_v = nk - 1
        bigN = 2 * num_v + 1

        values = data.ravel().astype(np.float64)
        perm1 = np.argsort(-values, kind='stable')
        sorted_values = values[perm1]
        min_val = sorted_values[-1]
        max_val = sorted_values[0]

        second_diffs, canonical = self._pre_hash_state(data)

        # ------------------------------------------------------------------ #
        # Hash each canonical form to a bigN-dimensional point.               #
        # One MD5 seed + one numpy RNG call replaces the original's bigN     #
        # sequential MD5 updates per vertex.                                  #
        # ------------------------------------------------------------------ #
        points = np.empty((num_v, bigN), dtype=np.float64)
        for p in range(num_v):
            h = hashlib.md5()
            h.update(canonical[p].tobytes())
            h.update(np.array([p + 1], dtype=np.uint16).tobytes())
            seed = int.from_bytes(h.digest()[:4], 'little')
            points[p] = np.random.default_rng(seed).random(bigN)

        core = second_diffs @ points

        embedding = np.empty(bigN + 2, dtype=np.float64)
        embedding[:bigN] = core
        embedding[-2] = max_val
        embedding[-1] = min_val
        return embedding, None

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        return 2 * n * k + 1
