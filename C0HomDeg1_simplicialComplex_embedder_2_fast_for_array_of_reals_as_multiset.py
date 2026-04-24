"""
Fast vectorised re-implementation of Simplex-2 embedder.

Mathematically identical pipeline to
C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset.py
up to and including the canonical-form step.  See the Simplex-1 fast
embedder for a description of the shared changes.  The additional S2-
specific change is that per-column sorts are vectorised with np.argsort
on the whole data matrix at once, and MSV indicators for all k components
are built in a single loop over k (unavoidable since each column has its
own sort order).
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
        """
        assert MultisetEmbedder.is_generic_data(data)
        n, k = data.shape
        num_v = k * (n - 1)

        perms = np.argsort(-data, axis=0, kind='stable')
        sorted_data = data[perms, np.arange(k)].astype(np.float64)

        first_diffs_by_col = sorted_data[:-1, :] - sorted_data[1:, :]
        all_deltas = first_diffs_by_col.T.ravel()

        all_msv = np.zeros((num_v, n, k), dtype=np.uint16)
        for c in range(k):
            perm_c = perms[:, c]
            single_c = np.zeros((n, n, k), dtype=np.uint16)
            single_c[np.arange(n), perm_c, c] = 1
            all_msv[c * (n - 1):(c + 1) * (n - 1)] = np.cumsum(single_c, axis=0)[:n - 1]

        perm_sort = np.argsort(-all_deltas, kind='stable')
        sorted_deltas = all_deltas[perm_sort]
        sorted_msv = all_msv[perm_sort]

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
        num_v = k * (n - 1)
        bigN = 2 * (n - 1) * k + 1

        min_elements = data.min(axis=0).astype(np.float64)
        second_diffs, canonical = self._pre_hash_state(data)

        # ------------------------------------------------------------------ #
        # Hash each canonical form to a bigN-dimensional point.               #
        # ------------------------------------------------------------------ #
        points = np.empty((num_v, bigN), dtype=np.float64)
        for p in range(num_v):
            h = hashlib.md5()
            h.update(canonical[p].tobytes())
            h.update(np.array([p + 1], dtype=np.uint16).tobytes())
            seed = int.from_bytes(h.digest()[:4], 'little')
            points[p] = np.random.default_rng(seed).random(bigN)

        core = second_diffs @ points

        length = k + bigN
        assert length == self.size_from_n_k(n, k)
        embedding = np.empty(length, dtype=np.float64)
        embedding[:k] = min_elements
        embedding[k:] = core
        return embedding, None

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        return 2 * n * k + 1 - k
