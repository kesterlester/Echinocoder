"""Algorithm 2 implemented using EncDec.py as a service library.

This is the "2b" variant.  Steps A–E are delegated to
EncDec.simplex_2_preprocess_steps; Step F (hash each canonical basis matrix to a
point in the unit hypercube, then form the weighted sum) is implemented here.

The default configuration (preserve_scale_in_step_1=False,
preserve_scale_in_step_2=True) is designed to be numerically equivalent to the
"2a" direct implementation: both yield the same embedding vector (to floating-point
precision) for the same input.  See test_simplex2_ab_parity.py.

Other configurations are valid and produce different (but still correct) embeddings:
  preserve_scale_in_step_1=True  — per-row gaps are normalised before merging
  preserve_scale_in_step_2=False — second Abel step stores Δ/L instead of α/L̂

See also:
  _2a_ — original direct implementation (the reference for parity tests)
  _2_  — wrapper that selects between _2a_ and _2b_
  EncDec.py — service library that covers Steps A–E
"""

import hashlib
from fractions import Fraction

import numpy as np
from typing import Any

import EncDec
from MultisetEmbedder import MultisetEmbedder


class Embedder(MultisetEmbedder):

    #: Default config matches _2a_: raw first-Abel gaps, scale-preserved second Abel.
    DEFAULT_PRESERVE_SCALE_STEP_1 = False
    DEFAULT_PRESERVE_SCALE_STEP_2 = True

    def __init__(self,
                 preserve_scale_in_step_1: bool = DEFAULT_PRESERVE_SCALE_STEP_1,
                 preserve_scale_in_step_2: bool = DEFAULT_PRESERVE_SCALE_STEP_2):
        """
        Args:
            preserve_scale_in_step_1: If True, the first (per-row) Abel summation
                normalises prefix sums by their count, producing L̂^(i)_r = L^(i)_r/(r+1)
                as Fraction entries.  If False (default, matching _2a_), raw cumulative
                integer sums are used.
            preserve_scale_in_step_2: If True (default, matching _2a_), the second
                (global) Abel summation normalises, yielding α_s coefficients and
                L̂_{(s)} basis vectors.  If False, Δ_{(s)} coefficients and raw L_{(s)}
                matrices are stored instead.
        """
        self.preserve_scale_in_step_1 = preserve_scale_in_step_1
        self.preserve_scale_in_step_2 = preserve_scale_in_step_2

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        return MultisetEmbedder.embed_kOne_sorting(data), metadata

    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        assert MultisetEmbedder.is_generic_data(data)
        n, k = data.shape

        # ── Steps A–E via EncDec ──────────────────────────────────────────────
        lin_comb, offsets = EncDec.simplex_2_preprocess_steps(
            data,
            preserve_scale_in_step_1=self.preserve_scale_in_step_1,
            preserve_scale_in_step_2=self.preserve_scale_in_step_2,
            canonicalise=True,
            use_assertions=True,
            debug=debug,
        )

        # lin_comb holds k(n-1) = nk-k canonical (coeff, basis_vec) pairs:
        #   preserve_scale_step2=True  → coeff = α_s, basis_vec = L̂_{(s)} (Fraction array)
        #   preserve_scale_step2=False → coeff = Δ_{(s)}, basis_vec = L_{(s)} (int array)
        expected_n_vertices = n * k - k
        assert len(lin_comb) == expected_n_vertices

        # ── Step F: hash + weighted sum ───────────────────────────────────────
        bigN = 2 * (n - 1) * k + 1  # hash embedding dimension

        second_part_of_embedding = sum(
            [float(coeff) * _hash_basis_vec_to_point(
                basis_vec, s + 1, bigN, self.preserve_scale_in_step_2)
             for s, (coeff, basis_vec) in enumerate(zip(lin_comb.coeffs, lin_comb.basis_vecs))]
        ) + np.zeros(bigN)

        # ── Assemble output ───────────────────────────────────────────────────
        # Extract per-row minima from offsets.
        # offsets is a LinComb with k terms; offsets.coeffs[i] = m_i (with
        # preserve_scale_step1=False) or n·m_i (with True).
        if self.preserve_scale_in_step_1:
            min_elements = np.array([float(c) / n for c in offsets.coeffs])
        else:
            min_elements = np.array([float(c) for c in offsets.coeffs])

        length_of_embedding = self.size_from_n_k(n, k)
        assert length_of_embedding == bigN + k

        embedding = np.zeros(length_of_embedding, dtype=np.float64)
        embedding[:k] = min_elements
        embedding[k:k + bigN] = second_part_of_embedding

        if debug:
            print(f"min_elements = {min_elements}")
            print(f"second_part_of_embedding = {second_part_of_embedding}")
            print(f"embedding = {embedding}")

        metadata = None
        return embedding, metadata

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        return 2 * n * k + 1 - k


# ── Step-F helper ─────────────────────────────────────────────────────────────

def _hash_basis_vec_to_point(
        canonical_basis_vec, index: int, dimension: int, preserve_scale_step2: bool
) -> np.ndarray:
    """Hash a canonical k×n basis matrix to a point in the unit hypercube.

    Produces the same result as Eji_LinComb.hash_to_point_in_unit_hypercube
    would produce for the Eji_LinComb object with _index=index and
    _eji_counts = L_{(s)} (the raw integer cumulative-count matrix).

    Args:
        canonical_basis_vec: the basis vector from EncDec's canonicalised lin_comb.
            If preserve_scale_step2=True: entries are Fraction objects equal to
            L_{(s)}[j,i] / index.  Integer counts are recovered by multiplying by index.
            If preserve_scale_step2=False: entries are integers equal to L_{(s)}[j,i].
        index: s+1, i.e. the Eji_LinComb _index value.
        dimension: output vector length.
        preserve_scale_step2: controls how to recover the integer counts (see above).
    """
    bv = np.asarray(canonical_basis_vec)

    if preserve_scale_step2:
        # Entries are Fraction(count / index); multiply by index to get exact integers.
        vf = np.vectorize(lambda x: int(Fraction(x) * index))
        counts = vf(bv).astype(np.uint16)
    else:
        # Entries are already integer cumulative counts.
        counts = bv.astype(np.uint16)

    counts = np.ascontiguousarray(counts, dtype=np.uint16)

    m = hashlib.md5()
    m.update(counts.tobytes())
    m.update(np.uint16(index).tobytes())
    seed = int.from_bytes(m.digest(), 'big')
    return np.random.default_rng(seed).random(dimension)


# ── Self-test ─────────────────────────────────────────────────────────────────

def tost():
    """Smoke-test: embed a small matrix and check the output has the right shape."""
    data = np.asarray([[4, 2, 3], [-3, 5, 1], [8, 9, 2], [2, 7, 2]])
    embedder = Embedder()
    result, _ = embedder.embed(data, debug=False)
    n, k = data.shape
    assert result.shape == (embedder.size_from_n_k(n, k),), f"Unexpected shape: {result.shape}"


def run_unit_tests():
    tost()


if __name__ == "__main__":
    run_unit_tests()

    default_test_input = np.asarray([[4, 2, 3],
                                     [-3, 5, 1],
                                     [8, 9, 2],
                                     [2, 7, 2]])
    embedder = Embedder()
    output = embedder.embed(default_test_input, debug=True)
    print("Embedding:")
    print(default_test_input)
    print("leads to:")
    print(output)
