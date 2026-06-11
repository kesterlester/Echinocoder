# Demonstrates that Cinf_numpy_polynomial_embedder_for_array_of_reals_as_multiset
# is NOT injective, and therefore over-promises by inheriting from MultisetEmbedder.
#
# We exhibit two multisets M1 and M2 that are distinct (verified below as multisets,
# not merely as ordered arrays) yet map to the same encoding.
#
# The algebraic reason: the algorithm encodes by projecting onto every pair of columns
# and polynomial-embedding the resulting complex multiset.  M1 and M2 are constructed
# so that every such column-pair projection yields the same multiset of complex numbers,
# even though the full 3-vectors differ.

import numpy as np
import pytest
import Cinf_numpy_polynomial_embedder_for_array_of_reals_as_multiset as poly_array

# Six distinct integer values stand in for A, B, X, Y, u, v.
A, B, X, Y, u, v = 1, 2, 3, 4, 5, 6

# M1 = { (A,X,u), (A,Y,v), (B,X,v), (B,Y,u) }
M1 = np.array([
    [A, X, u],
    [A, Y, v],
    [B, X, v],
    [B, Y, u],
], dtype=float)

# M2 = { (A,X,v), (A,Y,u), (B,X,u), (B,Y,v) }
M2 = np.array([
    [A, X, v],
    [A, Y, u],
    [B, X, u],
    [B, Y, v],
], dtype=float)


def multiset_of_rows(arr):
    """Return the rows of a 2D array as a sorted tuple-of-tuples: a canonical multiset representation."""
    return tuple(sorted(tuple(row) for row in arr))


def test_M1_and_M2_are_distinct_multisets():
    """M1 and M2 are genuinely different multisets — not just different row-orderings of the same elements."""
    assert multiset_of_rows(M1) != multiset_of_rows(M2)


def test_polynomial_array_embedder_maps_M1_and_M2_to_same_encoding():
    """
    Despite M1 != M2 as multisets, the polynomial array embedder produces the same output
    for both, proving it is not injective and should not be called an embedder.
    """
    embedder = poly_array.Embedder()
    encoding_M1, _shape_M1, _meta_M1 = embedder.embed(M1)
    encoding_M2, _shape_M2, _meta_M2 = embedder.embed(M2)
    np.testing.assert_allclose(encoding_M1, encoding_M2, rtol=1e-12, atol=1e-12)
