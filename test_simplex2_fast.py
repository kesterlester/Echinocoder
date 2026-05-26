#!/usr/bin/env python
# Tests for the fast Simplex-2 embedder.
# Structure mirrors test_simplex2.py; input data is identical.
# Expected embedding values differ from the original because the fast
# implementation uses a seeded numpy RNG rather than sequential MD5 calls
# to map canonical forms to points in the unit hypercube.
# Metadata is None (the fast embedder does not return ascending_data).

import numpy as np
from C0HomDeg1_simplicialComplex_embedder_2_fast_for_array_of_reals_as_multiset import Embedder


def test_simplex2_fast():
    data_1 = np.asarray([[4, 2, 3],
                         [-3, 5, 1],
                         [8, 9, 2],
                         [2, 7, 2]])

    data_2 = np.asarray([[4, 3, 3],
                         [-3, 5, 1],
                         [8, 9, 2],
                         [2, 6, 2]])

    embedder = Embedder()
    output_1 = embedder.embed(data_1, debug=False)
    output_2 = embedder.embed(data_2, debug=False)

    expected_embedding_1 = np.array([
        -3.               ,  2.               ,  1.               ,
         6.516753543250672,  9.617185361165806,  8.203659342823979,
        11.294730972084668, 11.166881104593841, 15.571013536870456,
        10.695647291780032, 11.580449444011997,  7.400121089749579,
         5.591885621586032,  4.927988326265223, 10.54571021127005 ,
         2.769626951234436,  6.563177374293916,  9.126501445354105,
         9.983954338957364, 14.248726423417128, 11.447734758481666,
         7.022294520873803,
    ])

    expected_embedding_2 = np.array([
        -3.                ,  3.                ,  1.                ,
        10.707898770016362 ,  7.992476998084883 ,  6.219729235616787 ,
         6.276224107244218 ,  9.808735702706546 , 14.773359568663793 ,
        14.458370950062275 ,  9.95942661369931  ,  7.60838325963325  ,
         7.4368498191447125,  6.758314126254907 ,  9.523269204656478 ,
         6.057683816131951 ,  5.277094720587709 ,  5.165988996523576 ,
         9.115501796889328 , 14.066831759088668 ,  9.848538014149984 ,
         9.56249885662044  ,
    ])

    for output, expected_embedding in (
        (output_1, expected_embedding_1),
        (output_2, expected_embedding_2),
    ):
        embedding, shape, metadata = output

        print(f"output    embedding={embedding}")
        print(f"expected  embedding={expected_embedding}")

        np.testing.assert_array_equal(embedding, expected_embedding)
        assert shape == (4, 3)
        assert metadata is None

    # Different inputs must produce different embeddings.
    assert not np.array_equal(output_1[0], output_2[0])


if __name__ == "__main__":
    test_simplex2_fast()
