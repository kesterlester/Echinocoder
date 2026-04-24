#!/usr/bin/env python
# Tests for the fast Simplex-1 embedder.
# Structure mirrors the Simplex-1-specific test calls in test_tools_etc.py
# (the tost_multiset_embedder calls that include embedder_C0HomDeg1_simplex1).
# Random seeds are identical so data is the same as in the original tests.

import numpy as np
import data_sources
from test_tools_etc import tost_multiset_embedder, make_randoms_reproducable
import C0HomDeg1_simplicialComplex_embedder_1_fast_for_array_of_reals_as_multiset as s1f

embedder_s1_fast = s1f.Embedder()


def test_simplex1_fast_various():
    make_randoms_reproducable()

    tost_multiset_embedder(
        data=data_sources.random_real_array_data(mn=(4, 0)),
        embedder=embedder_s1_fast,
    )

    tost_multiset_embedder(
        data=data_sources.random_real_array_data(mn=(0, 4)),
        embedder=embedder_s1_fast,
    )

    tost_multiset_embedder(
        data=data_sources.random_real_array_data(mn=(4, 1)),
        embedder=embedder_s1_fast,
    )

    tost_multiset_embedder(
        data=data_sources.random_real_array_data(mn=(1, 4)),
        embedder=embedder_s1_fast,
    )

    tost_multiset_embedder(
        data=data_sources.random_real_array_data(mn=(2, 3)),
        embedder=embedder_s1_fast,
    )

    tost_multiset_embedder(
        data=data_sources.random_real_array_data(mn=(3, 3)),
        embedder=embedder_s1_fast,
    )

    tost_multiset_embedder(
        data=data_sources.random_real_array_data(mn=(4, 3)),
        embedder=embedder_s1_fast,
    )


if __name__ == "__main__":
    test_simplex1_fast_various()
