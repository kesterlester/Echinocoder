# Encode array of real numbers treated as multiset.
# E.g. this method can encode things like:
#
#        [[3,4],[4,2],[-2,1]]
#
# when used to represent this multiset:
#
#        {{ [3,4],[4,2],[-2,1] }}
#
# total_encoding_size = m*(m-1)*n
#
# where n is the number of vectors in the multiset and m is their dimension.
# E.g. in the example above m=2 and n=3.
#
# NOTE: This encoder is NOT injective and therefore NOT an embedder.
# See test_Cinf_numpy_polynomial_encoder_for_array_of_reals_as_multiset_not_injective.py
# for an explicit counterexample demonstrating the failure of injectivity.

from MultisetEncoder import MultisetEncoder
import numpy as np
import Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset
import itertools
import tools
from typing import Any

class Encoder(MultisetEncoder):

    def encode_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        embedding, _shape, _metadata = Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.embed(data.flatten())
        return embedding, metadata

    def encode_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata_out = None

        assert MultisetEncoder.is_generic_data(data)
        n, m = data.shape
        assert n > 1 and m > 1

        assert m != 1

        if m == 2:
            # The "vectors" are 2-long, so turn them into complex numbers, encode, and expand.
            complex_list = data[:,0] + complex(0,1)*data[:,1]
            complex_encoding, shape_, metadata_ = Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.embed(complex_list)
            real_encoding = tools.expand_complex_to_real_pairs(complex_encoding)
            assert len(real_encoding) == self.size_from_n_k(n, m)
            return real_encoding, metadata_out

        # The "vectors" are 3-or-more-long, so take all m-choose-2 pairings of elements
        # and encode each pair using the m==2 case above.
        number_of_pairs = m*(m-1)//2
        current_encoding_index = 0
        for row1_index, row2_index in itertools.combinations(range(m), 2):
            pair_of_rows = data[:, [row1_index, row2_index]]
            encoding_for_pair, metadata_3 = self.encode_generic(pair_of_rows, debug)
            if current_encoding_index == 0:
                pair_encoding_size = len(encoding_for_pair)
                total_encoding_size = pair_encoding_size * number_of_pairs
                real_encoding = np.empty(total_encoding_size, dtype=encoding_for_pair.dtype)
            real_encoding[current_encoding_index:current_encoding_index+pair_encoding_size] = encoding_for_pair
            current_encoding_index += pair_encoding_size

        assert len(real_encoding) == n*m*(m-1)
        assert len(real_encoding) == self.size_from_n_k(n, m)
        return real_encoding, metadata_out

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        return n*k*(k-1)
