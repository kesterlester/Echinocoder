import numpy as np
import Cinf_numpy_polynomial_encoder_for_array_of_reals_as_multiset as poly_encoder
import Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset     as burs_encoder
import Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset  as poly_list
from MultisetEncoder import MultisetEncoder
from typing import Any

class Encoder(MultisetEncoder):
    """
    This encoder uses whichever of the bursar or polynomial encoders would be most
    efficient for the multiset to be encoded.

    NOTE: This encoder is NOT injective and therefore NOT an embedder.
    The polynomial array encoder it delegates to is not injective (see
    Cinf_numpy_polynomial_encoder_for_array_of_reals_as_multiset.py), so neither is this.
    """

    def __init__(self):
        super().__init__()
        self._poly_encoder = poly_encoder.Encoder()
        self._burs_encoder = burs_encoder.Embedder()

    def encode_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        poly_size = self._poly_encoder.size_from_array(data)
        burs_size = self._burs_encoder.size_from_array(data)

        if burs_size == -1 and poly_size == -1:
            raise ValueError()
        elif poly_size <= burs_size:
            if debug: print(f"Hybrid uses poly encoder for data of shape {data.shape}.")
            encoding, size_, metadata = self._poly_encoder.encode(data)
        else:
            if debug: print(f"Hybrid uses burs encoder for data of shape {data.shape}.")
            encoding, size_, metadata = self._burs_encoder.encode(data)

        assert len(encoding) == self.size_from_array(data)
        return encoding, metadata

    def encode_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        embedding, _shape, _metadata = poly_list.embed(data.flatten())
        return embedding, metadata

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        poly_size = self._poly_encoder.size_from_n_k(n, k)
        burs_size = self._burs_encoder.size_from_n_k(n, k)

        if burs_size == -1 and poly_size == -1:
            return -1
        elif poly_size <= burs_size:
            return poly_size
        else:
            return burs_size


def tost():  # Renamed from test -> tost to avoid pytest mis-detection.
    encoder = Encoder()

    poly_input = np.asarray([[4,2],[-3,5],[8,9],[2,7],[3,2]])
    _ = encoder.encode(poly_input, debug=True)

    burs_input = np.asarray([[4,2,-3,5,8],[9,2,7,3,2]])
    _ = encoder.encode(burs_input, debug=True)


def run_unit_tests():
    tost()

def main():
    run_unit_tests()

    encoder = Encoder()
    good_input = np.asarray([[4,2],[-3,5],[8,9],[2,7]])
    output = encoder.encode(good_input, debug=False)

    print("Encoding:")
    print(f"{good_input}")
    print("leads to:")
    print(f"{output}")

if __name__ == "__main__":
    main()
