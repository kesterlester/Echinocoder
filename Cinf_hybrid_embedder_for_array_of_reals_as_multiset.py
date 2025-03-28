import numpy as np
from MultisetEmbedder import MultisetEmbedder
import Cinf_numpy_polynomial_embedder_for_array_of_reals_as_multiset as poly_encoder
import Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset     as burs_encoder
from typing import Any

class Embedder(MultisetEmbedder):
    """
    This encoder uses whichever of the bursarial or polynomial encoders would be modt efficient
    for the set to be encoded.
    """

    def __init__(self):
        super().__init__()
        self._poly_encoder = poly_encoder.Embedder()
        self._burs_encoder = burs_encoder.Embedder()

    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        poly_size = self._poly_encoder.size_from_array(data)
        burs_size = self._burs_encoder.size_from_array(data)
        
        if burs_size == -1 and poly_size == -1:
            raise ValueError()
        elif poly_size <= burs_size:
            if debug: print(f"Hybrid uses poly embedder for data of shape {data.shape}.")
            embedding, size_, metadata = self._poly_encoder.embed(data)
        else:
            if debug: print(f"Hybrid uses burs embedder for data of shape {data.shape}.")
            embedding, size_, metadata = self._burs_encoder.embed(data)

        assert len(embedding) == self.size_from_array(data)
        return embedding, metadata
   
    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        return MultisetEmbedder.embed_kOne_polynomial(data), metadata

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        poly_size = self._poly_encoder.size_from_n_k(n,k)
        burs_size = self._burs_encoder.size_from_n_k(n,k)
        
        if burs_size == -1 and poly_size == -1:
            return -1
        elif poly_size <= burs_size:
            return poly_size
        else:
            return burs_size


def tost(): # Renamed from test -> tost to avoid pycharm mis-detecting / mis-running unit tests!
    embedder = Embedder()

    poly_input = np.asarray([[4,2],[-3,5],[8,9],[2,7],[3,2]])
    _ = embedder.embed(poly_input, debug=True)

    burs_input = np.asarray([[4,2,-3,5,8],[9,2,7,3,2]])
    _ = embedder.embed(burs_input, debug=True)


def run_unit_tests():
    tost() # Renamed from test -> tost to avoid pycharm mis-detecting / mis-running unit tests!

def main():
    run_unit_tests()


    embedder = Embedder()
    good_input = np.asarray([[4,2],[-3,5],[8,9],[2,7]])
    output = embedder.embed(good_input, debug=False)

    print("Embedding:")
    print(f"{good_input}")
    print("leads to:")
    print(f"{output}")

if __name__ == "__main__":
    main()
