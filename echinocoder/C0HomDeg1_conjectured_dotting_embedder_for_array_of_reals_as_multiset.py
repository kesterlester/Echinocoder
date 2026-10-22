import numpy as np
from tools import sort_each_np_array_row
#import hashlib
from MultisetEmbedder import MultisetEmbedder
import C0HomDeg1_dotting_encoder_for_array_of_reals_as_multiset as dotting_encoder
import math
from typing import Any

class Embedder(MultisetEmbedder):
    """
    This encoder encodes via sorted length-n lists of dot products.
    The first k-such lists are with standard coordinate axes.
    After that come "extra_dots" extra dot products. extra_dots is set in the constructor.
    """

    def __init__(self, n: int, k: int):
        super().__init__()
        if k<0 or n<0:
            raise ValueError() # It makes no sense to claim we will code vectors or sets for k<0.

        if k==0 or n==0:
            raise ValueError() # We choose not to undertake to embed non-data! We could, but we would need to manage many special cases internally which really ought to the the responsibility for external users. They themselves have no need to embed non-data!

        M_min = (k-1)*(math.floor(math.log2(n) ) + 1) + 1 # Conjuectured in OneNote.Symmetries.DOT_CONFUSABLE to be sufficient
        extra_dots = M_min - k 
        assert extra_dots >= 0
        assert k==1 or extra_dots>0

        self.n = n
        self.k = k
        self._dotting_encoder  = dotting_encoder.Encoder(k=k, extra_dots=extra_dots)

    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any): # was embed(...).  fixed 21st Feb 2025.
        expected_size = self.size_from_array(data)
        if expected_size == -1:
            raise ValueError(f"We do not undertake to embed data of shape {data.shape}")
        if debug:
            print(f"data is {data}")
        embedding, size, metadata = self._dotting_encoder.encode(data)
        if debug:
            print(f"Embedding is {embedding} with length {len(embedding)}")
        assert len(embedding) == expected_size
        metadata = None
        return embedding, metadata
    
    def size_from_n_k_generic(self, n: int, k: int) -> int:
        if n!=self.n or k!=self.k:
            return -1 # We are optimised to work with a certain n and a certain k
        return self._dotting_encoder.size_from_n_k(n,k)

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        return MultisetEmbedder.embed_kOne_sorting(data), metadata


def tost(): # Renamed from test -> tost to avoid pycharm mis-detecting / mis-running unit tests!

    embedder = Embedder(n=5, k=2)
    print("Embedder matrix is\n",embedder._dotting_encoder.matrix)
    assert embedder.size_from_n_k(5,2) == 5*(2+2) # For k=2 and n=4 a single dot can be ambiguous so need at least two dots.  Conjecture is that two dots suffices for k=2 when n<8

    try:
        bad_input = np.asarray([[]]) # non data
        _ = embedder.embed(bad_input)
        assert True
    except ValueError:
        assert False

    try:
        bad_input = np.asarray([[4,2,2],[-3,5,1],[8,9,0],[2,7,4],[3,2,1]]) # k=3 data but k=2 embedder
        _ = embedder.embed(bad_input)
        assert False # We should not get here! Last line should throw Value Error
    except ValueError:
        assert True

    try:
        bad_input = np.asarray([[4,2],[-3,5],[8,9],[2,7],[3,2],[7,7]]) # n=6 data but n=5 embedder
        _ = embedder.embed(bad_input)
        assert False # We should not get here! Last line should throw Value Error
    except ValueError:
        assert True

    try:
        good_input = np.asarray([[4,2],[-3,5],[8,9],[2,7],[3,2]])
        _ = embedder.embed(good_input)
        assert True # We should not get here! Last line should throw Value Error
    except ValueError:
        assert False

def run_unit_tests():
    tost() # Renamed from test -> tost to avoid pycharm mis-detecting / mis-running unit tests!

def main():
    run_unit_tests()


    embedder = Embedder(n=4, k=2)
    good_input = np.asarray([[4,2],[-3,5],[8,9],[2,7]])
    output = embedder.embed(good_input, debug=True)

    print("Embedding:")
    print(f"{good_input}")
    print("leads to:")
    print(f"{output}")

if __name__ == "__main__":
    main()
