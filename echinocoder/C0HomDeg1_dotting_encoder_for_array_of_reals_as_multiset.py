import numpy as np
from tools import sort_each_np_array_row
#import hashlib
from MultisetEncoder import MultisetEncoder
from typing import Any


class Encoder(MultisetEncoder):
    """
    This encoder encodes via sorted length-n lists of dot products.
    The first k-such lists are with standard coordinate axes.
    After that come "extra_dots" extra dot products. extra_dots is set in the constructor.
    """

    def __init__(self, k: int, extra_dots: int):
        super().__init__()
        if k<0 or extra_dots<0:
            raise ValueError() # It makes no sense to claim we will code vectors or sets for k<0 or extra_dots<0

        self.k = k
        self.extra_dots = extra_dots
        self.total_dots = k + extra_dots

        rng = np.random.Generator(np.random.PCG64(0))

        self.matrix = np.concatenate((np.identity(k), rng.standard_normal((self.extra_dots, self.k))), axis=0) # Technically this matrix should then be checked -- to see that none of its kxk minors have zero determinant. But for now we "assume" that that is true.

    def encode(self, data: np.ndarray, debug=False) -> (np.ndarray, (int, int), Any):
        if debug:
            print(f"data is {data}")
    
        n,k = data.shape

        # Catch a few special cases:
        if n==0 or k==0:
            encoding = np.array([], dtype=np.float64)
            assert len(encoding) == self.size_from_n_k(n, k)
            return encoding, (n,k), None

        if k != self.k:
            assert self.size_from_n_k(n,k) == -1
            raise ValueError(f"This encoder is setup for k={self.k} so does not like data having k={k}")

        if debug:
            print(f"About to multiply {self.matrix} by {data.T}")

        ##### HERE IS THE ACTUAL ENCODING:     ###########
        encoding = sort_each_np_array_row(self.matrix @ data.T).flatten()
        ##### THE ACTUAL ENCODING IS COMPLETE! ###########

        if debug:
            print(f"Encoding is {encoding} with length {len(encoding)} which should be {self.size_from_n_k(n,k)}")
        assert len(encoding) == self.size_from_n_k(n, k)
        return encoding, (n,k), None

    def decode(self, encoding, encoding_n_and_k):
        return None

    def size_from_n_k(self, n: int, k: int) -> int:
        if k==0 or n==0:
            return 0
        if k != self.k:
            return -1
        return n*(k+self.extra_dots)



## def hash_to_128_bit_md5_int(md5):
##     return int.from_bytes(md5.digest(), 'big') # 128 bits worth.
## 
## def hash_to_64_bit_reals_in_unit_interval(md5):
##     """
##     An md5 sum is 64 bits long so we get two such reals.
##     N.B. This hash is of self._eji_counts only -- i.e. it ignores self._index.
##     For the purposes to which this hash will be used, that is believed to be apporopriate.
##     """
## 
##     x = hash_to_128_bit_md5_int(md5)
##     bot_64_bits = x & 0xffFFffFFffFFffFF
##     top_64_bits = x >> 64
##     return np.float64(top_64_bits)/(1 << 64), np.float64(bot_64_bits)/(1 << 64)



def tost(): # Renamed from test -> tost to avoid pycharm mis-detecting / mis-running unit tests!
    encoder = Encoder(k=2,extra_dots=6)
    print("Encoder matrix is\n",encoder.matrix)
    assert encoder.size_from_n_k(5,2) == 5*(2+6)
    #calculated = np.array([2, 3, 4, 1, 0])
    #expected = np.array([2, 3, 4, 1, 0])
    #np.testing.assert_array_equal(calculated, expected)


def run_unit_tests():
    tost() # Renamed from test -> tost to avoid pycharm mis-detecting / mis-running unit tests!

def main():
    run_unit_tests()

    encoder = Encoder(k=2, extra_dots=10)
    good_input = np.asarray([[4,2],[-3,5],[8,9],[2,7]])
    output = encoder.encode(good_input, debug=True)

    print("Encoding:")
    print(f"{good_input}")
    print("leads to:")
    print(f"{output}")

if __name__ == "__main__":
    main()
