import numpy as np
from typing import Any

class MultisetEncoder:
    """
    This is a base class for objects which encode length-n multisets of k-vectors.

    Encoders are not necessarily injective (not necessarily embedders).
    MultisetEmbedder is a subclass of MultisetEncoder: every embedder IS-A encoder,
    but not every encoder is an embedder.

    Strictly speaking we encode "multisets" not "sets" since the containers can hold
    repeated objects and retain knowledge of the number of repeats.  However, we are
    sometimes guilty of abbreviating "multiset" to just "set".

    The multiset to be encoded should be passed in as a 2D numpy array with shape (n,k).
    The order of the vectors within the numpy array can be arbitrary.
    E.g. to encode a multiset containing the 2-vectors (2,2), (4,5) and (1,2) one could call

        encode(np.asarray([[2,2], [4,5], [1,2]]))

    or

        encode(np.asarray([[4,5], [2,2], [1,2]]))

    and both should have the same output -- at least up to numerical precision.
    
    This leeway (permission to have small deviations on account of floating point
    precision, rather than demanding bit-for-bit identical embeddings) is granted 
    to implementations in order to allow them to be faster (sometimes) than would
    be the case if they were all required to canonicalise their input sets.
    Users wanting bit-for-bit identical output under permutations of input
    vectors could easily sort their vectors (in any way) prior to using any 
    encoder.

    All encoders return a tuple comprising:

        (1) a one-dimensional array of real floats,
        (2) the (n,k) size of the data which was encoded,
        (3) None or a metadata packet describing how the data was encoded.

    Derived classes implement the encoding logic via:

        encode_kOne(data, debug)    — for inputs with k==1
        encode_generic(data, debug) — for inputs with n>1 and k>1

    The dispatch (edge cases for n=0, k=0, n=1) is handled by encode() in this base
    class so derived classes do not need to repeat it.

    Derived classes also implement:

        size_from_n_k_generic(n, k) — encoding length for n>1 and k>1

    The full size_from_n_k() (including edge cases) is handled here.
    """

    def encode(self, data: np.ndarray, debug=False) -> (np.ndarray, (int, int), Any):
        n, k = data.shape
        expected_size = self.size_from_n_k(n, k)

        if n < 0 or k < 0:
            raise ValueError("Numpy arrays should not have negative sizes!")
        if n == 0 or k == 0:
            return np.asarray([], dtype=np.float64), (n, k), None
        if n == 1:
            encoding = data.flatten()
            assert len(encoding) == k
            assert len(encoding) == expected_size
            return encoding, (n, k), None
        if k == 1:
            assert self.is_kOne_n_k(n, k)
            encoding, metadata = self.encode_kOne(data, debug)
            assert len(encoding) == n
            assert len(encoding) == expected_size
            return encoding, (n, k), metadata

        assert n > 1 and k > 1
        assert self.is_generic_n_k(n, k)
        encoding, metadata = self.encode_generic(data, debug)
        assert len(encoding) == expected_size
        return encoding, (n, k), metadata

    def encode_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """
        Derived classes should implement this method.
        Encode data for which n>=0 and k==1.  Output length must be exactly n.
        Implementations may assume (without checking) that the input satisfies k==1.

        Common implementations delegate to a static helper, e.g.:

            def encode_kOne(self, data, debug=False):
                from MultisetEmbedder import MultisetEmbedder
                return MultisetEmbedder.embed_kOne_polynomial(data), None  # Cinf

        or

            def encode_kOne(self, data, debug=False):
                from MultisetEmbedder import MultisetEmbedder
                return MultisetEmbedder.embed_kOne_sorting(data), None     # C0
        """
        raise NotImplementedError()

    def encode_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """
        Derived classes should implement this method.
        Encode data for which n>1 and k>1.
        Implementations may assume (without checking) that the input satisfies n>1 and k>1.
        """
        raise NotImplementedError()

    def size_from_n_k(self, n: int, k: int) -> int:
        """
        Returns the number of reals in the encoding of a multiset of n k-vectors.
        Returns -1 if this encoder cannot encode inputs of that shape.
        """
        if n < 0 or k < 0:
            return -1
        if n == 0 or k == 0:
            return 0
        if n == 1:
            return k
        if k == 1:
            return n
        return self.size_from_n_k_generic(n, k)

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        """
        Derived classes should implement this method.
        Returns the encoding size for inputs with n>1 and k>1.
        Implementations may assume (without checking) that n>1 and k>1.
        Return -1 if inputs of that shape are not encodable by this encoder.
        """
        raise NotImplementedError()

    def size_from_array(self, data: np.ndarray) -> int:
        """
        Returns the number of reals that would result from encoding the multiset
        represented by data.  Returns -1 if data of that shape is not encodable.
        """
        n, k = data.shape
        return self.size_from_n_k(n, k)

    @staticmethod
    def is_kOne_data(data: np.ndarray) -> bool:
        n, k = data.shape
        return MultisetEncoder.is_kOne_n_k(n, k)

    @staticmethod
    def is_kOne_n_k(n: int, k: int) -> bool:
        return n >= 0 and k == 1

    @staticmethod
    def is_generic_data(data: np.ndarray) -> bool:
        n, k = data.shape
        return MultisetEncoder.is_generic_n_k(n, k)

    @staticmethod
    def is_generic_n_k(n: int, k: int) -> bool:
        return n > 1 and k > 1
