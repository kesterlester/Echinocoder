import numpy as np
from typing import Any
import Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset as poly_list
from MultisetEncoder import MultisetEncoder

class MultisetEmbedder(MultisetEncoder):
    """
    This is a base class for objects which embed (injectivity-guaranteed encode)
    length-n multisets of k-vectors.

    MultisetEmbedder IS-A MultisetEncoder.  All embedders are encoders, but not all
    encoders are embedders.  Use MultisetEncoder as the base class for algorithms
    that do not guarantee injectivity.

    Strictly speaking these are "multisets" not "sets" since the sets can hold repeated
    objects and retain knowledge of the number of repeats.  However, we are sometimes
    guilty of abbreviating "multiset" to just "set".

    The set to be embedded should be passed in as a 2D numpy array with shape (n,k).
    The order of the vectors within the numpy array can be arbitrary.
    E.g. to embed a multiset containing the 2-vectors (2,2), (4,5) and (1,2) one could call

        embed(np.asarray([[2,2], [4,5], [1,2]]))

    or

        embed(np.asarray([[4,5], [2,2], [1,2]]))

    and both should have the same output -- at least up to numerical precision. This leeway
    (permission to have small deviations on account of floating point precision, rather than
    demanding bit-for-bit identical embeddings) is granted to implementations in order to
    allow them to be faster (sometimes) than would be the case if they were all required to
    canonicalise their input sets.  Someone wanting bit-for-bit identical output under
    permutations of input vectors could easily sort their vectors (in any way) prior to
    using any embedder.

    All embedders return a tuple, comprising:

        (1) a one-dimensional array of real floats,
        (2) the (n,k) size of the data which was encoded
        (3) None or a metadata packet describing how the data was encoded.

    embed() is an alias for encode(); the dispatch logic lives in MultisetEncoder.encode().
    Derived classes implement embedding via:

        embed_kOne(data, debug)    — for inputs with k==1
        embed_generic(data, debug) — for inputs with n>1 and k>1

    These are bridged automatically to the encode_kOne / encode_generic interface that
    MultisetEncoder.encode() expects.

    Derived classes implement:

        size_from_n_k_generic(n, k) — encoding/embedding length for n>1 and k>1

    In principle, a given embedder can embed sets of different sizes n and/or k.
    However, some embedders might wish to restrict themselves to certain fixed n or k at
    initialisation.  Derived classes report what inputs they can handle via size_from_n_k().
    A return value of -1 means that shape is not supported.
    """

    def embed(self, data: np.ndarray, debug=False) -> (np.ndarray, (int, int), Any):
        """embed() is an alias for encode().  The dispatch lives in MultisetEncoder."""
        return self.encode(data, debug)

    # ------------------------------------------------------------------
    # Bridges from MultisetEncoder's abstract interface to MultisetEmbedder's
    # embed_kOne / embed_generic interface.
    # ------------------------------------------------------------------

    def encode_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """Satisfies MultisetEncoder contract by delegating to embed_kOne."""
        return self.embed_kOne(data, debug)

    def encode_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """Satisfies MultisetEncoder contract by delegating to embed_generic."""
        return self.embed_generic(data, debug)

    # ------------------------------------------------------------------
    # Abstract methods for embedder subclasses to implement.
    # ------------------------------------------------------------------

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """
        Derived classes should implement this method.
        This method should OPTIMALLY embed data for which n>=0 and k==1.
        OPTIMALLY means that the embedding size must therefore be n.
        Implementations may assume (without checking) that data fed to it has k==1.
        It is likely that most implementations will implement this method either as:

            def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
                return MultisetEmbedder.embed_kOne_sorting(data), None   # Want C0 piecewise linear

        or as:

            def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
                return MultisetEmbedder.embed_kOne_polynomial(data), None  # Want Cinf

        depending on whether they want piecewise linearity or differentiability.
        """
        raise NotImplementedError()

    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """
        Derived classes should implement this method.
        This method should embed data for which n>1 and k>1.
        Implementations may assume (without checking) that the input satisfies n>1 and k>1.
        """
        raise NotImplementedError()

    # ------------------------------------------------------------------
    # Sanity-check helper (for use in unit tests of derived classes).
    # ------------------------------------------------------------------

    def test_me(self):
        _ = self.size_from_n_k_generic(2, 2)
        _ = self.embed_generic(np.array([[1, 2], [3, 4]], dtype=np.float64))
        _ = self.embed_kOne(np.array([[1,], [3,],], dtype=np.float64))

    # ------------------------------------------------------------------
    # Static helpers for use by embed_kOne implementations.
    # These live here rather than on MultisetEncoder because they
    # produce embeddings (injective encodings), not mere encodings.
    # ------------------------------------------------------------------

    @staticmethod
    def embed_kOne_sorting(data: np.ndarray) -> np.ndarray:
        assert MultisetEmbedder.is_kOne_data(data)
        return np.sort(data.flatten())

    @staticmethod
    def embed_kOne_polynomial(data: np.ndarray) -> np.ndarray:
        assert MultisetEmbedder.is_kOne_data(data)
        embedding, shape_, metadata_ = poly_list.embed(data.flatten())
        return embedding
