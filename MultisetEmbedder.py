import numpy as np
from typing import Any
import Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset as poly_list

class MultisetEmbedder:
    """
    This is a base class for objects which embed length-n sets of k-vectors.

    Encoders are not necessarily embedders since encoders do not need to be injective. 
    All embedders are encoders, however.

    Strictly speaking these are "multisets" not "sets" since the sets can hold repeated 
    objects and retain knowledge of the number of repeats.  However, we are sometimes guity 
    of abbreviating "multiset" to just "set".

    The set to be embedded should be passed in as a 2D numpy array with shape (n,k).
    The order of the vectors within the numpy array can be arbitrary
    E.g. to embed a multiset containing the 2-vectors (2,2), (4,5) and (1,2) one could call
    
        embed(np.asarray([[2,2], [4,5], [1,2]]))

    or

        embed(np.asarray([[4,5], [2,2], [1,2]]))

    and both should have the same output -- at least up to numerical precision. This leeway 
    (permission to have small deviations on account of floating point precision, rather than
    demanding bit-for-bit identical embeddings) is granted to implementations in order to 
    allow them to be faster (sometimes) than would be the case if they were all required to
    canonicalise their input sets.  Someone wanting bit-for-bit identical output under permu-
    tations of input vectors could easily sort their vectors (in any way) prior to using any 
    embedder.

    All embedders return a tuple, comprising:

        (1) a one-dimensional arrays of real floats,
        (2) the (n,k) size of the data which was encoded
        (3) None or a meta-data packet describing how the data was encoded. 

    In principle, a given embedder can embed sets of different sizes n and or k.
    However, some embedders might wish to restrict themselves to certain fixed n or k at
    initialisation (e.g. if an embedder were to need a significant amount of n and k dependent set-up
    cost that it wished to do once only).  Thus all embedders are expected to be able to tell 
    callers what sizes of input they can and cannot embed, and how long the possible embeddings 
    will be.  Derived classes do this by implementing the method

        size_from_n_k(n: int, k, int) -> int:

    for which:
    
        * a return value of >=0 is the number of reals in an embedding if sets of size (n,k) are encodable, and
        * a return value of -1 indicates that embedding for that n and k is impossible.

    No embedder should have to deal with $k<0$ or $n<0$ as these are crazy nonsense.

    """

    def embed(self, data: np.ndarray, debug=False) -> (np.ndarray, (int, int), Any):

        n,k = data.shape
        expected_order = self.size_from_n_k(n,k)

        if n<0 or k<0:
            raise ValueError("Numpy array's should not have negative sizes!!!")
        if n==0 or k==0:
            return np.asarray([], dtype=np.float64), (n,k), None
        if n==1:
            embedding = data.flatten() # This implmentation is a coverall.
            assert len(embedding) == k
            assert len(embedding) == expected_order
            return embedding, (n,k), None
        if k==1:
            assert k==1 and n>=0 # Preconditions for calling self.embed_kOne !
            assert self.is_kOne_n_k(n, k) # Precondition for calling self.embed_kOne !
            embedding, metadata = self.embed_kOne(data, debug) # Derived classes should implement this method!
            assert len(embedding) == n # Derived classes are required to meet this condition in their output!
            assert len(embedding) == expected_order
            return embedding, (n,k), metadata

        assert n>1 and k>1 # Preconditions for calling self.embed_generic !
        assert self.is_generic_n_k(n,k) # Precondition for calling self.embed_generic !
        embedding, metadata = self.embed_generic(data, debug) # Derived classes should implement this method!
        assert len(embedding) == expected_order
        return embedding, (n,k), metadata

    def size_from_n_k(self, n: int, k: int) -> int:
        """
        This function returns the number of reals that the embedding would contain if a set containing 
        n "k-vectors" were to be embedded. Derived classes implementing this method should return -1 if
        they are not able to embed sets for that n and k.
        """
        if n<0 or k<0:
            return -1
        if n==0 or k==0:
            return 0
        if n==1:
            return k
        if k==1:
            return n
        return self.size_from_n_k_generic(n, k) # Derived classes should implement this method!
        
    def size_from_array(self, data: np.ndarray) -> int:
        """
        This function returns the number of reals that the embedding would contain if the set represented 
        by "data" were to be embedded. -1 is returned if data of the supplied type is not encodable by this embedder.
        """
        n,k = data.shape
        return self.size_from_n_k(n,k)

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """
        Derived classes should implement this method.
        This method should OPTIMALLY embed data for which n>=0 and k==1. We call this "kOne data".
        OPTIMALLY means that the embedding size must therefore be n.
        Implementations may assume (without checking) that data fed to it has the above type.
        It is likely that most implementations will implement this method either as:

            def embed_kOne(self, data: np.ndarray, debug=False) -> np.ndarray, Any:
                return MultisetEmbedder.embed_kOne_sorting(data) # Want piecewise linear!

        or as:

            def embed_kOne(self, data: np.ndarray, debug=False) -> np.ndarray, Any:
                return MultisetEmbedder.embed_kOne_polynomial(data) # Want Cinf!

        depending on whether they want piecewise linearity or differentiability.
        """
        raise NotImplementedError()

    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        """
        Derived classes should implement this method.
        This method should embed data for which n>1 and k>1. We call this "generic data".
        Implementations may assume (without checking) that data fed to it has the above type.
        """
        raise NotImplementedError()

    def size_from_n_k_generic(self, n: int, k:int) -> int:
        """
        Derived classes should implement this method.
        This method should report the embedding size for data for which n>1 and k>1. We call this "generic data".
        Implementations may assume (without checking) that data fed to it has the above type.
        """
        raise NotImplementedError()

    def test_me(self):
        _ = self.size_from_n_k_generic(2,2) # Check implementation exists
        _ = self.embed_generic(np.array([[1,2],[3,4]], dtype=np.float64)) # Check implementation exits
        _ = self.embed_kOne(np.array([[1,],[3,],], dtype=np.float64)) # Check implementation exits
        

    @staticmethod
    def embed_kOne_sorting(data: np.ndarray) -> np.ndarray:
        assert MultisetEmbedder.is_kOne_data(data)
        return np.sort(data.flatten())

    @staticmethod
    def embed_kOne_polynomial(data: np.ndarray) -> np.ndarray:
        assert MultisetEmbedder.is_kOne_data(data)
        embedding, shape_, metadata_ =  poly_list.embed(data.flatten())
        return embedding

    @staticmethod
    def is_kOne_data(data: np.ndarray) -> bool:
        n,k = data.shape
        return MultisetEmbedder.is_kOne_n_k(n,k)

    @staticmethod
    def is_kOne_n_k(n: int, k: int) -> bool:
        return n>=0 and k==1

    @staticmethod
    def is_generic_data(data: np.ndarray) -> bool:
        n,k = data.shape
        return MultisetEmbedder.is_generic_n_k(n,k)

    @staticmethod
    def is_generic_n_k(n: int, k: int) -> bool:
        return n>1 and k>1
