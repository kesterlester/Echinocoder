# Embed array of real numbers treated as multiset.
# E.g. this method can embed things like:
#
#        [[3,4],[4,2],[-2,1]]
#
# when used to represent this multiset {{ }}:
#
#        {{ [3,4],[4,2],[-2,1] }}
#
# total_embedding_size = m*(m-1)*n
#
# where n is the number of vectors in the multiset and m is their dimension. E.g in the example above m=2 and n=3.

from MultisetEmbedder import MultisetEmbedder
import numpy as np
import Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset
import itertools
import tools

class Embedder(MultisetEmbedder):
    def embed(self, data: np.ndarray, debug=False) -> np.ndarray:
    
        n,m = data.shape
        #print(f"Size (m,n)=({m},{n})")
    
        # We need to do different things depending on the shape of the input:
    
        if m==1:
            # The "vectors" are 1-long, so reduce the thing to a list.
            real_embedding = Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.embed(data[:,0])
            assert len(real_embedding) == self.size_from_n_k(n,m)
            return real_embedding
    
        if m==2:
            # The "vectors" are 2-long, so turn them into complex numbers, embed, and expand:
            complex_list = data[:,0]+complex(0,1)*data[:,1]
            complex_embedding = Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.embed(complex_list) # Yes, this is explopiting a secret private feature ...
            real_embedding = tools.expand_complex_to_real_pairs(complex_embedding)
            assert len(real_embedding) == self.size_from_n_k(n,m)
            return real_embedding
    
        # The "vectors" are 3-or-more-long, so need to take all m-choose-2 pairings of elements in them, and send them out.
        # We can use the m==2 case code above for that.
        number_of_pairs = m*(m-1)//2
        # pair_embedding_size = 2*n 
        # total_embedding_size = m*(m-1)*n
        current_embedding_index = 0 
        for row1_index, row2_index in itertools.combinations(range(m), 2):
            pair_of_rows = data[:,[row1_index, row2_index]]
            #print("About to request coding for",pair_of_rows)
            embedding_for_pair = self.embed(pair_of_rows)
            if current_embedding_index == 0:
                # Allocate array of appropriate type for all the pairs that will come later:
                # Get embedding size.  This SHOULD be 2*n but safer to get from data for future proofing
                pair_embedding_size = len(embedding_for_pair)
                #print("pes",pair_embedding_size,"nop",number_of_pairs)
                total_embedding_size = pair_embedding_size * number_of_pairs
                real_embedding = np.empty(total_embedding_size, dtype = embedding_for_pair.dtype)
            real_embedding[current_embedding_index:current_embedding_index+pair_embedding_size] = embedding_for_pair
            current_embedding_index += pair_embedding_size
    
        assert len(real_embedding) == n*m*(m-1)
        assert len(real_embedding) == self.size_from_n_k(n,m)
        return real_embedding

    def size_from_n_k(self, n: int, k: int) -> int:
        if n==0 or k==0: return 0
        return n*k*(k-1) if k>1 else n

