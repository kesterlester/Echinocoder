# Embed list of n real numbers treated as an element of a realprojectivespace.
#
# This particular implentation explots a complex trick which (when n is even and greater than 0) allows the embedding save one real number over that required by trhe underlying regular implementation.
#
# I.e. the method embeds to:
#
#       0 reals when n is 0
#       2n-2 reals when n is even and n>0     (A SAVING OF ONE REAL!!)
#       2n-1 reals when n is odd.
#
# E.g. this method can embed an ordered list like this
#
#        x = [3,4,-2,5]
#
# which is deemed to be equivalent to its negation:
#
#        -x = [-3,-4,2,-5].
#
# Here is a table of the first few embedding orders (first column is n, second column is the order):
#
#  0  0
#  1  1
#  2  2
#  3  5
#  4  6
#  5  9
#  6 10


import numpy as np
import tools
import Cinf_numpy_regular_embedder_for_list_of_realsOrComplex_as_realOrComplexprojectivespace as underlying_embedder

def embed(data, debug=False):
    # data=np.asarray(data)
    complex_data = tools.real_pairs_to_complex_zip(data)
    complex_embedding = underlying_embedder.embed(complex_data)
    real_embedding = tools.expand_complex_to_real_pairs(complex_embedding)
    # Now trim off last entry but ONLY IF we know it will be zero by construction -- this is when the data has an odd length
    if len(data)%2 == 1:
        real_embedding = real_embedding[:-1]
    return real_embedding

unit_test_input_output_pairs = [
                    ( np.asarray([  ]),       [ ] ),
                    ( np.asarray([-3]),       [9] ), 
                    ( np.asarray([3,4,-2,5]), [ -7, 24, -52, 14, -21, -20] ),
                    ( np.asarray([3,4,-2  ]), [ -7, 24, -12, -16, 4] ),
                    ]

def main():
    for data in (np.asarray([3,4,-2,5]), np.asarray([3,4,-2])):
        print()
        complex_data = tools.real_pairs_to_complex_zip(data)
        complex_embedding = underlying_embedder.embed(complex_data)
        real_embedding = tools.expand_complex_to_real_pairs(complex_embedding)
        # Now trim off last entry but ONLY IF we know it will be zero by construction -- this is when the data has an odd length
        if len(data)%2 == 1:
            real_embedding = real_embedding[:-1]
    
        print("data ",data)
        print("complex_data ",complex_data)
        print("complex_embedding",complex_embedding)
        print("real_embedding ",real_embedding)
        # Should print
        """
        input  [ 3  4 -2  5]
        reals  [ 3 -2]
        images  [4 5]
        data  [ 3  4 -2  5]
        complex_data  [ 3.+4.j -2.+5.j]
        complex_embedding [ -7.+24.j -52.+14.j -21.-20.j]
        real_embedding  [ -7.  24. -52.  14. -21. -20.]
        
        input  [ 3  4 -2]
        reals  [ 3 -2]
        images  [4 0]
        data  [ 3  4 -2]
        complex_data  [ 3.+4.j -2.+0.j]
        complex_embedding [ -7.+24.j -12.-16.j   4. +0.j]
        real_embedding  [ -7.  24. -12. -16.   4.]
        """


# Just for testing/debug:
if __name__ == "__main__":
    main()

