#!/usr/bin/env python3

"""
Running this script should generate output resembling the following:

The set of vectors 
[[ 8 -1 -4  3]
 [-8 -5  9  7]
 [ 8  2  7 -7]]
embeds to 
(array([   8.,   -4.,  -57.,  -80., -488., -394.,    8.,   12.,  -63.,
        144., -952.,  636.,    8.,    3.,  -15.,  112., -456.,  851.,
         -4.,   12.,   -6.,  -21.,    5.,  309.,   -4.,    3.,   42.,
         40., -186.,   68.,   12.,    3.,   48.,   34., -406.,  392.]), (3, 4), None).
We expected
[8, -4, -57, -80, -488, -394, 8, 12, -63, 144, -952, 636, 8, 3, -15, 112, -456, 851, -4, 12, -6, -21, 5, 309, -4, 3, 42, 40, -186, 68, 12, 3, 48, 34, -406, 392]


Another set of vectors 
[[ 8  2  7 -7]
 [ 8 -1 -4  3]
 [-8 -5  9  7]]
differs from the first by a permutation of the vectors.
Nonetheless, it embeds to 
(array([   8.,   -4.,  -57.,  -80., -488., -394.,    8.,   12.,  -63.,
        144., -952.,  636.,    8.,    3.,  -15.,  112., -456.,  851.,
         -4.,   12.,   -6.,  -21.,    5.,  309.,   -4.,    3.,   42.,
         40., -186.,   68.,   12.,    3.,   48.,   34., -406.,  392.]), (3, 4), None)
too.


If we change one of the vectors by a small amount, the embedding changes by only a small amount. E.g.
[[ 8.    -1.    -4.     3.   ]
 [-8.    -5.     9.     7.   ]
 [ 8.     2.     7.    -7.001]]
embeds to:
(array([   8.   ,   -4.   ,  -57.   ,  -80.   , -488.   , -394.   ,
          8.   ,   12.   ,  -63.   ,  144.   , -952.   ,  636.   ,
          8.   ,    2.999,  -14.99 ,  112.   , -455.968,  851.085,
         -4.   ,   12.   ,   -6.   ,  -21.   ,    5.   ,  309.   ,
         -4.   ,    2.999,   42.01 ,   40.006, -186.022,   68.016,
         12.   ,    2.999,   48.01 ,   33.995, -406.001,  392.057]), (3, 4), None)

Now let's try some other embedders

Using the simplex1 embedder ... 
... we find that the ouput is (array([ 4.71175355,  5.97716044,  2.27870903,  6.31338771,  4.50981973,
       10.55305146,  7.39127809,  9.11676996,  5.40387331,  5.29204317,
       12.18165909,  7.91533274,  6.34160765,  6.40703423,  6.30875449,
        8.19956912,  3.62518034,  7.35944864,  4.40777783, 12.99478808,
       12.77917328, 10.71712844,  7.73758228,  9.        , -8.        ]), (3, 4), None).

Using the simplex2 embedder ... 
... we find that the ouput is (array([-8.        , -5.        , -4.        , -7.        , 28.33984331,
       36.83771553, 18.1290216 , 41.57434862, 24.52992354, 19.80764839,
       29.21590987, 25.43954949, 22.59115854, 17.53657363, 15.48742714,
       35.84843865, 23.83542322, 32.88810865, 13.79639461, 17.94517704,
       24.16444411]), (3, 4), {'ascending_data': array([[-8, -5, -4, -7],
       [ 8, -1,  7,  3],
       [ 8,  2,  9,  7]]), 'input_data': array([[ 8, -1, -4,  3],
       [-8, -5,  9,  7],
       [ 8,  2,  7, -7]])}).

Using the conjectured_dotting embedder ... 
... we find that the ouput is (array([ -8.        ,   8.        ,   8.        ,  -5.        ,
        -1.        ,   2.        ,  -4.        ,   7.        ,
         9.        ,  -7.        ,   3.        ,   7.        ,
        -1.10904362,   4.49028978,   6.15278722,  -7.02170733,
        -1.0637313 ,  20.84294686, -12.81292792,  -1.74738463,
         6.63680093]), (3, 4), None).

"""

import numpy as np


set_of_vectors_to_embed = np.array([
      (8,-1,-4,3),
      (-8,-5,9,7),
      (8,2,7,-7)])
another_set_of_vectors_to_embed = np.array([
      (8,2,7,-7),
      (8,-1,-4,3),
      (-8,-5,9,7)])
tweaked_set_of_vectors_to_embed = np.array([
      (8,-1,-4,3),
      (-8,-5,9,7),
      (8,2,7,-7.001)])

import Cinf_numpy_polynomial_encoder_for_array_of_reals_as_multiset as Cinf_np_ar
embedder = Cinf_np_ar.Encoder()

embedding = embedder.encode(set_of_vectors_to_embed)
another_embedding = embedder.encode(another_set_of_vectors_to_embed)

expected_embedding = [   8,   -4,  -57,  -80, -488, -394,    8,   12,  -63,  144, -952,  636, 8,    3,  -15,  112, -456,  851,   -4,   12,   -6,  -21,    5,  309, -4,    3,   42,   40, -186,   68,   12,    3,   48,   34, -406,  392]

print(f"The set of vectors \n{set_of_vectors_to_embed}\nembeds to \n{embedding}.\nWe expected\n{expected_embedding}\n\n")
print(f"Another set of vectors \n{another_set_of_vectors_to_embed}\ndiffers from the first by a permutation of the vectors.\nNonetheless, it embeds to \n{another_embedding}\ntoo.\n\n")
print(f"If we change one of the vectors by a small amount, the embedding changes by only a small amount. E.g.\n{tweaked_set_of_vectors_to_embed}\nembeds to:")
print(embedder.encode(tweaked_set_of_vectors_to_embed))

print("")
print("Now let's try some other embedders")
print("")

import C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset as simplex1
import C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset as simplex2
import C0HomDeg1_conjectured_dotting_embedder_for_array_of_reals_as_multiset as conjectured_dotting

embedders_to_try = (
                     ("simplex1", simplex1.Embedder()),
                     ("simplex2", simplex2.Embedder()),
                     ("conjectured_dotting", conjectured_dotting.Embedder(n=3, k=4)), # This embedder needs to specialise on particular n and k
                   )

for name, embedder in embedders_to_try:
    print(f"Using the {name} embedder ... ")
    embedding = embedder.embed(set_of_vectors_to_embed)
    print(f"... we find that the ouput is {embedding}.")
    print()


   

