#!/usr/bin/env python3

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "echinocoder"))  # the library modules

import numpy as np

debug = False

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

n,k = set_of_vectors_to_embed.shape

test_sets_of_vectors = [ 
        set_of_vectors_to_embed,
        another_set_of_vectors_to_embed,
        tweaked_set_of_vectors_to_embed, 
    ]

import C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset as simplex1
import C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset as simplex2
import C0HomDeg1_conjectured_dotting_embedder_for_array_of_reals_as_multiset as conjectured_dotting

some_embedders = [ 
    simplex1.Embedder(),
    simplex2.Embedder(),
    conjectured_dotting.Embedder(n, k),
]


for arr in test_sets_of_vectors:
    print("")
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print("")
    
    print(f"\n==========================\nVectors to be encoded are:|")
    for vec in arr:
       print(f"{vec}, ") 

    concat = [ n, k]
    for embedder in some_embedders:
       embedding, (n_out, k_out), metadata = embedder.embed(arr)
       assert (n_out, k_out) == (n, k)
       if debug:
           print(f"\nembedding for embedder {embedder} is {embedding}\n")
       concat.extend( embedding )
    
    print(f"\n==========================\nconcatenated embedding is:")
    for number in concat:
        print(number)
    

    print("")
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print("")

