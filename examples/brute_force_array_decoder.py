#!/usr/bin/env python 

"""
The purpose of this script is to provide a method that will attempt a brute force decoding of any array encoder which outputs "input_data" within its metadata assuming that all it has to do is re-order ascending_data(input_data).

Note that this is not a true decoder, as it does not confirm that the encoding could generate the ascending_data. However, any encoder could pacakge ascending_data alongside any other encoding without increasing its order beyond (nk) if already (nk) .... 
"""

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "echinocoder"))  # the library modules

import numpy as np
import C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset as simplex2
import tools

def do_test_on(embedder, data):
    encoding, (n,k), metadata = embedder.embed(data)
    
    print(f"data = {data}")
    print(f"ecoding = {encoding}")
    print(f"metadata = {metadata}")
    
    def canonical_form(data):
        return tools.sort_np_array_rows_lexicographically(data)
    
    # See if we can attempt brute force decode:
    if "ascending_data" in metadata:
        # yes, it seems we can try
        ascending_data = metadata["ascending_data"]
        print(f"ascending data = {ascending_data}")
    
        if "input_data" in metadata:
            assert np.equal(ascending_data, tools.ascending_data(metadata["input_data"])).all()
    
        bad_matches = set()
        good_matches = set()
        bad_cover_factor = 0
        good_cover_factor = 0

        canonical = canonical_form(data)

        next_print = 1
        for counter, hypothesis in enumerate(tools.permute_columns_except_first(ascending_data)):
            #print("Might be:")
            #print(hypothesis)
            encoding_of_hypothesis, _, _ = embedder.embed(hypothesis)
            if counter+1 == next_print:
                next_print = int((next_print+1)*1.1)
                print(".", end="", flush=True)
    
            exact_match = np.array_equal(encoding_of_hypothesis, encoding)
            approx_match = np.allclose(encoding_of_hypothesis, encoding)
         
            # Currenrly use approx match
            if approx_match:
                print("\n\nWe found an encoding match between")
                print(data)
                print("and")
                print(hypothesis)
                canonical_hypothesis = canonical_form(hypothesis)
    
                ch = tuple(tuple(i) for i in canonical_hypothesis.tolist())
    
                if (canonical == canonical_hypothesis).all():
                    good_matches.add(ch)
                    print("This is a GOOD DESIRED match.")
                    good_cover_factor += 1
                else:
                    bad_matches.add(ch)
                    print("This match is BAD.")
                    bad_cover_factor += 1
                print(f"So far: #GOOD={len(good_matches)} and #BAD={len(bad_matches)}.") 
    
        print(f"At end: GOOD={good_matches} and BAD={bad_matches}.") 
        print(f"At end: #GOOD={len(good_matches)} and #BAD={len(bad_matches)}.") 
        print(f"At end: #GOOD_COVER_FACTOR={good_cover_factor}")
        print(f"At end: #BAD_COVER_FACTOR={bad_cover_factor}")
    
        assert len(bad_matches) == 0 
        assert len(good_matches) == 1

        print("Passed test\n=============\n")




def run_the_decoder_a_few_times():
    embedder = simplex2.Embedder()
    
    
    do_test_on(embedder, np.array([[1,2],
                                   [1,3],
                                   [1,4],
                                   [1,5],
                                   [1,6]]))
    
    do_test_on(embedder, np.array([[1,2,3],
                                   [4,-5,6],
                                   [-7,8,9]]))
    
    do_test_on(embedder, np.array([[1,2,3],
                                   [4,-5,6],
                                   [4,-5,6],
                                   [-7,8,9]]))
    
    
    do_test_on(embedder, np.array([[1,2,3],
                                   [4,-5,6],
                                   [4.1,-5.1,6.1],
                                   [4.2,-5.2,6.2],
                                   [-7,8,9]]))
    
    do_test_on(embedder, np.array([[1,2,3],
                                   [4,-5,6],
                                   [4,-5,6],
                                   [4,-5,61],
                                   [4,-5,61],
                                   [4,-5,61],
                                   [4,-5,61],
                                   [4,-5,61],
                                   [-7,8,9]]))

if __name__ == "__main__":
    run_the_decoder_a_few_times()
