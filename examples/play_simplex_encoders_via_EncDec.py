#!/usr/bin/env python

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "echinocoder"))  # the library modules

import EncDec
import numpy as np
from itertools import pairwise
import sys
import tools

def pretty_print_lin_comb(lin_comb: EncDec.LinComb):
    for coeff, basis_elt in zip(lin_comb.coeffs, lin_comb.basis_vecs):
        print(" + ", float(coeff), "*", tools.numpy_array_of_frac_to_str(basis_elt), f" (index sum {np.sum(basis_elt)})")

def print_first_part_of_simplex_encoding(set_array : np.array,
                                         simplex_method = EncDec.simplex_1_preprocess_steps,
                                         preserve_scale_in_step_1 = False,
                                         preserve_scale_in_step_2 = False,
                                         canonicalise = True,
                                         debug = False):

    lin_comb_2_second_diffs, offsets = simplex_method(set_array,
                                                      preserve_scale_in_step_1 = preserve_scale_in_step_1,
                                                      preserve_scale_in_step_2 = preserve_scale_in_step_2,
                                                      canonicalise = canonicalise,
                                                      use_assertions = True,
                                                      debug = debug)

    lin_comb_3 = lin_comb_2_second_diffs + offsets

    print(f"=======================\n{simplex_method} encoded")
    print(EncDec.numpy_array_of_frac_to_str(set_array))
    print(f"{'with' if canonicalise else 'without'} canonicalisation to second diff(s) ")
    pretty_print_lin_comb(lin_comb_2_second_diffs)
    print("and offset(s)")
    pretty_print_lin_comb(EncDec.LinComb(offsets))
    print("========================")
    if len(lin_comb_2_second_diffs)>=1:
        print("The first diff basis vector in the above is")
        print(EncDec.numpy_array_of_frac_to_str(tmp:=lin_comb_2_second_diffs.basis_vecs[0]), " with one-norm ", np.sum(tmp))
    if len(lin_comb_2_second_diffs)>=2:
        print("and the differences between the subsequent (non-offset) basis vectors are:")
        [ print(EncDec.numpy_array_of_frac_to_str(tmp:=b-a), " with one-norm ", np.sum(tmp)) for a,b in list(pairwise( lin_comb_2_second_diffs.basis_vecs )) ]





def loc(the_list, query):
    for i, entry in enumerate(the_list):
        if entry == query:
            return i
    return None

if __name__ == "__main__":

    # Defaults:
    arr = np.array([[2,3,4], [4,7,1], [3,-2,1], [9,8,2]])
    preserve_scale_in_step_1 = False
    preserve_scale_in_step_2 = False
    simplex_method = EncDec.simplex_1_preprocess_steps

    # Override defaults:
    if (l := loc(sys.argv, "array")) is not None and l+1 < len(sys.argv): 
        arr = np.array(eval(sys.argv[l+1]))
    
    low=9
    high=30

    n, k = arr.shape

    if (l := loc(sys.argv, "n")) is not None and l+1 < len(sys.argv): 
        n = int(sys.argv[l+1])
    if (l := loc(sys.argv, "k")) is not None and l+1 < len(sys.argv): 
        k = int(sys.argv[l+1])
    if (l := loc(sys.argv, "low")) is not None and l+1 < len(sys.argv): 
        low = int(sys.argv[l+1])
    if (l := loc(sys.argv, "high")) is not None and l+1 < len(sys.argv): 
        high = int(sys.argv[l+1])
    if (l := loc(sys.argv, "method")) is not None and l+1 < len(sys.argv): 
        tmp = int(sys.argv[l+1])
        if tmp==1:
            simplex_method = EncDec.simplex_1_preprocess_steps
        elif tmp==2:
            simplex_method = EncDec.simplex_2_preprocess_steps
        else:
            raise ValueError(f"Don't know method {tmp}. Only know methods 1 and 2.")
            
    if (l := loc(sys.argv, "scales")) is not None and l+2 < len(sys.argv):
        preserve_scale_in_step_1 = bool(int(sys.argv[l+1]))
        preserve_scale_in_step_2 = bool(int(sys.argv[l+2]))
    
    randomise_array = "random" in sys.argv
    if randomise_array:
         #print (f"n is {n} and k is {k}")
         arr = np.random.randint(low=low, high=high, size=(n,k))
   
    print()
    for canonicalise in (False, True):
        print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        print_first_part_of_simplex_encoding(arr, simplex_method=simplex_method,  canonicalise=canonicalise, preserve_scale_in_step_1=preserve_scale_in_step_1, preserve_scale_in_step_2=preserve_scale_in_step_2)
        print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        print()
    

