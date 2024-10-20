#!/usr/bin/env python3

import      Cinf_python_polynomial_encoder_for_list_of_reals_as_multiset as encoder_Cinf_py_li
import      Cinf_numpy_polynomial_encoder_for_array_of_reals_as_multiset as encoder_Cinf_np_ar
import          Cinf_sympy_bursar_encoder_for_array_of_reals_as_multiset as encoder_Cinf_sp_bur_ar
import      Cinf_sympy_evenBursar_encoder_for_array_of_reals_as_multiset as encoder_Cinf_sp_evenBur_ar
import       Cinf_numpy_polynomial_encoder_for_list_of_reals_as_multiset as encoder_Cinf_np_li
import                  C0_sorting_encoder_for_list_of_reals_as_multiset as encoder_C0_li
import C0_numpy_simplicialComplex_encoder_for_array_of_reals_as_multiset as encoder_C0_np_simplex

import data_sources
import numpy as np

fail_count = 0

def test_multiset_encoder(data, encoder=None, encoders=None, number_of_shuffled_copies=3, expected_encoding=None, relative_tolerance=0, absolute_tolerance=0):
    global fail_count
    print("ORIGINAL DATA is",data)
    shuffled_data = data.copy()
    if encoders is None:
        encoders = [ encoder ]
    for i in range(number_of_shuffled_copies):
        for encoder in encoders:
            encoding = encoder.encode(shuffled_data)
            #encoding_fails = expected_encoding is not None and not np.array_equal(np.asarray(encoding),np.asarray(expected_encoding))
            encoding_fails = expected_encoding is not None and not np.allclose(np.asarray(encoding),np.asarray(expected_encoding),rtol=relative_tolerance,atol=absolute_tolerance)
            if encoding_fails:
                fail_count += 1
                print("FAIL! Expected "+str(expected_encoding)+" ... ", end='')
            print("ENCH",encoder.__name__,"generates",encoding,"of length",len(encoding),"when encoding",shuffled_data,"with (m,n)=",shuffled_data.shape[::-1],". Expectation was "+str(expected_encoding))
            
        np.random.shuffle(shuffled_data) 
    print()

def self_test_realprojectivespace_encoder(encoder):
    for inp, out in encoder.unit_test_input_output_pairs:
       test_realprojectivespace_encoder(inp, encoder=encoder, expected_encoding=out)

def test_realprojectivespace_encoder(data, encoder=None, encoders=None, expected_encoding=None):
    global fail_count
    print("ORIGINAL DATA is",data)
    data_copy = data.copy()
    if encoders is None:
        encoders = [ encoder ]
    for sign in [+1,-1]:
        data_copy *= sign
        for encoder in encoders:
            encoding = encoder.encode(data_copy)
            encoding_fails = expected_encoding is not None and not np.array_equal(np.asarray(encoding),np.asarray(expected_encoding))
            if encoding_fails:
                fail_count += 1
                print("FAIL! Expected "+str(expected_encoding)+" ... ", end='')
            print("ENCH",encoder.__name__,"generates",encoding,"of length",len(encoding),"when encoding",data_copy,"having n=",len(data_copy),". Expectation was "+str(expected_encoding))
    print()

def make_randoms_reproducable():
    import random
    random.seed(1)
    np.random.seed(1)

def test_tools():


    global fail_count
    import tools


    data = np.asarray([3+4j, 4+2j])
    out_expected = [3, 4, 4, 2]
    out_real = tools.expand_complex_to_real_pairs(data)
    print("tools.expand_complex_to_real_pairs")
    print("data in ",data)
    print("out_expected ",out_expected)
    print("out_real     ",out_real)
    if not np.array_equal(out_expected,out_real):
        fail_count += 1
        print("FAIL in expand_complex_to_real_pairs")
    print()

    for data, expected in ( ([1,2,3,4,5,], [1+2j, 3+4j, 5,]),  ([1,2,3,4,], [1+2j, 3+4j, ]),  ):
       data = np.asarray(data)
       expected = np.asarray(expected)
       out = tools.real_pairs_to_complex_zip(data)
       print("tools.real_pairs_to_complex_zip")
       print("data in ",data)
       print("out_expected ",expected)
       print("out_real     ",out)
       if not np.array_equal(expected,out):
           fail_count += 1
           print("FAIL in real_pairs_to_complex_zip")
       print()

    


def test_various_encoders():
    make_randoms_reproducable()

    test_multiset_encoder(
       data=np.asarray([9,-4,21,-8,5]),
       encoder=encoder_C0_li,
       expected_encoding = [-8,-4,5,9,21]
    )

    test_multiset_encoder(
       data=data_sources.random_real_linear_data(n=4),
       encoder=encoder_C0_li,
    )

    test_multiset_encoder(
       data=data_sources.random_complex_linear_data(n=4),
       encoders=[ encoder_Cinf_np_li, encoder_Cinf_py_li, ],
    )

    test_multiset_encoder(
       data=data_sources.random_real_linear_data(n=4),
       encoders=[ encoder_Cinf_np_li, encoder_Cinf_py_li, ],
    )

    test_multiset_encoder(
       data=data_sources.random_real_array_data(mn=(1,4)),
       encoders=[ encoder_Cinf_np_ar, ],
    )

    test_multiset_encoder(
       data=data_sources.random_real_array_data(mn=(2,3)),
       encoders=[ encoder_Cinf_np_ar, ],
    )

    test_multiset_encoder(
       data=data_sources.random_real_array_data(mn=(3,3)),
       encoders=[ encoder_Cinf_np_ar, ],
    )

    test_multiset_encoder(
       data=data_sources.random_real_array_data(mn=(4,3)),
       encoders=[ encoder_Cinf_np_ar, encoder_Cinf_sp_bur_ar, ],
    )

    test_multiset_encoder(
       data=np.array(((-7,-8,-1,-9),(9,-7,-6,5),(-9,4,9,-7))),
       encoders=[ encoder_Cinf_sp_bur_ar, ],
       expected_encoding = [-11, 2, -11, -7, -17, 62, 55, -202, 110, 120, -81, 315, -748, 58, 1289, -1497, 457, 1139, -1460, -45, 567],
    )

    test_multiset_encoder(
       data=np.array(((8,-1,-4,3),(-8,-5,9,7),(8,2,7,-7))),
       encoders=[ encoder_Cinf_np_ar, ],
       expected_encoding = [   8,   -4,  -57,  -80, -488, -394,    8,   12,  -63,  144, -952,  636, 8,    3,  -15,  112, -456,  851,   -4,   12,   -6,  -21,    5,  309, -4,    3,   42,   40, -186,   68,   12,    3,   48,   34, -406,  392], 
    )

    test_multiset_encoder(
       data=np.array(((3,1,4),(2,2,5))),
       encoders=[ encoder_Cinf_sp_bur_ar, ],
       expected_encoding =  [9, 3, 5, 20, 13, 25, 8, 6],
    )

    test_multiset_encoder(
       data=np.array(((3,1,4),(2,2,5))),
       encoders=[ encoder_Cinf_sp_evenBur_ar, ],
       number_of_shuffled_copies=10,
       expected_encoding = [9, 20, 3, 13, 2, 5, 23, 8, 6],
    )

    import Cinf_numpy_regular_encoder_for_list_of_realsOrComplex_as_realOrComplexprojectivespace
    self_test_realprojectivespace_encoder(Cinf_numpy_regular_encoder_for_list_of_realsOrComplex_as_realOrComplexprojectivespace)

    import Cinf_numpy_complexPacked_encoder_for_list_of_reals_as_realprojectivespace
    self_test_realprojectivespace_encoder(Cinf_numpy_complexPacked_encoder_for_list_of_reals_as_realprojectivespace)

    test_multiset_encoder(
       data=np.array(((1,2),(1,0),(5,2))),
       encoders=[ encoder_C0_np_simplex, ],
       number_of_shuffled_copies=10,
       expected_encoding = [1.74545455e+00, 1.13000000e+01, 8.73272727e+01, 7.28500000e+02, 6.30183636e+03, 5.54913000e+04, 4.93020782e+05],
       relative_tolerance=1e-8,
    )

def test_everything():
    test_tools()
    test_various_encoders()
    print(str(fail_count)+" failures")
    import tuple_rank
    tuple_rank.unit_test_tuple_rank()

test_everything()
