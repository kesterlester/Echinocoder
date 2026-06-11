#!/usr/bin/env python3

import Historical.C0_simplicialComplex_embedder_1_for_array_of_reals_as_multiset as C0_np_simplex_historical
embedder_C0_np_simplex_historical = C0_np_simplex_historical.Embedder()

import C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset as C0HomDeg1_simplex1
embedder_C0HomDeg1_simplex1 = C0HomDeg1_simplex1.Embedder()
import C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset as C0HomDeg1_simplex2
embedder_C0HomDeg1_simplex2 = C0HomDeg1_simplex2.Embedder()
import C0HomDeg1_conjectured_dotting_embedder_for_array_of_reals_as_multiset as conjectured_dotting
embedder_conjectured_dotting = conjectured_dotting.Embedder(n=2, k=2)
import  Cinf_numpy_polynomial_encoder_for_array_of_reals_as_multiset as Cinf_np_ar
embedder_Cinf_np_ar = Cinf_np_ar.Encoder()
import      Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset as Cinf_sp_bur_ar
embedder_Cinf_sp_bur_ar = Cinf_sp_bur_ar.Embedder()
import  Cinf_sympy_evenBursar_embedder_for_array_of_reals_as_multiset as Cinf_sp_evenBur_ar
embedder_Cinf_sp_evenBur_ar = Cinf_sp_evenBur_ar.Embedder()
import Cinf_hybrid_encoder_for_array_of_reals_as_multiset as hybrid
embedder_hybrid = hybrid.Encoder()

import  Cinf_python_polynomial_embedder_for_list_of_reals_as_multiset as embedder_Cinf_py_li
import   Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset as embedder_Cinf_np_li
import              C0_sorting_embedder_for_list_of_reals_as_multiset as embedder_C0_li
import data_sources
import numpy as np
from tools import __line__, permute_columns_except_first
import tuple_ize

def self_test_realprojectivespace_embedder(embedder):
    for inp, out in embedder.unit_test_input_output_pairs:
        external_test_realprojectivespace_embedder(inp, embedder=embedder, expected_embedding=out)

def external_test_realprojectivespace_embedder(data, embedder=None, embedders=None, expected_embedding=None):
    print("ORIGINAL DATA is",data)
    data_copy = data.copy()
    if embedders is None:
        embedders = [ embedder ]
    for sign in [+1,-1]:
        data_copy *= sign
        for embedder in embedders:
            embedding = embedder.embed(data_copy)
            print("ENCH",embedder.__name__,"generates",embedding,"of length",len(embedding),"when embedding",data_copy,"having n=",len(data_copy),". Expectation was "+str(expected_embedding))
            assert expected_embedding is not None and np.array_equal(np.asarray(embedding),np.asarray(expected_embedding))
    print()

def make_randoms_reproducable():
    import random
    random.seed(1)
    np.random.seed(1)

def test_tools():
    import tools

    data = np.asarray([3+4j, 4+2j])
    out_expected = [3, 4, 4, 2]
    out_real = tools.expand_complex_to_real_pairs(data)
    print("tools.expand_complex_to_real_pairs")
    print("data in ",data)
    print("out_expected ",out_expected)
    print("out_real     ",out_real)
    assert np.array_equal(out_expected,out_real)
    print()

    for data, expected in ( ([1,2,3,4,5,], [1+2j, 3+4j, 5,]),  ([1,2,3,4,], [1+2j, 3+4j, ]),  ):
        data = np.asarray(data)
        expected = np.asarray(expected)
        out = tools.real_pairs_to_complex_zip(data)
        print("tools.real_pairs_to_complex_zip")
        print("data in ",data)
        print("out_expected ",expected)
        print("out_real     ",out)
        assert np.array_equal(expected,out)
        print()

    perm = [3,0,1,2]
    inverse_perm_calculated = tools.invert_perm(perm)
    inverse_perm_expected = [1,2,3,0]

    assert np.array_equal(inverse_perm_calculated, inverse_perm_expected)
        
    perm = [3,5,0,4,1,2]
    inverse_perm_calculated = tools.invert_perm(perm)
    inverse_perm_expected = [2, 4, 5, 0, 3, 1]

    assert np.array_equal(inverse_perm_calculated, inverse_perm_expected)


    #### CHECK ARRAY ROW SORT ####
    inp = np.array([[1, 0, 2],
                    [0, 5, 2],
                    [3, 0, 8]])

    expected = np.array([[0, 1, 2],
                         [0, 2, 5],
                         [0, 3, 8]])
    got = tools.sort_each_np_array_row(inp)
    assert (got == expected).all()

    #### CHECK ARRAY COLUMN SORT ####
    inp = np.array([[1, 0, 3],
                    [0, 5, 2],
                    [3, 0, 8]])
    expected = np.array([[0, 0, 2],
                         [1, 0, 3],
                         [3, 5, 8]])
    got = tools.sort_each_np_array_column(inp)
    assert (got == expected).all()


    ### CHECK COLUMN PERMUTER ###
    arr = np.array([[3, 4, 5], 
                    [6, 7, 8]])
    perm_gen = permute_columns_except_first(arr)
    expected_answers = [
                         np.array([[3, 4, 5],
                                   [6, 7, 8]]),
                         np.array([[3, 4, 8],
                                   [6, 7, 5]]),
                         np.array([[3, 7, 5],
                                   [6, 4, 8]]),
                         np.array([[3, 7, 8],
                                   [6, 4, 5]]),
                       ]
    i=-1
    for perm in perm_gen:
       i+=1
       print(perm, '\n')
       assert (expected_answers[i] == perm).all()
    assert i+1 == 4


def tost_multiset_embedder(data, embedder=None, embedders=None, number_of_shuffled_copies=3, expected_embedding=None, relative_tolerance=0, absolute_tolerance=0):
    exact = relative_tolerance == 0 and absolute_tolerance == 0
    print("ORIGINAL DATA is",data)
    shuffled_data = data.copy()
    if embedders is None:
        embedders = [ embedder ]
    first_embedding = dict()
    for i in range(number_of_shuffled_copies):
        for embedder in embedders:
            print(f"embedder is {embedder}")
            
            try_to_encode = True
            if hasattr(embedder, "size_from_array"):
                siz = embedder.size_from_array(shuffled_data)
                if siz==-1:
                    try_to_encode = False

            if not try_to_encode:
                print(f"Skipping test of data on {embedder} as shape ({shuffled_data.shape}) is wrong.")
            else:
                print(f"Trying embedder {embedder} on {shuffled_data}.")
                if hasattr(embedder, 'encode') and callable(embedder.encode):
                    embedding, size_, metadata_ = embedder.encode(shuffled_data, debug=True)
                else:
                    embedding, size_, metadata_ = embedder.embed(shuffled_data, debug=True)
                print(f"got {embedding}, {size_}, {metadata_}")
                #embedding_fails = expected_embedding is not None and not np.array_equal(np.asarray(embedding),np.asarray(expected_embedding))
                # Check subsequent embeddings are same as first embedding. I.e. check for permutation invariance.
                if i==0:
                    first_embedding[embedder]=embedding
                else:
                    #print("MOOOCOWFIRST",embedder, first_embedding[embedder])
                    #print("MOOOCOW__NOW",embedder, embedding)
                    if exact:
                        np.testing.assert_equal(embedding, first_embedding[embedder], strict=True)
                    else:
                        np.testing.assert_allclose(np.array(embedding, dtype=float), np.array(first_embedding[embedder], dtype=float), atol=absolute_tolerance, rtol=relative_tolerance, strict=True, equal_nan=False)

                # Also checks embedding against expected, if given
                if expected_embedding is not None:
                    #print("MOOO1",embedding)
                    #print("MOOO2",expected_embedding)
                    if exact:
                        np.testing.assert_equal(embedding, expected_embedding, strict=True)
                    else:
                        np.testing.assert_allclose(np.array(embedding, dtype=float), np.array(expected_embedding, dtype=float), atol=absolute_tolerance, rtol=relative_tolerance, strict=True, equal_nan=False)

        np.random.shuffle(shuffled_data) 
    print()

def test_various_embedders():

    # with lists as inputs:
    # self.assertEqual(ell( [(1,5),(2,42),(3,100)], 101),
    #                  ell( [(7,5),(1,42),(3,100)], 101))
    # ( (172) maps 1->7 and 2->1 in S(8) )


    make_randoms_reproducable()
    all_ar_embedders=[
        embedder_C0HomDeg1_simplex1,
        embedder_C0HomDeg1_simplex2,
        embedder_conjectured_dotting,
        embedder_Cinf_np_ar,
        embedder_Cinf_sp_bur_ar,
        embedder_hybrid,
        embedder_C0_np_simplex_historical,
        ]

    print(__file__, __line__)
    tost_multiset_embedder(
       data=np.asarray([9,-4,21,-8,5]),
       embedder=embedder_C0_li,
       expected_embedding = [-8,-4,5,9,21],
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_linear_data(n=4),
       embedder=embedder_C0_li,
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_complex_linear_data(n=4),
       embedders=[ embedder_Cinf_np_li, embedder_Cinf_py_li, ],
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_linear_data(n=4),
       embedders=[ embedder_Cinf_np_li, embedder_Cinf_py_li, ],
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_array_data(mn=(4,0)),
       embedders=all_ar_embedders,
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_array_data(mn=(0,4)),
       embedders=all_ar_embedders,
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_array_data(mn=(4,1)),
       embedders=all_ar_embedders,
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_array_data(mn=(1,4)),
       embedders=all_ar_embedders,
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_array_data(mn=(2,3)),
       embedders=all_ar_embedders,
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_array_data(mn=(3,3)),
       embedders=[ 
       embedder_Cinf_np_ar, 
       embedder_Cinf_sp_bur_ar, 
       embedder_C0HomDeg1_simplex1, 
       embedder_C0HomDeg1_simplex2, 
       embedder_C0_np_simplex_historical,
       ],
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=data_sources.random_real_array_data(mn=(4,3)),
       embedders=[ 
         embedder_Cinf_np_ar, 
         embedder_Cinf_sp_bur_ar, 
         embedder_C0HomDeg1_simplex1, 
         embedder_C0HomDeg1_simplex2, 
         embedder_C0_np_simplex_historical,
         ],
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=np.array(((-7,-8,-1,-9),(9,-7,-6,5),(-9,4,9,-7))),
       embedders=[ embedder_Cinf_sp_bur_ar, ],
       expected_embedding = [-11, 2, -11, -7, -17, 62, 55, -202, 110, 120, -81, 315, -748, 58, 1289, -1497, 457, 1139, -1460, -45, 567],
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=np.array(((8,-1,-4,3),(-8,-5,9,7),(8,2,7,-7))),
       embedders=[ embedder_Cinf_np_ar, ],
       expected_embedding = [   8.0,   -4,  -57,  -80, -488, -394,    8,   12,  -63,  144, -952,  636, 8,    3,  -15,  112, -456,  851,   -4,   12,   -6,  -21,    5,  309, -4,    3,   42,   40, -186,   68,   12,    3,   48,   34, -406,  392], 
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=np.array(((3,1,4),(2,2,5))),
       embedders=[ embedder_Cinf_sp_bur_ar, ],
       expected_embedding =  [9, 3, 5, 20, 13, 25, 8, 6],
    )

    print(__file__, __line__)
    tost_multiset_embedder(
       data=np.array(((3,1,4),(2,2,5))),
       embedders=[ embedder_Cinf_sp_evenBur_ar, ],
       number_of_shuffled_copies=10,
       expected_embedding = [9, 20, 3, 13, 2, 5, 23, 8, 6],
    )

    print(__file__, __line__)
    import Cinf_numpy_regular_embedder_for_list_of_realsOrComplex_as_realOrComplexprojectivespace
    self_test_realprojectivespace_embedder(Cinf_numpy_regular_embedder_for_list_of_realsOrComplex_as_realOrComplexprojectivespace)

    print(__file__, __line__)
    import Cinf_numpy_complexPacked_embedder_for_list_of_reals_as_realprojectivespace
    self_test_realprojectivespace_embedder(Cinf_numpy_complexPacked_embedder_for_list_of_reals_as_realprojectivespace)

    print(__file__, __line__)
    tost_multiset_embedder(
       data=np.array(((1,2),(1,0),(5,2))),
       embedders=[ 
         embedder_C0_np_simplex_historical,
       ],
       number_of_shuffled_copies=100,
       expected_embedding = [ 0.41666667, -0.28108724,  0.26858269, -0.34436233,  0.16483989, -0.15923648,
                             0.21872145, -0.05275334,  0.10091283, -0.17636017,  0.04387883, -0.13616001,
                             0.19830351, -0.11750502,  0.04748888, -0.16028752,  0.22837654, -0.1166443,
                             0.19825184, -0.2125061 ,  0.07986435, -0.16658635,  0.17372507, -0.07108562,
                             0.20191163, ],
       relative_tolerance=1e-7,
       absolute_tolerance=1e-7,
    )

def test_tuple_ize():
    for a_in, a_out in tuple_ize.tuple_ize.unit_test_input_output_pairs:
         assert tuple_ize.tuple_ize(a_in) == a_out
