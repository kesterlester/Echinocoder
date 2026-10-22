from fractions import Fraction
from itertools import pairwise
import numpy as np

from EncDec import array_to_lin_comb
from EncDec import barycentric_subdivide
from EncDec import LinComb, MonoLinComb
from tools import numpy_array_of_frac_to_str

def test_BarycentricSubdivide_split_not_preserving_scale():
    print("###########################################")

    pIn = LinComb([MonoLinComb(-1, np.array([1, 0, 0])),
                   MonoLinComb(-7, np.array([0, 1, 0])),
                   MonoLinComb(10, np.array([0, 0, 1]))])

    pOut, qOut = barycentric_subdivide(pIn, return_offset_separately=True, preserve_scale=False, debug=True, use_assertion_self_test=True)

    print(f"pOut\n{pOut}\nqOut{qOut}")

    pOut_expected = LinComb([MonoLinComb(11, [0, 0, 1]),
                             MonoLinComb(6, [1, 0, 1])])
    qOut_expected = MonoLinComb(-7, [1, 1, 1])

    assert pOut == pOut_expected
    assert qOut == qOut_expected

    in_array = pIn.to_numpy_array()
    out_array = (pOut+qOut).to_numpy_array() 

    assert np.array_equal(in_array, out_array)

def test_BarycentricSubdivide_no_split_not_preserving_scale():
    print("###########################################")

    pIn = LinComb([MonoLinComb(-1, np.array([1, 0, 0])),
                   MonoLinComb(-7, np.array([0, 1, 0])),
                   MonoLinComb(10, np.array([0, 0, 1]))])

    pOut = barycentric_subdivide(pIn, return_offset_separately=False, preserve_scale=False, debug=True, use_assertion_self_test=True)

    print(f"pOut\n{pOut}")

    pOut_expected = LinComb([MonoLinComb(11, [0, 0, 1]),
                             MonoLinComb(6, [1, 0, 1]),
                             MonoLinComb(-7, [1, 1, 1])])

    assert pOut == pOut_expected

def test_BarycentricSubdivide_split_preserve_scale():
    print("###########################################")


    pIn = LinComb([MonoLinComb(-1, np.array([1, 0, 0])),
                   MonoLinComb(-7, np.array([0, 1, 0])),
                   MonoLinComb(10, np.array([0, 0, 1]))])

    pOut, qOut = barycentric_subdivide(pIn, return_offset_separately=True, preserve_scale=True, debug=True, use_assertion_self_test=True)
    print(f"pOut\n{pOut}\nqOut\n{qOut}")
    pOut_expected = LinComb([MonoLinComb(11, np.array([Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)], dtype=object)), 
                             MonoLinComb(12, np.array([Fraction(1, 2), Fraction(0, 1), Fraction(1, 2)], dtype=object))])
    qOut_expected = MonoLinComb(-21, [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)])

    assert pOut == pOut_expected
    assert qOut == qOut_expected

def test_BarycentricSubdivide_no_split_preserve_scale():
    print("###########################################")


    pIn = LinComb([MonoLinComb(-1, np.array([1, 0, 0])),
                   MonoLinComb(-7, np.array([0, 1, 0])),
                   MonoLinComb(10, np.array([0, 0, 1]))])

    pOut  = barycentric_subdivide(pIn, return_offset_separately=False, preserve_scale=True, debug=True, use_assertion_self_test=True)
    print(f"pOut\n{pOut}")
    pOut_expected = LinComb([MonoLinComb(11, np.array([Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)], dtype=object)), 
                             MonoLinComb(12, np.array([Fraction(1, 2), Fraction(0, 1), Fraction(1, 2)], dtype=object)), 
                             MonoLinComb(-21, np.array([Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)], dtype=object))])

    assert pOut == pOut_expected


def test_MonoLinComb():
    a = MonoLinComb(4, [[1,2],[3,4]])
    b = MonoLinComb(4, np.array([[1,2],[3,4]]))
    c = MonoLinComb(4, np.array([[1,2],[3,5]]))
    d = MonoLinComb(4, np.array([[1,2],[3,4],[5,6]]))
    
    e = MonoLinComb(-21, [Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)])
    f = MonoLinComb(-21, np.array([Fraction(1, 3), Fraction(1, 3), Fraction(1, 3)], dtype="O"))

    assert a == b
    assert a != c
    assert a != d
    assert b != c
    assert b != d
    assert c != d

    assert e == f

def test_LinComb():
    a = LinComb(MonoLinComb(4, [[1,2],[3,4]]))
    b = LinComb(MonoLinComb(4, np.array([[1,2],[3,4]])))
    c = LinComb(MonoLinComb(4, np.array([[1,2],[3,5]])))
    d = LinComb(MonoLinComb(4, np.array([[1,2],[3,4],[5,6]])))

    assert a == b
    assert a != c
    assert a != d
    assert b != c
    assert b != d
    assert c != d

    l1 = LinComb([MonoLinComb(11, np.array([0, 0, 1])),
                  MonoLinComb(6, np.array([1, 0, 1]))])
    l2 = LinComb([MonoLinComb(11, np.array([0, 0, 1])),
                  MonoLinComb(6, np.array([1, 0, 1]))])

    assert l1 == l2


def test_array_to_lin_comb():
    print("###########################################")

    arr = np.asarray([[ 4, 2],
                      [-3, 5],
                      [ 8, 9],
                      [ 2 ,7]])

    enc = array_to_lin_comb(arr)

    expected = LinComb([MonoLinComb(4, [[1,0]
  ,  [0,0]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(2, [[0,1]
  ,  [0,0]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(-3, [[0,0]
  ,  [1,0]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(5, [[0,0]
  ,  [0,1]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(8, [[0,0]
  ,  [0,0]
  ,  [1,0]
  ,  [0,0]]), MonoLinComb(9, [[0,0]
  ,  [0,0]
  ,  [0,1]
  ,  [0,0]]), MonoLinComb(2, [[0,0]
  ,  [0,0]
  ,  [0,0]
  ,  [1,0]]), MonoLinComb(7, [[0,0]
  ,  [0,0]
  ,  [0,0]
  ,  [0,1]])])
    not_expected = LinComb([MonoLinComb(4, [[1,0]
  ,  [0,0]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(2, [[0,1]
  ,  [0,0]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(-3, [[0,11110]
  ,  [1,0]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(5, [[0,0]
  ,  [0,1]
  ,  [0,0]
  ,  [0,0]]), MonoLinComb(8, [[0,0]
  ,  [0,0]
  ,  [1,0]
  ,  [0,0]]), MonoLinComb(9, [[0,0]
  ,  [0,0]
  ,  [0,1]
  ,  [0,0]]), MonoLinComb(2, [[0,0]
  ,  [0,0]
  ,  [0,0]
  ,  [1,0]]), MonoLinComb(7, [[0,0]
  ,  [0,0]
  ,  [0,0]
  ,  [0,1]])])
    print(f"=======================\nArray to lin comb made encoded")
    print(arr)
    print("to")
    print(f"{enc}")
    assert enc == expected
    assert enc != not_expected

    reconstructed_array = enc.to_numpy_array()
    assert np.array_equal(reconstructed_array, arr)



