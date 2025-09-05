import confusable_multisets
import sympy as sp

Mat = sp.Matrix
siz = confusable_multisets.size_of_confusable_multiset_formed_from

def test():

    assert siz(Mat([(1,0),(0,1)])) == 2
    assert siz(Mat([(1,0),(0,1),(1,1)])) == 3
    assert siz(Mat([(1,0),(0,1),(1,2)])) == 4
    assert siz(Mat([(1,0),(0,1),(0,0)])) == 0
