import confusable_multisets
import sympy as sp

Mat = sp.Matrix
con_siz = confusable_multisets.size_of_confusable_multiset_formed_from
EE_siz = confusable_multisets.size_of_EE_multiset_formed_from


def test():

    assert EE_siz(Mat([(1,0),(0,1)])) == 2
    assert EE_siz(Mat([(1,0),(0,1),(1,1)])) == 3
    assert EE_siz(Mat([(1,0),(0,1),(1,2)])) == 4
    assert EE_siz(Mat([(1,0),(0,1),(0,0)])) == 0

    assert con_siz(Mat([(1,0),(0,1)])) == 2
    assert con_siz(Mat([(1,0),(0,1),(1,1)])) == 3
    assert con_siz(Mat([(1,0),(0,1),(1,2)])) == 4
    try:
        _ = con_siz(Mat([(1,0),(0,1),(0,0)]))
        assert False
    except:
        assert True
