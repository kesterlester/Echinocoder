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

def test2():

    scaled_bat_matrix = sp.Matrix([[-2, -1], [0, -1], [1, 1], [-2, 1], [-1, 3], [1, -1]])
    E, O, C, EE, OO = confusable_multisets.analyze_B(scaled_bat_matrix, plot_if_2d=False, show_C_if_plotting = False)

    assert EE.total() == OO.total()
    assert len(EE) != len(OO) # This result might be surprising until you realise that EE and OO are objects of type collection.Counter !
     
