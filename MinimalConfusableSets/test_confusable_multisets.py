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
    E, O, C, EE, OO = confusable_multisets.mitm_compute_E_O_C_EE_OO(scaled_bat_matrix)

    assert EE.total() == OO.total()
    assert len(EE) != len(OO) # This result might be surprising until you realise that EE and OO are objects of type collection.Counter !
     
def test_alpha_attacking_matrix():
    """
    def alpha_attacking_matrix(
                    L_matrix : sp.Matrix,  # Any number of rows, but M columns.
                    bat_matrix : sp.Matrix, # M rows, and k columns, so that each row is a bat vector.
                    ) -> sp.Matrix :
        # Modified from the Tom Ruane version in decider_function.
        # This method generates the matrix A for which the solns of A.(vec of alphas)
        # are the same as the solutions to L.(alpha1 w1, alpha2 w2, ... , alphaM, wM)
        # where w1 is the first bat (i.e. first row of bat matrix) and w2 the second,
        # and so on.

        R, M = L_matrix.shape
        M_B, k = bat_matrix.shape
    """

    M = 3
    rows = 2
    k = 4

    import sympy as sp

    L = sp.randMatrix(r=rows, c=M, min=-100, max=100)
    L_rows, L_cols = L.shape
    assert L_cols == M

    B = sp.randMatrix(r=M, c=k, min=-100, max=100)
    B_rows, B_cols = B.shape
    assert B_cols == k

    ans = confusable_multisets.alpha_attacking_matrix(L, B)

    effective_row = 0
    for i in range(rows):
        for kk in range(k):
            for j in range(M):
                assert ans[effective_row, j] == L[i, j] * B[j, kk]
            effective_row += 1
    assert effective_row > 0 or rows == 0 or k == 0 or M == 0

    L = sp.Matrix([[3, 4, 2],
                   [-7, 8, 5]])
    B = sp.Matrix([[3, 1],
                   [2, 9],
                   [-8, 6]])

    assert (

            sp.Matrix([
                [L[0, 0] * B[0, 0], L[0, 1] * B[1, 0], L[0, 2] * B[2, 0], ],
                [L[0, 0] * B[0, 1], L[0, 1] * B[1, 1], L[0, 2] * B[2, 1], ],
                [L[1, 0] * B[0, 0], L[1, 1] * B[1, 0], L[1, 2] * B[2, 0], ],
                [L[1, 0] * B[0, 1], L[1, 1] * B[1, 1], L[1, 2] * B[2, 1], ],
            ]) == confusable_multisets.alpha_attacking_matrix(L, B)

    )


