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

    from deciders import alpha_attacking_matrix

    ans = alpha_attacking_matrix(L, B)

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
            ]) == alpha_attacking_matrix(L, B)

    )



def test_M5_k2():

    import deciders
    from sympy import Matrix
    
    L_matrix = Matrix([
    [-1, 0, 0,  1, 1],
    [ 0, 0, 1, -1, 1]])
    bat_matrices = [Matrix([
    [-3, -2],
    [ 1, -1],
    [ 1,  2],
    [-3,  2],
    [-1,  4]]), Matrix([
    [ 2,  1],
    [-2,  2],
    [ 2,  0],
    [ 1, -3],
    [ 3,  4]]), Matrix([
    [ 0, -1],
    [ 1,  3],
    [-2, -1],
    [ 1,  1],
    [ 1, -1]]), Matrix([
    [ 2, -2],
    [-1,  0],
    [ 0, -1],
    [-1,  3],
    [-1, -1]])]
    
    has_True = False
    has_False = False

    for i, bat_matrix in enumerate(bat_matrices):
        decider = deciders.Rational_Decider(M=5, k=2, unscaled_bad_bat_matrix=bat_matrix)
        
        vote_for_collapse, null_space = decider.vote_for_collapse_and_null_space(L_matrix)

        print("vote_for_collapse is: ",vote_for_collapse)
        
        # Expect votes_for_collapse=(True, True, True, False) 
        expected = (True, True, True, False)
        assert vote_for_collapse == expected[i]
