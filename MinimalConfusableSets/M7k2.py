#!/usr/bin/env python3

import confusable_multisets
from sympy import Matrix, MutableDenseMatrix, ImmutableDenseMatrix, Rational, Integer

def demo(L_matrix_string=None):

    M=7
    k=2

    if L_matrix_string==None:
        L_matrix = Matrix([
            [-1, -1, -1, -1, -1,  0,  0],
            [-1, -1,  0,  0,  0, -1,  0],
            [ 0,  0, -1, -1,  0,  0, -1],
            ])
    else:
        L_matrix = Matrix(eval(L_matrix_string))

    assert ( L_matrix.shape[0] == 0 # no rows
             or L_matrix.shape[1] == M )# K cols

    assert M == 7  # unscaled bad bats currently suppose this!
    assert k == 2  # unscaled bad bats currently suppose this!

    unscaled_bad_bats= Matrix([[-11575,  2898],
                               [  7809,  5440],
                               [ -9614, 10710],
                               [  7015,  7050],
                               [  7451, 11043],
                               [ 22430, -6115],
                               [   472, 17542]])

    corn = confusable_multisets.confusable_sets_or_None(L_matrix, unscaled_bad_bats, M) 

    if corn is not None:
        EE, OO, point_in_null_space, scaled_bad_bats = corn
        confusable_multisets.plot(title_add=f"L_mat = {repr(L_matrix)}", M=M, k=k, EE=EE, OO=OO, C=None)

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        demo(sys.argv[1])
    else:
        demo()

