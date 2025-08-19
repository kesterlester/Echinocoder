"""
A viable RRE-form match-matrix "RRE" is one that can be solved for

    RRE . ( alpha1 w1, alpha2 w2, ... , alphaM wM)^T = 0

for fixed general position w1, w2, ... , wM (each in k-dims) with
all alphas being non-zero.  For short "without causing collapse".

Any viable RRE-form match-matrix with R rows and M columns, when still including pivot
columns, but with zero-rows removed, must have a pivot for row r (in zero-based
position p(r)) satisfying:

    p(r) <= M - (R-r)k - 1     for 0 <= r < R.     (*)

Further more, a viable RRE matrix so trimmed cannot have more than (M-1)/k rows, i.e.:

    R <= (M-1)/k.                                  (**)

Proofs:

The "best case scenario" for solving a vertex-match matric in zero-trimmed RRE-form
is that every element (other than those in the pivot columns, or below the staircase)
are non-zero.

Taking as an example:

    RRE_MATRIX(R=3 rows, M=10 columns), with k=2 dimesions.

In this example, the TIGHTEST (i.e. only just viable) RRE_MATRIX that could be solved
would be as follows (in which "x" means something non-zero, and "a" means anything):

0001x0x0aa     Row r=0, Pivot pos here = 3 = M-3k-1 = M-(3-0)k-1 = M - (R-r)k - 1.
000001x0xa     Row r=1, Pivot pos here = 5 = M-2k-1 = M-(3-1)k-1 = M - (R-r)k - 1. 
00000001xx     Row r=2, Pivot pos here = 7 = M-1k-1 = M-(3-2)k-1 = M - (R-r)k - 1. 
==========
0123456789

Hence in general, p(r) <= M - (R-r)k - 1 as stated in (*).

Note that such pivot positions are necessary but not sufficient for viability.
E.g. the following matrix satisfies the condition but is not solvable as row 0
causes collapse if k=2:

0001000000     Row r=0, Pivot pos here = 3 = M-3k-1 = M-(3-0)k-1 = M - (R-r)k - 1.
000001x0xa     Row r=1, Pivot pos here = 5 = M-2k-1 = M-(3-1)k-1 = M - (R-r)k - 1. 
00000001xx     Row r=2, Pivot pos here = 7 = M-1k-1 = M-(3-2)k-1 = M - (R-r)k - 1.
==========
0123456789

Here is an example for M=10, k=3

1xx0xa0aaa     Row r=0, Pivot pos here = 0 = M-3k-1 = M-(3-0)k-1 = M - (R-r)k - 1.
0001xx0xaa     Row r=1, Pivot pos here = 3 = M-2k-1 = M-(3-1)k-1 = M - (R-r)k - 1. 
0000001xxx     Row r=2, Pivot pos here = 6 = M-1k-1 = M-(3-2)k-1 = M - (R-r)k - 1. 
==========
0123456789

A consequence of this is that an RRE matrix which is not wide enough to accommodate a 
pivot at the required position in the first row MUST cause collapse.

I.e. all matrices would cause collapse if

    p(0) < 0

i.e. if

    M - (R-0)k - 1 < 0

i.e. if

    M - 1 < Rk

i.e. if

    R > (M-1)/k.

Hence a necessary condition for "non-collapse" is

    R <= (M-1)/k.                   (**)

"""

import sympy as sp

def strip_zero_rows(M: sp.Matrix) -> sp.Matrix:
    """
    This function returns a copy of M with all-zero rows removed.
    This is needed because of the requirements explained in the docstring at the top of this file.
    """
    nonzero_rows = [i for i in range(M.rows) if any(M[i, j] != 0 for j in range(M.cols))]
    return M[nonzero_rows, :]

def max_RRE_matrix_pivot_positions(M, k, number_of_pivots):
    p(r) <= M - (R-r)k - 1     for 0 <= r < R.     (*)
    return 

def max_height_for_viable_RRE_matrix_without_zero_rows(M, k):
    """
    This is the maximum number of rows that an RRE matrix (with zero rows stripped!)
    can have if it is to be viable.
    """
    return math.floor((M-1)/k) # See (**) in docstring at top of file.


def some_row_causes_collapse(mat: sp.Matrix, k: int):
        """
        Return True if there exists a row in M that has:
        - at least one non-zero entry, and
        - fewer than k+1 non-zero entries.
        """
        k_plus_1 = k+1
        rows, cols = mat.shape

        for r in range(rows):
            nonzero_count = 0
            for c in range(cols):
                if mat[r, c] != 0:
                    nonzero_count += 1
                    if nonzero_count == k_plus_1:
                        # This row is a good row as it has >= k+1 non-zero entries.
                        break # i.e. stop scanning the columns of this row!
            # We reached the end of the row, so assess what to do:
            if 1 <= nonzero_count < k_plus_1:
                return True #  This is a bad row! It has a number of non-zero entries in {1, 2, ... , k}
            assert nonzero_count==0 or nonzero_count==k_plus_1
            # This row was good, so try the next row.
        # We finished trying rows, so if we got here all rows are good!
        return False
