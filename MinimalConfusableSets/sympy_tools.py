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
import math

def strip_zero_rows(M: sp.Matrix) -> sp.Matrix:
    """
    This function returns a copy of M with all-zero rows removed.
    This is needed because of the requirements explained in the docstring at the top of this file.
    """
    nonzero_rows = [i for i in range(M.rows) if any(M[i, j] != 0 for j in range(M.cols))]
    return M[nonzero_rows, :]

def max_pivot_positions_for_viable_stripped_RRE_matrix(RRE_matrix_shape, k,):
    # Assumes RRE matrix has been stripped of rows with all zeros.
    # RRE matrix has R rows and M columns.
    # k is the number of dimensions of space.
    #
    # p(r) <= M - (R-r)k - 1     for 0 <= r < R.   See (*) in docstring at top of this file.
   
    R, M = RRE_matrix_shape

    return (M - (R-r)*k - 1 for r in range(R))

def pivot_positions_are_all_viable_for_stripped_RRE_matrix(RRE_matrix_shape, pivot_positions, k,):
    # Assumes RRE matrix has been stripped of rows with all zeros.
    # k is the number of dimensions of space.
    
    max_pivot_positions = max_pivot_positions_for_viable_stripped_RRE_matrix(RRE_matrix_shape, k)

    return all( pivot_pos <= max_pivot_pos for pivot_pos, max_pivot_pos in zip(pivot_positions, max_pivot_positions))

def max_rows_for_viable_stripped_RRE_matrix(M, k):
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


from typing import Optional, Literal, Union
import numpy as np
from sympy import ImmutableMatrix

def int_tuple_without_zeros(size: int,
                            lam: float = 0,
                            ) -> tuple:
    """
    Create a tuple of length "size" containing only non-zero integers whose closeness to zero is controlled by the parameter "lam".
    """

    if size < 0:
        raise ValueError(f"Vector's length must be a non-negative integer, not {size}.")

    if lam<0:
        raise ValueError(f"Length scale lam must be >=0, not {la,}.")

    return tuple( int(np.random.choice( (-1,+1) )*(poisson+1)) for poisson in np.random.poisson(lam=lam, size=size) )

def normal_int_matrix(
    rows: int,
    cols: int,
    sigma: float = 1.0,
    mean: float = 0.0,
    seed: Optional[Union[int, np.random.SeedSequence]] = None,
    rng: Optional[Union[np.random.Generator, np.random.RandomState]] = None,
    rounding: Literal["half_to_even", "half_away_from_zero"] = "half_to_even",
) -> ImmutableMatrix:
    """
    Create an rows x cols sympy.ImmutableMatrix of integers sampled by:
      1) draw from Normal(mean, sigma)
      2) round to nearest integer
    RNG control:
      - Pass `seed` to get a fresh reproducible generator.
      - Pass an existing `rng` (numpy Generator/RandomState) to control/statefully reuse it.
      - Pass neither for non-deterministic draws.

    Rounding:
      - "half_to_even" uses NumPy's rint (banker's rounding).
      - "half_away_from_zero" rounds 0.5 magnitudes up, preserving sign.

    Raises:
      ValueError if sigma <= 0 or rows/cols not positive.
    """
    if rows <= 0 or cols <= 0:
        raise ValueError("rows and cols must be positive integers.")
    if sigma <= 0:
        raise ValueError("sigma must be > 0.")

    # Choose the random number generator
    if rng is None:
        rng = np.random.default_rng(seed)  # seed=None -> non-deterministic; seed=int -> reproducible

    # Sample
    samples = rng.normal(loc=mean, scale=sigma, size=(rows, cols))

    # Round
    if rounding == "half_to_even":
        rounded = np.rint(samples)               # ties to even
    elif rounding == "half_away_from_zero":
        rounded = np.sign(samples) * np.floor(np.abs(samples) + 0.5)
    else:
        raise ValueError("rounding must be 'half_to_even' or 'half_away_from_zero'.")

    rounded = rounded.astype(int)
    return ImmutableMatrix(rounded.tolist())

import itertools
from sympy import Matrix

def rows_are_vectors_in_general_position(B: Matrix) -> bool:
    """
    Check if the rows of B (vectors) are in general linear position.
    That means every subset of up to d vectors is linearly independent.
    This is not the same thing as them being POINTS in general position
    (see implementation of rows_are_points_in_general_position(B: Matrix) -> bool:).

    Args:
        B: sympy Matrix of shape (R, d), entries in rationals or integers

    Returns:
        True if points are vectors in general position, False otherwise.
    """
    print("checking rows are vectors in general position")

    R, d = B.shape
    max_k = min(R, d)
    
    for k in range(1, max_k + 1):
        for indices in itertools.combinations(range(R), k):
            subset = B[list(indices), :]  # k x d
            if subset.rank() < k:
                return False
    return True


def rows_are_points_in_general_position(B: Matrix) -> bool:
    """
    Check if the rows of matrix B (points in d-dimensional space)
    are point in general linear position. This is the *affine* notion:
    no k points lie in a (k-2)-dimensional affine flat, for k=2..d+1.

    Note this is not the same
    thing as vectors in general position. See other implementation, 
    i.e. implementation of rows_are_vectors_in_general_position(B: Matrix) -> bool.

    Args:
        B: sympy Matrix of shape (R, d), entries in rationals or integers

    Returns:
        True if points are in general position, False otherwise.
    """
    R, d = B.shape

    # Check subsets of size k = 2 ... min(d+1, R)
    for k in range(2, min(d+1, R) + 1):
        for indices in itertools.combinations(range(R), k):
            subset = B[list(indices), :]   # k x d
            base = subset[0, :]            # (1, d)
            rest = subset[1:, :]           # (k-1, d)
            # replicate base row to shape (k-1, d)
            base_repeated = Matrix.vstack(*([base] * (rest.rows)))
            diffs = rest - base_repeated   # (k-1, d)
            if diffs.rank() < k - 1:
                return False

    return True

from warnings import warn
def general_position_integer_bat_matrix(
                                        M: int, # Number of bats
                                        k: int, # Dimension of each bat
                                        seed = 0,
                                        starting_sigma: float = 1,
                                        sigma_growth_factor: float = 1.05, # Must be > 1.0
                                        ) -> ImmutableMatrix:
    """
    Gives an Mxk matrix, each row is a bat vector having k dims.
    Has the special requirement that every k of the bats be in general position.
    """

    sigma = starting_sigma 
    trial_matrix = normal_int_matrix(rows=M, cols=k, seed=seed, sigma=sigma)
    print(f"We have a trial bat matrix with sigma={sigma}")
    while  not rows_are_vectors_in_general_position(trial_matrix):

        seed += 1
        trial_matrix = normal_int_matrix(rows=M, cols=k, seed=seed, sigma=sigma)
        if (seed % 1000) == 0:
            warn(f"Warning, {seed} attempts at making a general position integer bat matrix. (sympy_tools.py)")
        sigma *= sigma_growth_factor # Exponential growth makes matrices tend to have small entries, only growing large when really needed.
        print(f"We have a trial bat matrix with sigma={sigma}")


    #warn(f"seed={seed}, sigma={sigma}")
    print("passed check")
    return trial_matrix

import sympy as sp

def lex_sort_sympy_matrix_by_rows(M: sp.Matrix) -> sp.Matrix:
    return lex_sort_sympy_matrix(M=M, by_cols=False)

def lex_sort_sympy_matrix_by_cols(M: sp.Matrix, by_cols : bool = True) -> sp.Matrix:
    return lex_sort_sympy_matrix(M=M, by_cols=True)

def lex_sort_sympy_matrix(M: sp.Matrix, by_cols: bool) -> sp.Matrix:
    """
    Lexicographically sort the columns or rows of a Sympy Matrix.

    Parameters:
    - M: sympy.Matrix
    - by_cols: True -> sort by columns, else False -> sort by rows.

    Returns:
    - sympy.Matrix with reordered columns or rows
    """

    if by_cols:
        # Extract columns as tuples for comparison
        cols = [tuple(M[:, j]) for j in range(M.cols)]
        # Sort column indices by lex order
        order = sorted(range(M.cols), key=lambda j: cols[j])
        return M[:, order]

    else:  # by rows
        # Extract rows as tuples for comparison
        rows = [tuple(M[i, :]) for i in range(M.rows)]
        # Sort row indices by lex order
        order = sorted(range(M.rows), key=lambda i: rows[i])
        return M[order, :]

def demo():
    rows=3
    cols=6
    mat = normal_int_matrix(rows, cols, seed=0)
    print(f"Here is a random normal int {rows}x{cols} matrix:")
    print(repr(mat))

    mat = normal_int_matrix(rows, cols, seed=1, sigma=100)
    print(f"Here is a random normal int {rows}x{cols} matrix:")
    print(repr(mat))

if __name__ == "__main__":
    demo()
