import sympy as sp
import sympy_tools as spt

def canonical_nullspace_key(M: sp.Matrix):
    """
    Compute a canonical key representing the nullspace of M
    up to coordinate permutations.

    M: sympy.Matrix (rectangular allowed), exact rationals.

    Returns: ImmutableMatrix (hashable canonical matrix).
    """
    M = spt.lex_sort_sympy_matrix_by_cols(M)
    rref = M.rref()[0]
    stripped = spt.strip_zero_rows(rref)
    return sp.ImmutableMatrix(stripped)

