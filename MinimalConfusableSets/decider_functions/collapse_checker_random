import numpy as np
import sympy as sp
from scipy.linalg import null_space


def construct_A(L, B):
    """
    Construct the matrix A from matrix L (R x M) and vector 
    matrix B (M x k) where B[i] is the ith vector in R^k
    """
    if isinstance(L, sp.Matrix): 
        L = np.array(L).astype(float)  # convert sympy to float array

    R, M = L.shape
    M_B, k = B.shape
    assert M == M_B, "L and B must work with the same no. of vectors"

    A = np.zeros((R * k, M))
    for r in range(R):
        for s in range(k): 
            for t in range(M):
                A[r * k + s, t] = L[r, t] * B[t, s]
    return A


def check_collapse_random(A, tol=1e-12, num_trials=1000):
    """
    Checks whether the nullspace of A contains any vector 
    (not just basis vectors) with all nonzero entries using random combinations.

    Parameters
    ----------
    A : numpy.ndarray
        The matrix A constructed from L and B
    tol : float
        Tolerance to treat small entries as zero
    num_trials : int
        Number of random linear combinations to test

    Returns
    -------
    True if any linear combination of nullspace basis vectors has all nonzero entries,
    False otherwise
    """
    ns_basis = null_space(A, rcond=tol)

    if ns_basis.shape[1] == 0:
        print("nullity = 0")
        return False
    
    if ns_basis.shape[1] == 1:
        alpha = ns_basis[:, 0]
        return np.all(np.abs(alpha) > tol)
    
    # For each trial, generate a random combination of nullspace basis vectors    
    n = ns_basis.shape[1]
    for _ in range(num_trials):
        coeffs = np.random.randn(n)  # step 1: random weights
        candidate = ns_basis @ coeffs  # step 2: form linear combination
        if np.all(np.abs(candidate) > tol):  # step 3: check for all non-zero entries
            return True
    
    return False
