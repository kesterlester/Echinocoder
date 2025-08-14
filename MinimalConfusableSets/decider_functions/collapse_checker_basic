"""
CHECK COLLAPSE FUNCTION NEEDS UPDATING! It does not check whether a lin comb of the ns basis vectors 
contains a vector with no zero components!

ALSO MUST GENERATE B
"""

import numpy as np
import sympy as sp
from scipy.linalg import null_space

def construct_A(L, B):
    """
    Construct the matrix A from matrix L (R x M) and vector 
    matrix B (M x k) where B[i] is the ith vector in R^k

    Parameters
    ----------
    L : sympy.Matrix
        MATRIX OF MATCHES
    B : numpy.ndarray
        MATRIX OF O AND E GENERATORS 

    Returns
    -------
    A: numpy.ndarray (R*k x M)

    """
    
    if isinstance(L, sp.Matrix): 
        L = np.array(L).astype(float) #convery sympy to float array
        
    R, M = L.shape
    M_B, k = B.shape
    assert M == M_B, "L and B must work with the same no. of vectors"
    
    A = np.zeros((R * k, M))
    
    # CHECK THIS WORKS
    for r in range(R):
        for s in range(k): 
            for t in range(M):
                A[r * k + s, t] = L[r, t] * B[t,s]
                
    return(A)


def check_collapse(A, tol=1e-12):
    """
    Checks if L collapses O and E by computing the nullspace of A and 
    testing whether it contains a vector with no zero entries. 

    Parameters
    ----------
    A : numpy.ndarray (R*k x M)
        MATRIX CONTAINING L AND B INFORMATION
    tol : SET TO DEFAULT??????

    Returns
    -------
    True if no collapse, False if collapse or nullity = 0

    """    
    ns_basis = null_space(A, rcond=tol)
    
    if ns_basis.shape[1] == 0: 
        print("nullity = 0")
        return False
    
    for i in range(ns_basis.shape[1]):
        alpha = ns_basis[:,i]
        if np.all(np.abs(alpha) > tol):
            return True 
        
    return False
