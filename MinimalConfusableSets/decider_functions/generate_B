"""
Generates B (size M) from D(size m) produced by eg sobol_2.
"""
import numpy as np
from scipy.linalg import null_space

def generate_B(D):
    """
    Given m vectors in k dimensions, grouped in blocks of (k-1),
    construct M vectors each orthogonal to a block of (k-1) vectors.

    Parameters
    ----------
    D : np.ndarray, shape (m, k)
        The input set of m vectors in R^k.
    k : int
        The dimension of each vector.

    Returns
    -------
    B : np.ndarray, shape (M, k)
        The set of M vectors, each orthogonal to a (k-1)-sized block of vectors from D.
    """
    m = D.shape[0]
    k = D.shape[1]
    B = []   
    step = k - 1

    for i in range(0, m, step):
        block = D[i:i + step]  # Get a block of k-1 vectors (or fewer at the end)
        block = np.array(block)
        ns = null_space(block)  # Find the orthogonal direction(s)
        if ns.size == 0:
            raise ValueError(f"No orthogonal vector found for block starting at index {i}.")
        w = ns[:, 0]  # Take the first basis vector in the nullspace
        B.append(w / np.linalg.norm(w))  #Do i really need to normalise here (D normalised already)?? Inefficient to norm twice

    return np.array(B)
