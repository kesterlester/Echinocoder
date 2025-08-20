from decider_functions.decider_function import generate_electrostatic, generate_sobol
import itertools
import numpy as np


def test_electrostatic():
    for M,k,tol in (
        (6,2,1e-4),
        (50,2,1e-4),
        (50,5,1e-4),
        ):
        B = generate_electrostatic(M=M, k=k, iters=500, learning_rate=0.01, power=2.0)
        for combo in itertools.combinations(range(M), k):
            matrix = B[list(combo), :]   # shape (k, k)
            det = np.linalg.det(matrix)
            assert abs(det) > tol, "Bad bats not lin ind"
                
    
