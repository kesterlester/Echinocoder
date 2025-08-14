from sobol import generate_sobol
from generate_B import generate_B
from collapse_checker_random import construct_A, check_collapse_random
import sympy as sp 

# SET PARAMETERS---------------------
m = 40
k = 5

# GENERATE D
D = generate_sobol(m,k)

# GENERATE B
B = generate_B(D)

# DEFINE L (for now). CHECK L CONTAINS k+1 NON ZERO ENTRIES

L = sp.Matrix([
    [1, -1, 0, 1, 0, 1, -1, -1, 1, 0],
    [0, 1, -1, 1, 1, 0, -1, -1, 0, 1]
])

# BUILD MATRIX A 

A = construct_A(L,B)

# CHECK COLLAPSE

result = check_collapse_random(A)

print("Collapse detected?", not result)
