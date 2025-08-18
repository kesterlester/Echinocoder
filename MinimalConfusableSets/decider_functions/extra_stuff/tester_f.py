from generate_B import generate_B
from collapse_checker_random import construct_A, check_collapse_random
import sympy as sp 
from importlib import import_module
    
L = sp.Matrix([
    [1, -1, 0, 1, 0, 1, -1, -1, 1, 0],
    ])

def f(m, k, method, L):
    """
    Approaching a decider function format. Checks whether the matches stored
    in L can be executed without collapsing O <-> E. 
    NOTE: ideally we remove m by assuming m = #cols of L*(k-1) to match decider
    function specification ???

    Parameters
    ----------
    m : integer
        NO. OF PROJECTION VECTORS IN SET D.
    k : integer
        DIMENSION.
    method : string
        METHOD TO DEFINE m EQUALLY SPACED VECTORS TO GENERATE D.
        OPTIONS (so far) ARE: "sobol", "gaussian", "electrostatic"
    L : sympy.Matrix
        MATRIX OF MATCHES (R x M)

    Returns
    -------
    Collapse detected? True OR False
    """
    method_module = import_module(method)
    generate_method = getattr(method_module, f"generate_{method}")

    D = generate_method(m, k)
    B = generate_B(D)
    A = construct_A(L, B)
    
    result = check_collapse_random(A)
    
    print("Collapse detected?", not result)

# example
f(40, 5, "sobol", L)
    
    
 
