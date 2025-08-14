from generate_B import generate_B
from collapse_checker_random import construct_A, check_collapse_random
import sympy as sp 
from importlib import import_module
    
L = sp.Matrix([
    [1, -1, 0, 1, 0, 1, -1, -1, 1, 0],
    [0, 1, -1, 1, 1, 0, -1, -1, 0, 1]
])

def f_multi(m, k, method_1, method_2, method_3, L):
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
    method_1 : string
        1st METHOD TO DEFINE m EQUALLY SPACED VECTORS TO GENERATE D.
        OPTIONS (so far) ARE: "sobol", "gaussian", "electrostatic"
    method_2 : string
        2nd METHOD TO DEFINE m EQUALLY SPACED VECTORS TO GENERATE D.
        OPTIONS (so far) ARE: "sobol", "gaussian", "electrostatic"
    etc. 
    
    L : sympy.Matrix
        MATRIX OF MATCHES (R x M)

    Returns
    -------
    Collapse detected? True OR False
    
    """
    methods = [method_1, method_2, method_3]
    
    for method_name in methods: 
    
        method_module = import_module(method_name)
        generate_method = getattr(method_module, f"generate_{method_name}")
    
        D = generate_method(m, k)
        B = generate_B(D)
        A = construct_A(L, B)
    
        result = check_collapse_random(A)
    
        print(f"Collapse detected by {method_name}?", not result)

# example
f_multi(40, 5, "sobol", "electrostatic", "gaussian",  L)
    
    
 
