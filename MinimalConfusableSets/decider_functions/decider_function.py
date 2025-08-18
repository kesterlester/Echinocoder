# this file is to be called decider_function
import numpy as np
from scipy.stats import qmc, norm
from scipy.linalg import null_space
import sympy as sp 

def generate_gaussian(M,k):

    gaussian_vectors = np.random.randn(M,k)
    norms = np.linalg.norm(gaussian_vectors, axis=1, keepdims= True)
    B = gaussian_vectors/norms
    
    return(B)
  
def generate_sobol(M, k, sample = "std", spread = "lin"):
    sampler = qmc.Sobol(d=k, scramble = True)
    q = (np.log2(M))

    if q.is_integer():
        positive_vectors = sampler.random_base2(int(q))
    else:
        if sample == "std":
            positive_vectors = sampler.random(M)
        elif sample == "rr":
            q_ceil = np.ceil(q).astype(int)
            all_vectors = sampler.random_base2(q_ceil)
            num_removals = 2**q_ceil - M
            indices = np.random.choice(all_vectors.shape[0], size = M, replace = False)
            positive_vectors = all_vectors[indices]
            #print(f"{num_removals} vectors were randomly removed")
        else:
            raise ValueError(f"unknown sample='{sample}' (use 'std' or 'rr')")

    if spread == "lin":
        vectors = 2 * positive_vectors - 1
    elif spread == "gauss":
        clip = 1e-12
        vectors_clipped = np.clip(positive_vectors, clip, 1-clip)
        vectors = norm.ppf(vectors_clipped)
    else:
        raise ValueError(f"unknown spread='{spread}' (use 'lin' or 'gauss')")

    norms = np.linalg.norm(vectors, axis=1, keepdims=True) 
    B = vectors / norms
    return(B)
    
def generate_electrostatic(M, k, iters=500, learning_rate=0.01, power=2.0):
    B = generate_gaussian(M, k)
    eps = 1e-12

    for q in range(1,iters+1):
        tot_force_on_B = np.zeros((M, k))

        for i in range(M):
            diff = B[i] - B                   
            diff[i] = 0.0
            dist = np.linalg.norm(diff, axis=1) + eps
            denom = dist**(power + 1)         
            force = (diff / denom[:, None]).sum(axis=0)
            tot_force_on_B[i] = force

        B_update = (B + tot_force_on_B * learning_rate)
        B = B_update / np.linalg.norm(B_update, axis = 1, keepdims = True)

    return(B)


def construct_A(L, B):
    if isinstance(L, sp.Matrix): 
        L = np.array(L).astype(float)

    R, M = L.shape
    M_B, k = B.shape
    assert M == M_B, "L and B must work with the same no. of vectors"

    L_rep  = np.repeat(L, k, axis=0)            # (R*k, M)
    BT_rep = np.tile(B.T, (R, 1))               # (R*k, M)
    A = L_rep * BT_rep                          # elementwise
    return A

"""
def check_collapse_random(A,num_trials=1000,tol=1e-12,):
    ns_basis = null_space(A, rcond=tol)

    if ns_basis.shape[1] == 0:
        return False
    
    if ns_basis.shape[1] == 1:
        alpha = ns_basis[:, 0]
        return np.all(np.abs(alpha) > tol)
    
    n = ns_basis.shape[1]
    coeffs = np.random.randn(n, num_trials)          # (n, T)
    candidates = ns_basis @ coeffs                   # (M, T)
    good = (np.abs(candidates) > tol).all(axis=0)    # (T,)
    return bool(good.any())
"""
# POTENTIAL FIXES TO MATMUL ERROR (may be overcompensation): 
# nullspace should be dense with vectors with no non-zero entries, so 1000 trials is EXCESSIVE, but we /n
# want a high degree of certainty 

def check_collapse_random(A, num_trials=1000, tol=1e-12):

    ns_basis = null_space(A, rcond=tol)  # shape (M, n)

    if ns_basis.shape[1] == 0:
        return False

    # remove any ns basis vectors with non-finite entries (NaN or inf)
    finite_checker = np.all(np.isfinite(ns_basis), axis=0)   # length n
    if not finite_checker.all():
        num_inf = int((~finite_checker).sum()) # note 2 self ~ flips bool value
        print(f"Warning: {num_inf} nullspace basis vector(s) contained NaN/Inf, so were removed")
        ns_basis = ns_basis[:, finite_checker]  # take all rows but only finite cols 

    # if everything was removed
    if ns_basis.shape[1] == 0:
        return False

    # normalize cols and remove near-zero columns
    col_norms = np.linalg.norm(ns_basis, axis=0)
    good_norm = col_norms >= 1e-18
    if not good_norm.all():
        num_bad = int((~good_norm).sum())
        print(f"Warning: {num_bad} near-zero nullspace basis vector(s) → removed")
        ns_basis = ns_basis[:, good_norm]
        col_norms = col_norms[good_norm]
    
    # if everything was removed
    if ns_basis.shape[1] == 0:
        return False

    ns_basis = ns_basis / col_norms  

    # if dim=1 the calculation is simple 
    if ns_basis.shape[1] == 1:
        alpha = ns_basis[:, 0]
        # absolute tolerance based on size of largest alpha value 
        tol_abs = max(tol, 1e-12 * np.max(np.abs(alpha)))
        return np.all(np.abs(alpha) > tol_abs)
    
    # if dim > 1 use random comb of nullspace basis vectors to seek one with no non-zero entries
    n = ns_basis.shape[1] # n is number of vectors in nullspace basis (n>1)
    coeffs = np.random.randn(n, num_trials)        # (n, trials)
    candidates = ns_basis @ coeffs                 # (M, trials)

    # absolute tolerance based on largest alpha value (experimental)
    scales = np.max(np.abs(candidates), axis=0) + 1e-18 # check 1e-18 is reasonable 
    tol_abs = np.maximum(tol, 1e-12 * scales) 

    good = (np.abs(candidates) > tol_abs).all(axis=0)
    return bool(good.any())


def prepare_B(k, M, method, iters, learning_rate, power, sample, spread):
    """
    method: 
        - electrostatic [RECOMMENDED] (slower, well spaced out bad bats, lower k (0-8), lower M (0-100)):
            
            -> iters (default=500++ faster than I thought): number of charged particle movement steps towards eqm
            -> learning_rate (=0.01): strength of response of charge to its neighbours (low for high mass charge)
            -> power (=2): 2 for coulomb's law
        
        - sobol (fast, less well spaced out bad bats, higher k, higher M):  
             
                If m is a power of 2, it uses Sobol's well balanced random_base2 method. 
             Else: 
                 -> sample = std: uses STANDARD sobol random(m) method, which may reduce balance
                 -> sample = rr: uses random_base2 method to generate 2^q ≥ m vectors,
                              then RANDOMLY REMOVES excess vectors to obtain exactly m.
                     
             We need the unit vectors to be spaced out on the unit sphere, not the positive
             orthant: 
                   -> spread = lin: uses a LINEAR function 2 * vectors - 1. Bias towards 
                                 (hyper)cube corners
                   -> spread = gauss: uses a GAUSSIAN function. Hypothesised less bias. 
                                       
    """
    
    if method == "sobol":
        B = generate_sobol(M, k, sample, spread)
    elif method == "electrostatic":
        B = generate_electrostatic(M, k, iters, learning_rate, power)
    else:
        raise ValueError(f"unknown method='{method}' (use 'sobol' or 'electrostatic')")
    return B


def decide_collapse(L: sp.Matrix, B: np.ndarray, num_trials, tol) -> bool:
    """
    Returns True if NO collapse (i.e., found a fully nonzero null vector),
    False if collapse or nullity=0 (so no need to find nullspace).
    """
    A = construct_A(L, B)
    return check_collapse_random(A, num_trials, tol)



    
    
 


        
        
        
        
        
        
        
        

   
 
    
    
    
