#CAUTION: THIS IS AN OLDER, SLOWER, LESS ROBUST VERSION THAN THE FUNC IN decider_function.py

"""
Sobol Sampling
------------------------
Generates m k-dimensional unit vectors in the positive orthant using the 
Sobol low-discrepancy sequence.

If m is a power of 2, it uses Sobol's well balanced random_base2 method. 

Else: 
    sample = std: uses STANDARD sobol random(m) method, which may reduce balance
    sample = rr: uses random_base2 method to generate 2^q ≥ m vectors,
                 then RANDOMLY REMOVES excess vectors to obtain exactly m.
        
We need the unit vectors to be spaced out on the unit sphere, not the positive
orthant: 
      spread = lin: uses a LINEAR function 2 * vectors - 1. Bias towards 
                    (hyper)cube corners
      spread = gauss: uses a GAUSSIAN function. Hypothesised less bias. 

"""
from scipy.stats import qmc, norm
import numpy as np

def generate_sobol(m, k, sample = "rr", spread = "gauss"):
    sampler = qmc.Sobol(d=k, scramble = True)
    
    q = (np.log2(m))
    
    if q.is_integer():
        
        positive_vectors = sampler.random_base2(int(q))
        
    else: 
    
        if sample == "std":
            
            print("WARNING: Sobol balance only guaranteed for m a power of 2")
            positive_vectors = sampler.random(m)

        if sample == "rr":
         
            print("WARNING: Sobol vectors will be generated for the nearest power of \
        2, then vectors will be randomly removed to reduce to m")
            
            q_ceil = np.ceil(q).astype(int)
        
            all_vectors = sampler.random_base2(q_ceil)
            
            num_removals = 2**q_ceil - m
            
            indices = np.random.choice(all_vectors.shape[0], size = m, replace = False)
            positive_vectors = all_vectors[indices]
            
            print(f"{num_removals} vectors were randomly removed")
    
    if spread == "lin":
        vectors = 2 * positive_vectors - 1
    
    if spread == "gauss": 
        clip = 1e-12
        vectors_clipped = np.clip(positive_vectors, clip, 1-clip)
        vectors = norm.ppf(vectors_clipped)
    
    norms = np.linalg.norm(vectors, axis=1, keepdims=True) 
    D = vectors / norms
    
    return(D)
    
