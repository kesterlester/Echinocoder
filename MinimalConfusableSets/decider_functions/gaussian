import numpy as np

def generate_gaussian(m,k):
    """
    Generate m vectors in k dimensions whose components are sampled from a 
    gaussian dist.
    
    """
    gaussian_vectors = np.random.randn(m,k)

    norms = np.linalg.norm(gaussian_vectors, axis=1, keepdims= True)

    D = gaussian_vectors/norms
    
    return(D)


    
    
