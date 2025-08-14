import numpy as np
from gaussian import generate_gaussian

def generate_electrostatic(m, k, iters=500, learning_rate=0.01, power=2.0):
    """
    Generate m k dimensional unit vectors, which are evenly spread out
    by placing a unit + charge on the tip of m vectors on S^(k-1). 
    Initially place vectors randomly on S^(k-1) using a gaussian distribution
    power = 2 for coulomb's law'
    """
    D = generate_gaussian(m, k)
    
    tot_force_on_D = np.zeros((m, k))
   
    for q in range(1,iters+1): 
         
       # Loop through each vector to calculate the total force on it
        for i, x in enumerate(D):
        # This will hold the sum of forces from other vectors on 'x'
            total_force_on_x = np.zeros(k) 
        
        # Loop through all other vectors to calculate the force they exert on 'x'
            for j, y in enumerate(D):
            # We need to make sure we don't calculate the force of a vector on itself
                if i == j:
                    continue
            
                force_on_x = (x - y) / np.linalg.norm(x - y)**(power + 1)
            
                total_force_on_x = total_force_on_x + force_on_x
            
            tot_force_on_D[i] = total_force_on_x
        D_update = (D + tot_force_on_D * learning_rate) 
        D_norms = np.linalg.norm(D_update, axis = 1, keepdims = True)
        D = D_update / D_norms

    return(D)  



        
        
        
        
        
        
        
        

   
 
    
    
    
