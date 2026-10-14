# Embed list of real numbers treated as multiset.
# E.g. this method can embed things like:
#
#        [[3,2],[4,1],[-2,1]]
#
# when the above represents the size n=3 multiset of k=2-vectors:
#
#        {{ [3,2],[4,1],[-2,1] }}
#
# Although this implementation claims to only operate on arrays of reals, 
# it might be able to embed complex arrays, albeit to complex outputs.
# This could be used by a complexly compressed method (see below).
# 
# If k==0 the number of outputs is also zero. If k>=1 then the number of outputs should be 
#
#          ORDER(n, k) == n + (k-1)*n*(n+1)/2  
#
# where n is the number of vectors in the set, and k is the dimension of each of those vectors. 
# This order has leading term 
#
#                k*n*n/2.
#
# Here is an $(n,k)=(2,3)$ example:
# Embedding  $\{(3,1,4),(2,2,5)\}$ should generate
# 
#       [9, 3, 5,     20, 13, 25, 8, 6]
#
# since the polynomial
#
#       $(y + (3+1x+4x^2))(y+(2+2x+5x^2))$
#
# is equal to 
#
#       $y^2 + y^1 (9 x^2 + 3 x + 5) + y^0 (20 x^4 + 13 x^3 + 25 x^2 + 8 x + 6)$
#
# .
# Note that if $k$ were even we could set up the vectors as being complex but of 
# dimension $k/2$.  This would lead to the COMPRESSED_ORDER being 
#
#       COMPRESSED_ORDER = 2*ORDER(n,k/2)
#                        = 2*( n + (k/2-1)*n*(n+1)/2 )
#                        = 2*n + (k-2)*n*(n+1)/2
#                        = 2*n + (k-1)*n*(n+1)/2 - n*(n+1)/2
#                        = ORDER(n,k) + n - n*(n+1)/2
#                        = ORDER(n,k) - ( n*(n+1)/2 - n )
#                        = ORDER(n,k) - n*( (n+1)/2 - 1 )
#                        = ORDER(n,k) - n*(n-1)/2
#
# which is less than ORDER(n,k) by n*(n-1)/2 ... but still has leading term k*n*n/2.

from math import prod
#from itertools import combinations
import numpy as np
#import tools
from sympy import Poly, abc
from MultisetEmbedder import MultisetEmbedder
from typing import Any

def embed_without_flattening(data):
    
    #print("Data is")
    #print(data)

    n = len(data)
    if n==0:
        return []

    k = len(data[0])
    if k==0:
        return []

    polys = [ Poly(vector, abc.x)+Poly(abc.y) for vector in data ]

    #print("polys",polys)

    product = prod(polys)
    #print ("Product",product)

    # # extract coeffs from product slowly:
    # for c in range(1,n+1):
    #      yPower = n-c
    #      for xPower in range(c*(k-1)+1):
    #          coeff = product.coeff_monomial((abc.x)**xPower * (abc.y)**yPower) # This could be made faster by using coeff_monomial((xPower,yPower)) or coeff_monomial((yPower,xPower)) however I have not figured out how to guarantee which order will work reliably.
    #          print(coeff)

    # extract coeffs from product fast:
    coeffs = [ [
             product.coeff_monomial((abc.x)**xPower * (abc.y)**(n-c) ) # n-c==yPower.  This line could be made faster by using coeff_monomial((xPower,yPower)) or coeff_monomial((yPower,xPower)) however I have not figured out how to guarantee which order will work reliably.
               for xPower in range(c*(k-1)+1) ]
               for c in range(1,n+1) ]

    #print("Coeffs",coeffs)
    
    return coeffs

    #xPolys = [ np.polynomial.Polynomial(vector) for vector in data ] 
    #print("x polynomials are: ")
    #for xPoly in xPolys:
    #    print (xPoly)
    # 
    #n=len(data)

    #return [ sum( [ prod(c)  for c in combinations(data, r+1) ] ) for r in range(len(data)) ]

    # return [ sum( [ prod(c)  for c in combinations(data, r+1) ] ) for r in range(len(data)) ]
    
    # Alternative
    # 
    # ans = np.array([ sum( [ prod(c)  for c in combinations(data, r+1) ] ) for r in range(len(data)) ])
    #
    # All embedders have to output lists of real numbers (at least for now) so:
    #if np.iscomplexobj(ans):
    #  ans=tools.expand_complex_to_real_pairs(ans)
    #
    # return ans

class Embedder(MultisetEmbedder):
    def size_from_n_k_generic(self, n: int, k: int) -> int:
        ans = n + (k-1)*n*(n+1)//2  # We only want ordinary division "/" but we use "//" to avoid promoting to float and result is same.
        #print("SIZE == ",ans)
        return ans

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        return MultisetEmbedder.embed_kOne_polynomial(data), metadata
 
    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
    
        n = len(data)
        assert n>1
    
        k = len(data[0])
        assert k>1
    
        flattened_coeffs = [ a for b in embed_without_flattening(data) for a in b ]
    
        EXPECTED_ORDER = self.size_from_n_k(n,k)
        ACTUAL_ORDER = len(flattened_coeffs)
        assert EXPECTED_ORDER == ACTUAL_ORDER
        if EXPECTED_ORDER != ACTUAL_ORDER:
            print("Expected Order",EXPECTED_ORDER)     
            print("Actual Order",ACTUAL_ORDER)     
            raise Exception("Bug in implementation of Bursar's method!")
   
        metadata = None
        return flattened_coeffs, metadata
