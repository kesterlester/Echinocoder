## Optimisation related to signs.

Suppose we are polynomial zipping the following  three multiset MS1, MS2 and MS3 of two-vectors, where I want you to think of a,b and x,y here as just real numbers, not as vector labels:

MS1 = { (a,x), (a,y), (b,x), (b,y), }
MS2= { (a,x), (a,-x), (a,y), (a,-y), (b,x), (b,-x), (b,y), (b,-y) }
MS3= { (a,x), (a,-x), (a,y), (a,-y), (b,x), (b,-x), (b,y), (b,-y),(-a,x), (-a,-x), (-a,y), (-a,-y), (-b,x), (-b,-x), (-b,y), (-b,-y)  }

MS1  has every top row element  paired with ONE sign of every bottom row element. 
MS2 has every top row element  paired with BOTH signs of every bottom row element. 
MS3  has every top row element  appearing with BOTH signs, each paired with BOTH signs of every bottom row element. 

If MS1, MS2 and MS3 were naively  zipped by existing echinocoder poly zipping embedders, those embedders would not be able to explot the known-existing symmetries that MS1, MS2 and MS3 have irrespective of the numerical values later plaved into a and b or x and  y, and as a consequence the embedding output size for MS3 would be four times the size of that of MS1 etc.

However, careful inspection of the polynomical being formed makes you realise that certain coeffs in the polynomial are FORCED to zero (or values that depend on other polynomials) by the above symmetries, and so actually fewer coefficients need to be recorded to recover the list.

E.g. the polynomial for MS2 = { (top_n, bot_n) | n in [0,7] } would be  

poly(z, MS2) = product_{n=0}^{7} ( z-(top_n + i bot_n) ) 
    = product_{n=0 to 7 modulo 2}   (z-(top_n + i bot_n))  (z-(top_n - i bot_n)) 
 = product_{n=0 to 7 modulo 2}  (z^2 - 2 top_n + |top_n + i bot_n|^2) 

and look carfully at that polynomial factor at the end:

 (z^2 - 2 top_n + |top_n + i bot_n|^2) 

... it has only two REAL coefficients, not two COMPLEX coefficients which it would have if MS2 had had  no symmetry at all.  This is a factor 2 saving in coefficeints.  

Similarly MS3 has a factor 4 saving in coefficients.

In fact, because of these, you can actually fit ALL of MS1 , MS2 and MS3 into the same output embedding space size.
