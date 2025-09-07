import sympy as sp

from same_nullspace_up_to_axis_perms import canonical_nullspace_key

def test1():
    M1 = sp.Matrix([[1, 2, 3],
                    [4, 5, 6]])
    
    M2 = sp.Matrix([[2, 1, 3],
                    [5, 4, 6]])
    
    print(M1.rref()[0])  # different RREF
    print(M2.rref()[0])
    
    k1 = canonical_nullspace_key(M1)
    k2 = canonical_nullspace_key(M2)
    
    print(f"key1 is {k1}")
    print(f"key2 is {k2}")

    assert k1==k2

def test2():
    M1 = sp.Matrix([[1, 2, 3],
                    [4, 5, 6]])
    
    M2 = sp.Matrix([[2, 4, 6], # a double of row 1 of M1, so does not change null space
                    [4, 5, 6]])
    
    k1 = canonical_nullspace_key(M1)
    k2 = canonical_nullspace_key(M2)
    
    print(f"key1 is {k1}")
    print(f"key2 is {k2}")

    assert k1==k2

def test3():
    M1 = sp.Matrix([[1, 2, 3],
                    [4, 5, 6]])
    
    M2 = sp.Matrix([[2, 4, 6], # a double of row 1 of M1, so does not change null space
                    [4, 5, 6], 
                    [4, 5, 6]]) # a second copy of row 2, sl does not change null space
    
    k1 = canonical_nullspace_key(M1)
    k2 = canonical_nullspace_key(M2)
    
    print(f"key1 is {k1}")
    print(f"key2 is {k2}")

    assert k1==k2
