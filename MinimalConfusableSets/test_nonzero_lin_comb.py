from nonzero_lin_comb import combine_many, combine_two
from sympy import Matrix


def check(vecs,coeffs, result):

    print(f"We are being asked to check coeffs {coeffs} are good for vecs {vecs}.")

    assert len(vecs) == len(coeffs)

    if len(vecs)==0:
        return

    assert len(vecs) >= 1

    for c in coeffs:
        assert c != 0

    running_total = vecs[0]*0

    for v,lam in zip(vecs, coeffs):
        running_total += lam*v

    assert result == running_total

    for i, val in enumerate(running_total):
        
        if val != 0:
            continue

        # Only permit val == 0 if all vecs were zero for that coord:
        assert val == 0
        for vec in vecs:
            assert vec[i] == 0

def test_two():

    u = Matrix([2, 3, 0, 0])
    v = Matrix([-2, 0, 5, 0])
    
    coeffs, w = combine_two(u, v)
    check([u,v], coeffs, w)

def test_three():

    a = Matrix([1, 2, 0])
    b = Matrix([0, 3, 4])
    c = Matrix([5, 0, 6])
    
    coeffs, w = combine_many([a, b, c])
    check([a,b,c], coeffs, w)

def test_tricky():

    a = Matrix([1, 1, 2, 2,])
    b = Matrix([-1, 1, 1, -1,])
    c = Matrix([0, 0, 0, 1,])
    
    coeffs, w = combine_many([a, b, c])
    check([a,b,c], coeffs, w)
