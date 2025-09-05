from sympy import Matrix, Rational

def find_lambda(u, v):
    """
    Given two vectors u, v (SymPy Matrix, column or row),
    find the smallest integer lambda such that
    w = u + lambda*v has zeros only where both u and v are zero.
    """
    assert u.shape == v.shape, "Vectors must have same shape"
    
    # Build forbidden set F = {-u_i/v_i : u_i and v_i both nonzero}
    forbidden = set()
    for ui, vi in zip(u, v):
        if ui != 0 and vi != 0:
            forbidden.add(Rational(-ui, vi))
    
    # Generate candidate integers
    candidate_lambda=1
    while True:

        #print(f"Trying CAND_LAM = {candidate_lambda}")
        if candidate_lambda not in forbidden:
            return candidate_lambda

        # We will escape the "while True" in a finite time as "forbidden" has finite length:
        assert 2*abs(candidate_lambda) <= len(forbidden)+1

        # Advance candidate_lambda in following sequence:
        # +1, -1, +2, -2, +3, -3, +4, ...
        if candidate_lambda>0:
            candidate_lambda = -candidate_lambda
        else:
            candidate_lambda = -candidate_lambda + 1

def combine_two(u, v):
    """Return coefficients (1, lambda) and the combined vector."""
    lam = find_lambda(u, v)
    w = Matrix(u) + lam * Matrix(v)
    return (1, lam), w


def combine_many(vectors):
    """
    Given a list of vectors, return a tuple of coefficients
    such that the combined vector has zeros only where all inputs have zeros.
    Also returns the combined vector.
    
    Uses iterative pairwise growth.
    """
    m = len(vectors)
    coeffs = [1] + [0]*(m-1)  # start with coeff for first vector = 1
    w = Matrix(vectors[0])
    
    for j in range(1, m):
        lam = find_lambda(w, vectors[j])
        coeffs[j] = lam
        w = w + lam * Matrix(vectors[j])
    
    return tuple(coeffs), w
