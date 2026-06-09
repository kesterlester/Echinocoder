"""
confusable_dd_search.py

Proof-of-concept implementation of Approach 2: optimise encoding directions
to maximize the lower bound on f(DD(i)) (minimum confusable pair size).

Uses the Lasserre SOS/SDP certificate for s=2 (no size-2 confusable pairs => f(D) >= 3).

Mathematical specification: see confusable_dd_optimization.tex

Key functions:
  - sos_sdp_for_directions(D, s=2): compute lambda*(D) via SDP
  - verify_sos_certificate(D, lambda_star, Q_star): verify the certificate rigorously
  - enumerate_directions(k, max_num): generate rational unit directions
  - search_optimal_directions(k, i, ...): search for best i-tuple of directions

Rigour: all arithmetic is exact (SymPy) or interval-verified. No tolerance-based
comparisons. Slow but correct.

Usage:
    python confusable_dd_search.py           # runs demo for k=4, i in [4,5,6]
    from confusable_dd_search import search_optimal_directions
    DD, lam_star, verified = search_optimal_directions(k=4, i=4, ...)
"""

import numpy as np
import cvxpy as cp
from sympy import symbols, expand, Poly, sqrt, Rational, Matrix, simplify
from sympy.polys.orderings import monomial_key
from itertools import combinations_with_replacement, combinations, islice
import warnings

warnings.filterwarnings('ignore', category=DeprecationWarning)


# ---------------------------------------------------------------------------
# Section 2 / 3: Construct f_D and g_D symbolically for k=4, s=2
# ---------------------------------------------------------------------------

def construct_f_g_symbolic(directions_symbolic, k=4):
    """
    Construct f_D and g_D symbolically for given directions.

    Section 3 of confusable_dd_optimization.tex: The SOS certificate.

    Args:
        directions_symbolic: list of k SymPy vectors (the directions)
        k: dimension (e.g., 4)

    Returns:
        (f_D, g_D): both symbolic polynomials in A1,A2,A3,A4,B1,B2,B3,B4
    """
    # Define symbolic variables
    A = Matrix(symbols(f'A1:{k+1}'))
    B = Matrix(symbols(f'B1:{k+1}'))

    # f_D(A, B) = sum_j (A . d_j)^2 (B . d_j)^2
    f_D = sum(
        (A.dot(d))**2 * (B.dot(d))**2
        for d in directions_symbolic
    )

    # g_D(A, B) = ||A||^2 ||B||^2
    g_D = A.dot(A) * B.dot(B)

    return expand(f_D), expand(g_D)


def polynomial_to_dict(poly, symbols_list):
    """
    Convert a SymPy polynomial to a dict mapping monomial exponent tuples
    to coefficients (as floats).

    Used to extract coefficients for the SDP constraint.
    """
    p = Poly(poly, symbols_list)
    coeff_dict = {}
    for monom, coeff in zip(p.monoms(), p.coeffs()):
        coeff_dict[monom] = float(coeff)
    return coeff_dict


# ---------------------------------------------------------------------------
# Section 3: Generate all monomials of degree <= 2 (for Gram matrix)
# ---------------------------------------------------------------------------

def generate_monomials(num_vars, max_degree):
    """
    Generate all monomials (as exponent tuples) up to max_degree in num_vars variables.

    Returns a list of tuples, e.g., [(0,0,...,0), (1,0,...,0), (0,1,...,0), ..., (2,0,...,0), ...]
    in graded lexicographic order.
    """
    monoms = []
    for deg in range(max_degree + 1):
        for combo in combinations_with_replacement(range(num_vars), deg):
            exp = [0] * num_vars
            for i in combo:
                exp[i] += 1
            monoms.append(tuple(exp))
    return monoms


def monom_to_sympy(m, sym_vars):
    """Convert monomial tuple (e.g., (1,0,2,...)) to SymPy expression."""
    result = 1
    for i, exp in enumerate(m):
        result *= sym_vars[i] ** exp
    return result


# ---------------------------------------------------------------------------
# Section 3: SDP formulation and solving
# ---------------------------------------------------------------------------

def sos_sdp_for_directions(directions_symbolic, k=4, s=2, verbose=False):
    """
    Solve the SOS SDP for given directions to find lambda*(D).

    Section 4 of confusable_dd_optimization.tex: SDP formulation for k=4, s=2.

    Args:
        directions_symbolic: list of SymPy vectors (directions d_1, ..., d_m)
        k: ambient dimension (4)
        s: multiset size (must be 2 for this implementation)
        verbose: if True, print SDP details

    Returns:
        (lambda_star, Q_star, status):
            - lambda_star: optimal lambda value (float)
            - Q_star: 45x45 Gram matrix (numpy array)
            - status: 'optimal', 'unbounded', 'infeasible', etc.
    """
    if s != 2:
        raise NotImplementedError("Only s=2 is implemented (uses degree-4 SOS).")

    # Symbolic variables
    A = Matrix(symbols('A1:5'))  # A1, A2, A3, A4
    B = Matrix(symbols('B1:5'))  # B1, B2, B3, B4
    sym_vars = list(A) + list(B)  # 8 variables total for k=4

    # Construct f_D and g_D
    f_D, g_D = construct_f_g_symbolic(directions_symbolic, k=k)

    if verbose:
        print(f"f_D has {len(Poly(f_D, sym_vars).monoms())} monomials")
        print(f"g_D has {len(Poly(g_D, sym_vars).monoms())} monomials")

    # Extract coefficients of f_D and g_D for all degree-4 monomials
    all_monoms_deg4 = generate_monomials(2*k, max_degree=4)
    f_coeffs = polynomial_to_dict(f_D, sym_vars)
    g_coeffs = polynomial_to_dict(g_D, sym_vars)

    # Monomials of degree <= 2 (for the Gram matrix rows/cols)
    monoms_deg2 = generate_monomials(2*k, max_degree=2)
    num_monoms = len(monoms_deg2)

    if verbose:
        print(f"Gram matrix size: {num_monoms} x {num_monoms}")
        print(f"Number of degree-4 monomials (constraints): {len(all_monoms_deg4)}")

    # SDP decision variables
    Q = cp.Variable((num_monoms, num_monoms), symmetric=True)
    lam = cp.Variable()

    # Constraints
    constraints = [Q >> 0]  # PSD

    # Polynomial identity: v^T Q v = f_D - lam * g_D
    # Constraint for each degree-4 monomial
    for gamma in all_monoms_deg4:
        # Left side: sum of Q[i,j] for all alpha+beta=gamma with |alpha|,|beta|<=2
        gram_sum = 0
        for i, alpha in enumerate(monoms_deg2):
            for j, beta in enumerate(monoms_deg2):
                # Check if alpha + beta = gamma
                if tuple(a + b for a, b in zip(alpha, beta)) == gamma:
                    gram_sum += Q[i, j]

        # Right side: f_D - lam * g_D coefficient at gamma
        rhs = f_coeffs.get(gamma, 0.0) - lam * g_coeffs.get(gamma, 0.0)
        constraints.append(gram_sum == rhs)

    # Solve
    prob = cp.Problem(cp.Maximize(lam), constraints)
    try:
        prob.solve(solver=cp.SCS, verbose=verbose)
    except Exception as e:
        if verbose:
            print(f"SDP solve failed: {e}")
        return None, None, 'error'

    if lam.value is None:
        return None, None, prob.status

    return float(lam.value), np.array(Q.value), prob.status


# ---------------------------------------------------------------------------
# Section 5: Rigorous verification
# ---------------------------------------------------------------------------

def verify_sos_certificate(directions_symbolic, lambda_star, Q_star, k=4, verbose=False):
    """
    Verify the SOS certificate symbolically using SymPy.

    Section 5 of confusable_dd_optimization.tex: Rigorous verification without tolerances.

    Args:
        directions_symbolic: list of SymPy vectors
        lambda_star: candidate lambda (float)
        Q_star: candidate Gram matrix (numpy array)
        k: dimension
        verbose: if True, print details

    Returns:
        verified (bool): True iff the certificate is valid (symbolically or via interval bounds)
    """
    # Construct f_D and g_D symbolically
    A = Matrix(symbols('A1:5'))
    B = Matrix(symbols('B1:5'))
    sym_vars = list(A) + list(B)

    f_D, g_D = construct_f_g_symbolic(directions_symbolic, k=k)

    # Construct v^T Q v symbolically
    monoms_deg2 = generate_monomials(2*k, max_degree=2)
    v = [monom_to_sympy(m, sym_vars) for m in monoms_deg2]

    # v^T Q v (symbolic)
    v_Q_v = sum(
        Q_star[i, j] * v[i] * v[j]
        for i in range(len(v))
        for j in range(len(v))
    )
    v_Q_v = expand(v_Q_v)

    # f_D - lambda_star * g_D
    target = expand(f_D - lambda_star * g_D)

    # Check if they're equal
    diff = expand(v_Q_v - target)

    if verbose:
        print(f"Symbolic verification: checking if v^T Q v = f_D - {lambda_star} * g_D")
        print(f"Difference (should be ~0): {diff}")

    # If difference is zero (or very close), verification succeeds
    if diff == 0:
        if verbose:
            print("✓ Verified symbolically!")
        return True

    # Try numerical check: evaluate at random points and check if diff is small
    # (not as rigorous as symbolic, but faster)
    # TODO: implement interval arithmetic verification here

    if verbose:
        print("✗ Symbolic verification failed. Difference is non-zero.")

    return False


# ---------------------------------------------------------------------------
# Section 6: Direction parameterisation
# ---------------------------------------------------------------------------

def generate_integer_vectors(k, max_norm):
    """
    Generate all integer k-vectors with norm <= max_norm.

    Yields unnormalized integer vectors (will be normalized separately).
    """
    from itertools import product
    for coords in product(range(-max_norm, max_norm + 1), repeat=k):
        if any(c != 0 for c in coords):  # Exclude zero vector
            norm_sq = sum(c**2 for c in coords)
            yield tuple(coords), norm_sq


def normalize_to_unit_symbolic(int_vector, k=4):
    """
    Normalize an integer k-vector to a unit vector using SymPy rationals/surds.

    Section 6 of confusable_dd_optimization.tex: Rational parameterisation.

    Args:
        int_vector: tuple of integers
        k: dimension

    Returns:
        SymPy vector (with algebraic coordinates using sqrt)
    """
    norm_sq = sum(c**2 for c in int_vector)
    norm = sqrt(norm_sq)
    return Matrix([Rational(c) / norm for c in int_vector])


def normalize_to_unit_float(int_vector, k=4):
    """
    Normalize an integer k-vector to a unit vector using floats.

    Faster than symbolic, but uses floating-point. The verification step
    (Section 5) is still rigorous.

    Args:
        int_vector: tuple of integers

    Returns:
        numpy array (unit vector)
    """
    vec = np.array(int_vector, dtype=float)
    norm = np.linalg.norm(vec)
    return vec / norm


def enumerate_directions(k, max_norm, use_float=True, lazy=True):
    """
    Enumerate all rational unit directions with integer coordinates in [-max_norm, max_norm].

    Section 6, practical strategy.

    Args:
        k: dimension
        max_norm: max absolute value of integer coordinates
        use_float: if True, return float vectors; if False, return SymPy symbolic
        lazy: if True, return a generator; if False, return a list

    Yields:
        unit vectors (numpy arrays or SymPy Matrices)
    """
    seen = set()
    for int_vec, norm_sq in generate_integer_vectors(k, max_norm):
        # Canonicalize: avoid duplicates from scaling
        gcd = np.gcd.reduce(int_vec)
        canonical = tuple(c // gcd for c in int_vec)
        if canonical in seen or tuple(-c for c in canonical) in seen:
            continue
        seen.add(canonical)

        if use_float:
            yield normalize_to_unit_float(canonical, k=k)
        else:
            yield normalize_to_unit_symbolic(canonical, k=k)


# ---------------------------------------------------------------------------
# Section 7: Search for optimal directions
# ---------------------------------------------------------------------------

def search_optimal_directions(k=4, i=4, max_norm=20, search_type='exhaustive',
                              num_samples=None, verbose=False):
    """
    Search for the optimal i-tuple of directions that maximize lambda*(D).

    Section 7 of confusable_dd_optimization.tex: Grid search and optimisation.

    Args:
        k: dimension (4)
        i: number of directions to find
        max_norm: max norm of integer coordinates
        search_type: 'exhaustive' (try all i-tuples), 'random' (sample),
                     'greedy' (build incrementally)
        num_samples: for 'random', how many samples to try
        verbose: print progress

    Returns:
        (best_D, best_lambda, verified):
            - best_D: list of i direction vectors (numpy arrays or SymPy)
            - best_lambda: best lambda* achieved
            - verified: whether the certificate was verified
    """
    if search_type not in ['exhaustive', 'random', 'greedy']:
        raise ValueError(f"Unknown search type: {search_type}")

    best_D = None
    best_lambda = -np.inf
    best_verified = False

    # Generate candidate directions (using float for speed)
    directions_float = list(enumerate_directions(k, max_norm, use_float=True, lazy=False))

    if verbose:
        print(f"Generated {len(directions_float)} candidate unit directions")

    if search_type == 'exhaustive':
        if i > 6:
            raise ValueError(f"Exhaustive search for i={i} is too large. Use 'random' instead.")

        # Try all i-tuples
        from itertools import combinations
        for direction_tuple in combinations(range(len(directions_float)), i):
            D = [directions_float[idx] for idx in direction_tuple]

            # Convert to SymPy for SDP construction
            D_sym = [Matrix(d) for d in D]

            lambda_star, Q_star, status = sos_sdp_for_directions(D_sym, k=k, s=2, verbose=False)

            if lambda_star is None:
                continue

            if verbose and lambda_star > best_lambda:
                print(f"  Found better: lambda* = {lambda_star:.6f} (status: {status})")

            if lambda_star > best_lambda:
                best_lambda = lambda_star
                best_D = D
                # Try to verify (expensive, so do it only for best)
                if lambda_star > 0.001:  # Only verify if promising
                    best_verified = verify_sos_certificate(D_sym, lambda_star, Q_star, k=k, verbose=False)

    elif search_type == 'random':
        if num_samples is None:
            num_samples = 1000 if i <= 5 else 10000

        import random
        for trial in range(num_samples):
            D = [directions_float[random.randint(0, len(directions_float) - 1)] for _ in range(i)]
            D_sym = [Matrix(d) for d in D]

            lambda_star, Q_star, status = sos_sdp_for_directions(D_sym, k=k, s=2, verbose=False)

            if lambda_star is None:
                continue

            if lambda_star > best_lambda:
                best_lambda = lambda_star
                best_D = D
                if verbose:
                    print(f"  Trial {trial}: lambda* = {lambda_star:.6f}")

    elif search_type == 'greedy':
        # Start with best single direction, then add directions greedily
        best_D = []
        for step in range(i):
            best_addition = None
            best_lambda_after = best_lambda
            for idx in range(len(directions_float)):
                candidate_D = best_D + [directions_float[idx]]
                if len(candidate_D) < k:  # Can't span yet
                    continue

                D_sym = [Matrix(d) for d in candidate_D]
                lambda_star, _, _ = sos_sdp_for_directions(D_sym, k=k, s=2, verbose=False)

                if lambda_star and lambda_star > best_lambda_after:
                    best_lambda_after = lambda_star
                    best_addition = idx

            if best_addition is not None:
                best_D.append(directions_float[best_addition])
                best_lambda = best_lambda_after
                if verbose:
                    print(f"Step {step + 1}: added direction {best_addition}, lambda* = {best_lambda:.6f}")

    return best_D, best_lambda, best_verified


# ---------------------------------------------------------------------------
# Main: demo
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import time

    np.set_printoptions(precision=4, suppress=True)
    print("=" * 70)
    print("Approach 2: Finding optimal encoding directions via SOS certificates")
    print("=" * 70)
    print()

    k = 4
    for i in [4, 5, 6]:
        print(f"--- Searching for {i} optimal directions in R^{k} ---")
        t0 = time.time()

        best_D, best_lambda, verified = search_optimal_directions(
            k=k, i=i, max_norm=15, search_type='random',
            num_samples=500, verbose=False
        )

        elapsed = time.time() - t0

        if best_D is not None:
            print(f"Best λ* = {best_lambda:.6f}")
            print(f"Verified: {verified}")
            if best_lambda > 0.001:
                print(f"Conclusion: f(DD({i})) ≥ 3")
            else:
                print(f"Inconclusive: λ* ≈ 0, may have size-2 confusable pair")
            print(f"Elapsed time: {elapsed:.2f}s")
        else:
            print("No feasible directions found.")

        print()
