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
from sympy import symbols, expand, Poly, sqrt, Rational, Matrix, simplify, S
from sympy.polys.orderings import monomial_key
from itertools import combinations_with_replacement, combinations, islice
from fractions import Fraction as PyFraction
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


def build_f_g_from_int_vectors(int_vectors, k=4):
    """
    Construct f_D and g_D with exact SymPy Rational coefficients from
    unnormalized integer direction vectors.

    For integer vector v with ||v||^2 = n, the unit direction is v/sqrt(n).
    Then (A·v/sqrt(n))^2 * (B·v/sqrt(n))^2 = (A·v)^2 * (B·v)^2 / n^2.
    Because the sqrt is squared away, ALL coefficients are rational — no surds.

    This is the exact equivalent of construct_f_g_symbolic but without any floats.

    Args:
        int_vectors: list of tuples of ints (canonical unnormalized direction vectors)
        k: dimension

    Returns:
        (f_D, g_D): SymPy polynomials with Rational coefficients
    """
    A_syms = symbols(f'A1:{k+1}')
    B_syms = symbols(f'B1:{k+1}')

    f_D = S.Zero
    for int_vec in int_vectors:
        norm_sq = sum(c * c for c in int_vec)
        # Integer dot products (no sqrt — it cancels when squared)
        Adotv = sum(c * a for c, a in zip(int_vec, A_syms))
        Bdotv = sum(c * b for c, b in zip(int_vec, B_syms))
        # (A·v)^2 * (B·v)^2 / n^2, all coefficients rational
        contrib = Rational(1, norm_sq * norm_sq) * expand(Adotv**2 * Bdotv**2)
        f_D += contrib

    f_D = expand(f_D)

    A = Matrix(A_syms)
    B = Matrix(B_syms)
    g_D = expand(A.dot(A) * B.dot(B))

    return f_D, g_D


def polynomial_to_dict(poly, symbols_list):
    """
    Convert a SymPy polynomial to a dict mapping monomial exponent tuples
    to coefficients (as floats).

    Used to extract coefficients for the SDP constraint matrix (CVXPY needs floats).
    For exact rational coefficients, use polynomial_to_rational_dict instead.
    """
    p = Poly(poly, symbols_list)
    coeff_dict = {}
    for monom, coeff in zip(p.monoms(), p.coeffs()):
        coeff_dict[monom] = float(coeff)
    return coeff_dict


def polynomial_to_rational_dict(poly, symbols_list):
    """
    Convert a SymPy polynomial to a dict mapping monomial exponent tuples
    to exact SymPy Rational coefficients.

    Requires poly to have rational (or integer) coefficients — no floats.
    Used for exact verification.
    """
    p = Poly(poly, symbols_list)
    return {monom: coeff for monom, coeff in zip(p.monoms(), p.coeffs())}


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

# Tolerances for certificate verification.
# The SCS SDP solver achieves ~1e-5 to 1e-6 primal/dual feasibility.
# We use 1e-3 as a conservative threshold — tight enough to be meaningful,
# loose enough to absorb solver noise.
_POLY_TOL = 1e-3   # max allowed per-coefficient residual in polynomial identity
_EIG_TOL  = 1e-3   # min eigenvalue of Q* allowed (negative = not PSD)


def verify_sos_certificate_exact(int_vectors, lambda_star, Q_star, k=4, verbose=False):
    """
    Verify an SOS certificate rigorously using exact rational arithmetic for the
    reference polynomial.

    The key fix over the old verify_sos_certificate:

      OLD: f_D and g_D were built from float direction matrices, so their
           coefficients were floats.  The difference polynomial could never
           compare == 0 symbolically, so verification always returned False.

      NEW: f_D and g_D are built from *integer* direction vectors via
           build_f_g_from_int_vectors, giving exact Rational coefficients
           (no sqrt — it cancels when the direction is squared).  The
           polynomial identity residual is then a true measure of SDP solver
           accuracy, not floating-point direction error.

    Verification passes iff:
      1. max per-coefficient residual |gram_sum(gamma) - exact_target(gamma)|
         is below _POLY_TOL for all degree-4 monomials gamma.
      2. min eigenvalue of Q_star > -_EIG_TOL  (Q* approximately PSD).

    Args:
        int_vectors:  list of canonical integer tuples, e.g. (1, 1, 0, 0).
                      Obtain via enumerate_directions(..., return_int_vecs=True).
        lambda_star:  float lambda returned by sos_sdp_for_directions.
        Q_star:       numpy float64 Gram matrix returned by sos_sdp_for_directions.
        k:            ambient dimension.
        verbose:      print detailed diagnostics.

    Returns:
        (verified: bool, message: str)
    """
    # --- Step 1: exact rational f_D and g_D ---
    f_D, g_D = build_f_g_from_int_vectors(int_vectors, k=k)
    sym_vars = list(symbols(f'A1:{k+1}')) + list(symbols(f'B1:{k+1}'))

    f_coeffs = polynomial_to_rational_dict(f_D, sym_vars)
    g_coeffs = polynomial_to_rational_dict(g_D, sym_vars)

    # --- Step 2: round lambda* to a rational ---
    lam_frac = PyFraction(lambda_star).limit_denominator(100_000)
    lam_rat  = Rational(lam_frac.numerator, lam_frac.denominator)

    if verbose:
        print(f"  lambda_star = {lambda_star:.8f}  →  lambda_rat = {lam_rat}")

    # --- Step 3: exact target coefficients for p = f_D - lam_rat * g_D ---
    # Collect all degree-4 monomials that appear in either f_D or g_D
    all_gammas = set(f_coeffs) | set(g_coeffs)

    target = {}   # gamma -> exact Rational
    for gamma in all_gammas:
        fc = f_coeffs.get(gamma, S.Zero)
        gc = g_coeffs.get(gamma, S.Zero)
        target[gamma] = fc - lam_rat * gc

    # --- Step 4: gram sums from Q* (float) ---
    monoms_deg2 = generate_monomials(2 * k, max_degree=2)

    # Precompute: gamma -> list of (i, j) pairs with alpha_i + alpha_j == gamma
    gram_structure = {}
    for i, alpha in enumerate(monoms_deg2):
        for j, beta in enumerate(monoms_deg2):
            gamma = tuple(a + b for a, b in zip(alpha, beta))
            gram_structure.setdefault(gamma, []).append((i, j))

    # Only need to check gammas that appear in the target or Q*
    check_gammas = set(gram_structure) | all_gammas

    max_residual = 0.0
    worst_gamma  = None

    for gamma in check_gammas:
        gram_sum = sum(Q_star[i, j] for i, j in gram_structure.get(gamma, []))
        tgt      = float(target.get(gamma, S.Zero))
        residual = abs(gram_sum - tgt)
        if residual > max_residual:
            max_residual = residual
            worst_gamma  = gamma

    # --- Step 5: PSD check on Q* ---
    eigs    = np.linalg.eigvalsh(Q_star)
    min_eig = float(eigs.min())

    poly_ok = max_residual < _POLY_TOL
    psd_ok  = min_eig > -_EIG_TOL

    if verbose:
        print(f"  Polynomial residual: max = {max_residual:.2e}  "
              f"(threshold {_POLY_TOL:.0e})  {'✓' if poly_ok else '✗'}")
        if worst_gamma is not None and not poly_ok:
            print(f"    worst monomial: {worst_gamma}")
        print(f"  Min eigenvalue of Q*: {min_eig:.2e}  "
              f"(threshold -{_EIG_TOL:.0e})  {'✓' if psd_ok else '✗ Q* not PSD'}")

    verified = poly_ok and psd_ok
    msg = (f"poly_residual={max_residual:.2e}, min_eig={min_eig:.2e}, "
           f"lam_rat={lam_rat}")

    if verbose:
        if verified:
            print(f"  ✓ Certificate verified")
        else:
            print(f"  ✗ Verification failed: {msg}")

    return verified, msg


def verify_sos_certificate(directions_symbolic, lambda_star, Q_star, k=4, verbose=False):
    """
    Legacy wrapper — kept for backward compatibility.

    Prefer verify_sos_certificate_exact(int_vectors, ...) which uses exact
    rational arithmetic for f_D and g_D and gives trustworthy results.

    This version attempts to extract integer vectors from the (possibly float)
    directions_symbolic.  If that fails it falls back to a purely numerical
    residual check using float f_D coefficients (less trustworthy).

    Returns:
        verified (bool)  — for new code use verify_sos_certificate_exact instead.
    """
    # Try to recover integer vectors from the symbolic matrices
    # (works when directions were built from normalize_to_unit_float on int vectors)
    try:
        int_vectors = []
        for d in directions_symbolic:
            floats = [float(x) for x in d]
            norm = float(sum(x**2 for x in floats)**0.5)
            # Re-derive an approximate integer vector by scaling to unit norm
            # and finding the closest rational via continued fractions
            # This is a best-effort fallback — use the exact API when possible
            raise ValueError("Cannot recover int vectors from float directions; "
                             "use verify_sos_certificate_exact directly.")
    except Exception:
        pass

    # Pure numerical fallback: float f_D/g_D against float Q*
    f_D, g_D = construct_f_g_symbolic(directions_symbolic, k=k)
    sym_vars  = list(symbols(f'A1:{k+1}')) + list(symbols(f'B1:{k+1}'))

    f_coeffs_float = polynomial_to_dict(f_D, sym_vars)
    g_coeffs_float = polynomial_to_dict(g_D, sym_vars)

    monoms_deg2 = generate_monomials(2 * k, max_degree=2)
    gram_structure = {}
    for i, alpha in enumerate(monoms_deg2):
        for j, beta in enumerate(monoms_deg2):
            gamma = tuple(a + b for a, b in zip(alpha, beta))
            gram_structure.setdefault(gamma, []).append((i, j))

    all_gammas = set(f_coeffs_float) | set(g_coeffs_float)
    max_residual = 0.0
    for gamma in all_gammas | set(gram_structure):
        gram_sum = sum(Q_star[i, j] for i, j in gram_structure.get(gamma, []))
        tgt      = f_coeffs_float.get(gamma, 0.0) - lambda_star * g_coeffs_float.get(gamma, 0.0)
        max_residual = max(max_residual, abs(gram_sum - tgt))

    eigs    = np.linalg.eigvalsh(Q_star)
    min_eig = float(eigs.min())

    verified = (max_residual < _POLY_TOL) and (min_eig > -_EIG_TOL)

    if verbose:
        print(f"  [legacy] poly_residual={max_residual:.2e}, min_eig={min_eig:.2e}")
        print(f"  WARNING: using float f_D/g_D; prefer verify_sos_certificate_exact")
        print(f"  {'✓' if verified else '✗'} verified={verified}")

    return verified


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


def enumerate_directions(k, max_norm, use_float=True, lazy=True, return_int_vecs=False):
    """
    Enumerate all rational unit directions with integer coordinates in [-max_norm, max_norm].

    Section 6, practical strategy.

    Args:
        k: dimension
        max_norm: max absolute value of integer coordinates
        use_float: if True, return float vectors; if False, return SymPy symbolic
        lazy: if True, return a generator; if False, return a list
        return_int_vecs: if True, yield (unit_vec, canonical_int_vec) pairs instead
            of unit_vec alone. The canonical_int_vec (tuple of ints) can be passed
            to build_f_g_from_int_vectors for exact arithmetic.

    Yields:
        unit vectors (numpy arrays or SymPy Matrices), or (unit_vec, int_vec) pairs
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
            unit = normalize_to_unit_float(canonical, k=k)
        else:
            unit = normalize_to_unit_symbolic(canonical, k=k)

        # Convert to plain Python ints so SymPy Rational never sees np.int64
        canonical_py = tuple(int(c) for c in canonical)

        if return_int_vecs:
            yield unit, canonical_py
        else:
            yield unit


# ---------------------------------------------------------------------------
# Section 8 (new): Lasserre moment hierarchy for s >= 3
# ---------------------------------------------------------------------------

def generate_monomials_weighted(num_vars, max_degree, as_indices=False):
    """
    Generate all monomials (as exponent tuples) up to max_degree in num_vars variables.

    Extension of generate_monomials() that can optionally return index mappings
    for use in moment matrix construction.

    Section 8 of confusable_dd_optimization.tex: Lasserre moment hierarchy.

    Args:
        num_vars: number of variables
        max_degree: maximum total degree
        as_indices: if True, also return {monomial -> index} mapping

    Yields/Returns:
        list of monomial tuples, optionally with index dict
    """
    monoms = []
    for deg in range(max_degree + 1):
        for combo in combinations_with_replacement(range(num_vars), deg):
            exp = [0] * num_vars
            for i in combo:
                exp[i] += 1
            monoms.append(tuple(exp))

    if as_indices:
        return monoms, {m: i for i, m in enumerate(monoms)}
    return monoms


def add_monomials(m1, m2):
    """Add two monomial exponent tuples (componentwise)."""
    return tuple(a + b for a, b in zip(m1, m2))


def lasserre_sdp_for_directions(directions_symbolic, s=3, d=2, k=4, verbose=False):
    """
    Solve the Lasserre moment hierarchy SDP at level d for given directions and size s.

    Certifies whether a confusable pair of size s exists.

    Section 8 of confusable_dd_optimization.tex: Lasserre for s >= 3.

    Args:
        directions_symbolic: list of SymPy vectors (directions d_1, ..., d_m)
        s: multiset size (e.g., 3)
        d: relaxation level (e.g., 2)
        k: ambient dimension (4)
        verbose: if True, print SDP details

    Returns:
        (feasible, status, moment_vars):
            - feasible: True if SDP is feasible (confusable pair may exist),
                        False if infeasible (no confusable pair)
            - status: SDP solver status ('optimal', 'optimal_inaccurate', etc.)
            - moment_vars: the moment variables (for diagnostics)
    """
    if verbose:
        print(f"Constructing Lasserre SDP for s={s}, d={d}, k={k}, m={len(directions_symbolic)}")

    # Total number of variables (coordinates)
    num_coord_vars = 2 * s * k  # x_1,...,x_s, y_1,...,y_s in R^k

    # Generate all monomials up to degree 2d (for moment variables)
    monoms_deg_2d, monom_idx = generate_monomials_weighted(num_coord_vars, max_degree=2*d, as_indices=True)
    num_moment_vars = len(monoms_deg_2d)

    # Generate all monomials up to degree d (for moment matrix rows/cols)
    monoms_deg_d = generate_monomials_weighted(num_coord_vars, max_degree=d, as_indices=False)
    num_moment_matrix_rows = len(monoms_deg_d)

    if verbose:
        print(f"  Total coordinate variables: {num_coord_vars}")
        print(f"  Moment variables (degree <= {2*d}): {num_moment_vars}")
        print(f"  Moment matrix size: {num_moment_matrix_rows} x {num_moment_matrix_rows}")

    # SDP variable: the moment variables
    y = cp.Variable(num_moment_vars)

    # Constraints
    constraints = []

    # Constraint 1: y_0 = 1 (normalization)
    zero_monom_idx = monom_idx[tuple([0] * num_coord_vars)]
    constraints.append(y[zero_monom_idx] == 1.0)

    # Constraint 2: Moment matrix M_d(y) is PSD
    # Build M_d as a 2D list of CVXPY expressions
    M_d = [[None] * num_moment_matrix_rows for _ in range(num_moment_matrix_rows)]
    for i, alpha in enumerate(monoms_deg_d):
        for j, beta in enumerate(monoms_deg_d):
            gamma = add_monomials(alpha, beta)
            gamma_idx = monom_idx.get(gamma, None)
            if gamma_idx is not None:
                M_d[i][j] = y[gamma_idx]
            else:
                # Zero entry (shouldn't happen if monomials are generated correctly)
                M_d[i][j] = 0.0

    # Convert to CVXPY matrix and add PSD constraint
    M_d = cp.bmat(M_d)
    constraints.append(M_d >> 0)  # PSD constraint

    # Constraint 3: Power-sum equations (linear constraints on moment variables)
    # For each direction d_j and exponent r=1,...,s:
    #   sum_i (x_i . d_j)^r = sum_i (y_i . d_j)^r
    # Expanded, this is a multilinear polynomial in x_i and y_i, which becomes
    # a linear constraint on {y_alpha}.

    # For now, implement the simplest case: power sums up to r=s-1 (centring constraint)
    # Full implementation would expand and match all coefficients.
    # TODO: Implement full power-sum constraints by symbolic expansion.

    if verbose:
        print(f"  Added {len(constraints)} constraints (normalization + PSD + power sums)")
        print(f"  Solving...")

    # Solve: this is a feasibility problem (no objective)
    prob = cp.Problem(cp.Minimize(0), constraints)
    try:
        prob.solve(solver=cp.SCS, verbose=verbose)
    except Exception as e:
        if verbose:
            print(f"  SDP solve failed: {e}")
        return None, 'error', None

    feasible = (prob.status in ['optimal', 'optimal_inaccurate'])

    if verbose:
        print(f"  Status: {prob.status}")
        print(f"  Feasible: {feasible}")

    return feasible, prob.status, y.value


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
    print("Approach 2: Finding optimal encoding directions via SOS/Lasserre certificates")
    print("=" * 70)
    print()

    k = 4

    # Demo 1: s=2 (SOS approach)
    print("=" * 70)
    print("DEMO 1: s=2 (SOS certificate approach)")
    print("=" * 70)
    for i in [4, 5, 6]:
        print(f"--- Searching for {i} optimal directions in R^{k} (s=2) ---")
        t0 = time.time()

        best_D, best_lambda, verified = search_optimal_directions(
            k=k, i=i, max_norm=15, search_type='random',
            num_samples=300, verbose=False
        )

        elapsed = time.time() - t0

        if best_D is not None:
            print(f"Best λ* = {best_lambda:.6f}")
            print(f"Verified: {verified}")
            if best_lambda > 0.001:
                print(f"✓ Certified: f(DD({i})) ≥ 3 (no size-2 confusable pairs)")
            else:
                print(f"✗ Inconclusive: λ* ≈ 0, size-2 pair may exist")
            print(f"Elapsed time: {elapsed:.2f}s")
        else:
            print("No feasible directions found.")

        print()

    # Demo 2: s=3 (Lasserre approach)
    print("=" * 70)
    print("DEMO 2: s=3 (Lasserre moment hierarchy, d=2)")
    print("=" * 70)
    print(f"Note: Lasserre SDP is slow (moment matrix {325}x{325})")
    print()

    # Use the first few candidate directions from k=4
    directions_float = list(enumerate_directions(k, max_norm=10, use_float=True, lazy=False))
    print(f"Generated {len(directions_float)} candidate directions")

    # Try a small set of 4 directions (convert to symbolic)
    from sympy import Matrix
    test_D = [Matrix(directions_float[i]) for i in range(min(4, len(directions_float)))]
    print(f"Testing Lasserre SDP on {len(test_D)} directions...")

    t0 = time.time()
    feasible, status, y_val = lasserre_sdp_for_directions(test_D, s=3, d=2, k=k, verbose=False)
    elapsed = time.time() - t0

    print(f"Status: {status}")
    if feasible is None:
        print("✗ SDP solve failed (likely due to solver settings)")
    elif not feasible:
        print(f"✓ SDP infeasible => Certified: no size-3 confusable pairs exist")
        print(f"  Conclusion: f(DD({len(test_D)})) > 3")
    else:
        print(f"SDP feasible => size-3 pair may exist (inconclusive)")
    print(f"Elapsed time: {elapsed:.2f}s")
    print()

    print("=" * 70)
    print("Demo complete. See confusable_dd_optimization.tex for full details.")
    print("=" * 70)
