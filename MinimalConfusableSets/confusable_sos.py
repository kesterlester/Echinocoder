"""
confusable_sos.py

SOS/SDP certificate for the absence of size-2 confusable multisets.

For k=2, implements the SDP from Section 4 of confusable_sdp.tex:

    maximise  lambda
    subject to  Q >> 0,
                v^T Q v  =  f(A,B) - lambda * g(A,B)  (polynomial identity)

where
    f(A,B) = sum_j (A.d_j)^2 (B.d_j)^2       [zero iff confusable pair exists]
    g(A,B) = ||A||^2 ||B||^2                   [normalisation]

If lambda* > 0, the certificate proves no size-2 confusable pair exists.
If lambda* = 0, a confusable pair may exist (investigate with confusable_search.py).

Also provides an exact O(2^m) combinatorial check for Proposition 1
(partitioning directions into two perpendicular-complement groups).

Usage:
    python confusable_sos.py              # runs built-in examples
    import confusable_sos as cs
    lam = cs.sos_lower_bound(directions)  # SDP certificate
    found, A, B, X, Y = cs.exact_check_s2(directions)  # exact check
"""

import numpy as np
import cvxpy as cp
from sympy import symbols, expand, Poly
from itertools import combinations_with_replacement, combinations


# ---------------------------------------------------------------------------
# Polynomial utilities
# ---------------------------------------------------------------------------

def monomials_upto(n_vars, max_deg):
    """
    Return all exponent tuples (multi-indices) for monomials of degree
    <= max_deg in n_vars variables, in graded lexicographic order.
    """
    result = []
    for d in range(max_deg + 1):
        for combo in combinations_with_replacement(range(n_vars), d):
            alpha = [0] * n_vars
            for i in combo:
                alpha[i] += 1
            result.append(tuple(alpha))
    return result


def poly_coefficients(expr, sym_vars, max_deg):
    """
    Extract the coefficient of each monomial from a sympy expression.
    Returns {exponent_tuple: float_coefficient} for monomials of degree
    <= max_deg.
    """
    p = Poly(expand(expr), sym_vars)
    return {
        tuple(m): float(c)
        for m, c in zip(p.monoms(), p.coeffs())
        if sum(m) <= max_deg
    }


# ---------------------------------------------------------------------------
# SOS SDP: certify absence of size-2 confusable pairs
# ---------------------------------------------------------------------------

def sos_lower_bound(directions, solver=cp.SCS):
    """
    Find the maximum lambda such that f(A,B) - lambda*g(A,B) is SOS,
    where the polynomials f and g are defined above.

    Works for k=2 (four variables: A1, A2, B1, B2).

    Args:
        directions: array of shape (m, 2), the encoding directions
        solver:     CVXPY solver to use (cp.SCS is free; cp.MOSEK is faster)

    Returns:
        lambda* (float).
        > 0  =>  certified: no confusable size-2 pair exists.
        = 0  =>  inconclusive (confusable pair may or may not exist).
        None =>  SDP solve failed.
    """
    D = np.array(directions, dtype=float)
    m, k = D.shape
    if k != 2:
        raise ValueError("sos_lower_bound is implemented for k=2 only.")

    # Symbolic setup
    A1, A2, B1, B2 = symbols('A1 A2 B1 B2')
    svars = (A1, A2, B1, B2)

    f_sym = sum(
        (D[j, 0] * A1 + D[j, 1] * A2) ** 2 *
        (D[j, 0] * B1 + D[j, 1] * B2) ** 2
        for j in range(m)
    )
    g_sym = (A1 ** 2 + A2 ** 2) * (B1 ** 2 + B2 ** 2)

    f_coeffs = poly_coefficients(f_sym, svars, 4)
    g_coeffs = poly_coefficients(g_sym, svars, 4)

    # Gram matrix for SOS certificate
    # Monomials of degree <= 2 in 4 variables: C(4+2,2) = 15
    monoms = monomials_upto(4, 2)   # degree-<=2 monomials (rows/cols of Q)
    N = len(monoms)                  # = 15 for k=2

    # All degree-<=4 monomials (used for the polynomial identity constraints)
    all_monoms = monomials_upto(4, 4)

    # SDP variables
    Q = cp.Variable((N, N), symmetric=True)  # Gram matrix, must be PSD
    lam = cp.Variable()                       # Lower bound to maximise

    # Build constraints: v^T Q v = f - lam*g  (one constraint per monomial)
    # For monomial gamma: sum_{alpha+beta=gamma} Q[i,j] = f_gamma - lam*g_gamma
    constraints = [Q >> 0]
    for gamma in all_monoms:
        gram_contribution = sum(
            Q[i, j]
            for i, alpha in enumerate(monoms)
            for j, beta in enumerate(monoms)
            if tuple(a + b for a, b in zip(alpha, beta)) == gamma
        )
        f_c = f_coeffs.get(gamma, 0.0)
        g_c = g_coeffs.get(gamma, 0.0)
        # Linear in (Q, lam)
        constraints.append(gram_contribution == f_c - lam * g_c)

    prob = cp.Problem(cp.Maximize(lam), constraints)
    prob.solve(solver=solver, verbose=False)

    if lam.value is None:
        return None
    return float(lam.value)


# ---------------------------------------------------------------------------
# Exact combinatorial check for s=2 (Proposition 1)
# ---------------------------------------------------------------------------

def exact_check_s2(directions, tol=1e-8):
    """
    Exact check for a size-2 confusable pair via Proposition 1:
    try all 2^m partitions of directions into two groups and check
    whether each group has a common perpendicular (null space dim >= 1).

    Works for any k and any m (exponential in m, practical for m <= 25).

    Returns:
        (True,  A, B, X, Y)  if a confusable pair exists, with explicit
                              vectors A, B and example pair X, Y.
        (False, None, None, None, None)  if none exists.
    """
    D = np.array(directions, dtype=float)
    m, k = D.shape

    def null_vector(rows):
        """Return a unit vector in the null space of matrix `rows`,
           or None if the null space is trivial."""
        if len(rows) == 0:
            return np.eye(k)[0]           # whole of R^k, pick e_1
        _, sv, Vt = np.linalg.svd(rows, full_matrices=True)
        # sv has min(len(rows), k) entries; trailing rows of Vt are null space
        n_sv = len(sv)
        threshold = tol * (sv[0] if n_sv > 0 and sv[0] > 0 else 1.0)
        n_null = k - np.sum(sv >= threshold)
        if n_null >= 1:
            return Vt[-1]     # last right singular vector
        return None

    # Enumerate all 2^m partitions (S, complement)
    for r in range(m + 1):
        for S_indices in combinations(range(m), r):
            S = list(S_indices)
            Sbar = [j for j in range(m) if j not in S]

            A = null_vector(D[S] if S else np.zeros((0, k)))
            if A is None:
                continue
            B = null_vector(D[Sbar] if Sbar else np.zeros((0, k)))
            if B is None:
                continue

            # Check A != +-B (non-trivial pair)
            if np.allclose(A, B, atol=tol) or np.allclose(A, -B, atol=tol):
                continue

            # Construct an explicit confusable pair
            # p = (A+B)/2, q = (A-B)/2; X = {p/2, -p/2}, Y = {q/2, -q/2}
            p = (A + B) / 2
            q = (A - B) / 2
            X = np.array([p / 2, -p / 2])
            Y = np.array([q / 2, -q / 2])
            return True, A, B, X, Y

    return False, None, None, None, None


# ---------------------------------------------------------------------------
# Main: worked examples from confusable_sdp.tex
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    np.set_printoptions(precision=4, suppress=True)

    angles = np.linspace(0, np.pi, 4, endpoint=False)
    D4 = np.column_stack([np.cos(angles), np.sin(angles)])
    D3 = D4[:3]
    D2 = D4[:2]

    # ------------------------------------------------------------------
    # SOS SDP certificate (Section 4 table)
    # ------------------------------------------------------------------
    print("=== SOS lower bound lambda* for size-2 injectivity ===")
    for label, D in [('2 dirs (0°, 45°)', D2),
                     ('3 dirs (0°, 45°, 90°)', D3),
                     ('4 dirs (equally spaced)', D4)]:
        lam = sos_lower_bound(D)
        if lam is None:
            verdict = 'SDP failed'
        elif lam > 1e-4:
            verdict = f'CERTIFIED injective on size-2 (lambda*={lam:.4f})'
        else:
            verdict = f'inconclusive (lambda*={lam:.4f}), confusable pair may exist'
        print(f"  {label}: {verdict}")
    print()

    # ------------------------------------------------------------------
    # Exact combinatorial check
    # ------------------------------------------------------------------
    print("=== Exact combinatorial check (Proposition 1) ===")
    for label, D in [('2 dirs', D2), ('3 dirs', D3), ('4 dirs', D4)]:
        found, A, B, X, Y = exact_check_s2(D)
        if found:
            print(f"  {label}: confusable pair EXISTS")
            print(f"    A = {A}")
            print(f"    B = {B}")
            print(f"    X = {X}")
            print(f"    Y = {Y}")
        else:
            print(f"  {label}: no confusable size-2 pair")
    print()

    # ------------------------------------------------------------------
    # Cross-check: the two methods should agree
    # ------------------------------------------------------------------
    print("=== Agreement between SDP and combinatorial check ===")
    for label, D in [('2 dirs', D2), ('3 dirs', D3), ('4 dirs', D4)]:
        lam = sos_lower_bound(D)
        found_exact, *_ = exact_check_s2(D)
        sdp_says_exists = (lam is not None and lam < 1e-4)
        # Both should agree: pair exists iff lambda*=0 iff exact check finds one
        agree = (found_exact == sdp_says_exists)
        print(f"  {label}: exact={found_exact}, sdp_inconclusive={sdp_says_exists}, agree={agree}")
