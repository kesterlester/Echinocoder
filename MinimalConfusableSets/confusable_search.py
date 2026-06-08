"""
confusable_search.py

Gradient-based numerical search for confusable multisets.

For given encoding directions and multiset size s, searches for two
distinct s-element multisets of k-vectors whose sorted projections onto
every encoding direction are identical (the confusability condition).

This is a heuristic: if a confusable pair is found it is genuine (the
residual is verified), but failure to find one does not prove injectivity.

Usage:
    python confusable_search.py          # runs built-in examples
    import confusable_search as cs
    result = cs.search_confusable_pair(s=4, k=2, directions=my_dirs)
"""

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import linear_sum_assignment


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def power_sums(X, directions):
    """
    Compute the power-sum moments of multiset X for each encoding direction.

    For a multiset X of s vectors in R^k and m encoding directions, the
    confusability condition (equation (1) in confusable_sdp.tex) requires
    the power sums to match:

        sum_{v in X} (v . d_j)^r,  for j = 1..m and r = 1..s.

    Returns a 1-D array of length m*s.
    """
    s = len(X)
    return np.array(
        [[np.sum((X @ d) ** r) for r in range(1, s + 1)]
         for d in directions]
    ).ravel()


def match_distance(X, Y):
    """
    Hungarian-matching distance between multisets X and Y.

    Returns sum of matched Euclidean distances under the optimal
    assignment.  Zero iff X and Y are the same multiset.
    """
    cost = np.array([[np.linalg.norm(x - y) for y in Y] for x in X])
    row, col = linear_sum_assignment(cost)
    return cost[row, col].sum()


def search_confusable_pair(s, k, directions,
                           n_restarts=500,
                           res_tol=1e-10,
                           sep_tol=1e-2,
                           seed=0,
                           verbose=False):
    """
    Search for a confusable pair of size s in R^k.

    Strategy: for each trial, fix X randomly (zero centroid), then
    optimise Y (random initialisation, NOT near X) to minimise the
    confusability residual ||power_sums(Y) - power_sums(X)||^2.
    If the residual falls below res_tol and the Hungarian distance
    between X and Y exceeds sep_tol, a genuine confusable pair has
    been found.

    Args:
        s:          multiset size to search for
        k:          ambient dimension (k=2 is the primary use case)
        directions: list or array of shape (m, k), the encoding directions
        n_restarts: number of random (X, Y_0) initialisations to try
        res_tol:    residual threshold below which the pair is accepted
        sep_tol:    minimum Hungarian distance required (rejects X≈Y)
        seed:       numpy random seed for reproducibility
        verbose:    if True, print progress

    Returns:
        (X, Y, residual) if a confusable pair is found, else None.
    """
    directions = np.array(directions, dtype=float)
    np.random.seed(seed)

    for trial in range(n_restarts):
        # Fixed X with zero centroid
        X = np.random.randn(s, k)
        X -= X.mean(axis=0)
        target = power_sums(X, directions)

        # Random Y_0, also zero centroid — crucially NOT initialised near X
        Y0 = np.random.randn(s, k)
        Y0 -= Y0.mean(axis=0)

        # Capture target by value to avoid closure issues across iterations
        def loss(y_flat, t=target):
            Y = y_flat.reshape(s, k)
            return np.sum((power_sums(Y, directions) - t) ** 2)

        result = minimize(
            loss, Y0.ravel(), method='L-BFGS-B',
            options={'ftol': 1e-22, 'gtol': 1e-15, 'maxiter': 10000}
        )

        if result.fun < res_tol:
            Y = result.x.reshape(s, k)
            dist = match_distance(X, Y)
            if dist > sep_tol:
                if verbose:
                    print(f"  Found at trial {trial}: residual={result.fun:.1e}, "
                          f"match_dist={dist:.3f}")
                return X, Y, result.fun

        if verbose and (trial + 1) % 100 == 0:
            print(f"  trial {trial + 1}/{n_restarts} ...")

    return None


# ---------------------------------------------------------------------------
# Verification helper
# ---------------------------------------------------------------------------

def verify_confusable_pair(X, Y, directions, atol=1e-4):
    """
    Verify a claimed confusable pair by checking sorted projections.

    Returns True iff, for every direction, the sorted dot-product tuples
    of X and Y agree to within atol.
    """
    for d in directions:
        pX = sorted(X @ d)
        pY = sorted(Y @ d)
        if not np.allclose(pX, pY, atol=atol):
            return False
    return True


# ---------------------------------------------------------------------------
# Main: worked examples from confusable_sdp.tex
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    np.set_printoptions(precision=4, suppress=True)

    angles = np.linspace(0, np.pi, 4, endpoint=False)
    D4 = np.column_stack([np.cos(angles), np.sin(angles)])
    D3 = D4[:3]
    D2 = D4[:2]

    print("Encoding directions:")
    print(f"  D2 (0°, 45°)             = {D2}")
    print(f"  D3 (0°, 45°, 90°)        = {D3}")
    print(f"  D4 (0°, 45°, 90°, 135°) = {D4}")
    print()

    # ------------------------------------------------------------------
    # Exhibit the cyclic s=3 confusable pair for D3 (Section 3 of doc)
    # ------------------------------------------------------------------
    print("=== Cyclic s=3 pair for 3 directions (Section 3) ===")
    a, b, c = 5.5, -5.7, 0.2
    X_cyclic = np.array([[a, b], [c, a], [b, c]])
    Y_cyclic = np.array([[a, c], [b, a], [c, b]])
    ok = verify_confusable_pair(X_cyclic, Y_cyclic, D3)
    res = np.sum((power_sums(X_cyclic, D3) - power_sums(Y_cyclic, D3)) ** 2)
    print(f"  (a,b,c) = ({a}, {b}, {c})")
    print(f"  X = {X_cyclic}")
    print(f"  Y = {Y_cyclic}")
    print(f"  Verified confusable: {ok}  (residual={res:.1e})")
    print()

    # ------------------------------------------------------------------
    # Numerical search table from Section 4
    # ------------------------------------------------------------------
    print("=== Numerical search results ===")
    print(f"  {'s':>3}  {'2 dirs':>12}  {'3 dirs':>12}  {'4 dirs':>12}")
    for s in [2, 3, 4, 5]:
        row = [f"  {s:>3}"]
        for D in [D2, D3, D4]:
            r = search_confusable_pair(s, k=2, directions=D, n_restarts=500)
            if r is not None:
                row.append(f"{'FOUND':>12}")
            else:
                row.append(f"{'not found':>12}")
        print("".join(row))
    print()
    print("  (Note: 'not found' is empirical evidence, not a proof of injectivity.)")
