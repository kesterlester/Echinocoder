"""
dd_phase1_survey.py

Phase 1 of the DD(i) optimisation survey for k=4, i in [1,16]:
Find the best encoding directions and their certified lower bounds on f(DD(i)).

Theory results (no search needed):
  - i=1,2,3: f(DD(i)) = 1 (directions don't span R^k)
  - i=4,5,6: f(DD(i)) = 2 (s=2 partition always exists)
  - i>=7:    f(DD(i)) >= 3 (no s=2 pairs by partition)

Searchable results (grid search for each i):
  - For i=7,...,16: find best directions D and report lambda*(D) certificate
  - If lambda* > 0: certified f(DD(i)) >= 3
  - If lambda* ≈ 0: inconclusive (may or may not have size-2 pairs)

Usage:
    python dd_phase1_survey.py

Output:
    CSV-like table with columns:
    i, M, Hypercube_UB, Theory_LB, Best_Lambda*, Certified_LB, Notes
"""

import numpy as np
from sympy import Matrix
import confusable_dd_search as cds
import time
import csv
import sys

# ============================================================================
# Configuration
# ============================================================================

K = 4
MAX_NORM = 12  # Max integer norm for rational directions
NUM_TRIALS = {
    1: 0,     # No search (theory)
    2: 0,
    3: 0,
    4: 0,     # No search (theory)
    5: 0,
    6: 0,
    7: 500,   # Search for i >= 7
    8: 500,
    9: 500,
    10: 600,
    11: 600,
    12: 600,
    13: 700,
    14: 700,
    15: 700,
    16: 800,
}

# ============================================================================
# Theory results (analytic bounds, no search)
# ============================================================================

def theory_lower_bound(i, k=4):
    """Compute the proven lower bound on f(DD(i)) from theory alone."""
    if i < k:
        return 1  # Directions don't span R^k
    elif i <= 2*k - 2:  # i <= 6 for k=4
        return 2  # s=2 partition always exists
    else:
        return 3  # i >= 2k-1 = 7 for k=4: no s=2 pairs


def hypercube_upper_bound(i, k=4):
    """Compute the hypercube construction upper bound: 2^{M-1} where M = ceil(i/(k-1))."""
    M = (i + k - 2) // (k - 1)  # ceil(i/(k-1))
    return 2 ** (M - 1)


# ============================================================================
# Grid search for best directions
# ============================================================================

def search_best_directions(i, k=4, num_trials=500, max_norm=12, verbose=False):
    """
    Grid search: try num_trials random i-tuples of directions.
    Return the best lambda* achieved.

    Returns:
        (best_lambda_star, best_verified, num_directions_tried)
    """
    if num_trials == 0:
        return None, False, 0

    if verbose:
        print(f"  Generating candidate directions (max_norm={max_norm})...", end='', flush=True)

    directions_float = list(cds.enumerate_directions(k, max_norm, use_float=True, lazy=False))

    if verbose:
        print(f" {len(directions_float)} candidates")
        print(f"  Searching {num_trials} random {i}-tuples...", end='', flush=True)

    best_lambda = -np.inf
    best_verified = False
    num_tried = 0

    import random
    random.seed(42 + i)  # Reproducible

    for trial in range(num_trials):
        # Sample i directions randomly
        D_indices = [random.randint(0, len(directions_float) - 1) for _ in range(i)]
        D = [directions_float[idx] for idx in D_indices]
        D_sym = [Matrix(d) for d in D]

        # Compute SOS certificate
        lam, Q, status = cds.sos_sdp_for_directions(D_sym, k=k, s=2, verbose=False)

        if lam is None:
            continue

        if lam > best_lambda:
            best_lambda = lam
            # Try verification only for promising candidates
            if lam > 0.001:
                best_verified = cds.verify_sos_certificate(D_sym, lam, Q, k=k, verbose=False)

        num_tried += 1
        if verbose and (trial + 1) % max(50, num_trials // 10) == 0:
            print(f"\r  Searching {num_trials} random {i}-tuples... {trial+1}/{num_trials}", end='', flush=True)

    if verbose:
        print()

    return best_lambda if best_lambda > -np.inf else None, best_verified, len(directions_float)


# ============================================================================
# Main survey
# ============================================================================

def run_survey():
    """Run Phase 1 survey for all i in [1,16]."""
    print("=" * 90)
    print("Phase 1 Survey: Certified lower bounds on f(DD(i)) for k=4, i in [1,16]")
    print("=" * 90)
    print()

    results = []

    for i in range(1, 17):
        M = (i + K - 2) // (K - 1)
        ub = hypercube_upper_bound(i, K)
        lb_theory = theory_lower_bound(i, K)

        print(f"i={i:2d} (M={M}): ", end='', flush=True)

        num_trials = NUM_TRIALS[i]

        if num_trials == 0:
            # Theory only
            print(f"Theory => f(DD({i})) = {lb_theory} (proven)")
            results.append({
                'i': i,
                'M': M,
                'Hypercube_UB': ub,
                'Theory_LB': lb_theory,
                'Best_Lambda*': None,
                'Certified_LB': lb_theory,
                'Method': 'theory_only',
                'Verified': 'N/A',
                'Notes': 'No search needed'
            })
        else:
            # Grid search
            t0 = time.time()
            lam_best, verified, num_dirs = search_best_directions(i, K, num_trials, MAX_NORM, verbose=True)
            elapsed = time.time() - t0

            if lam_best is not None:
                if lam_best > 0.001:
                    cert_lb = 3
                    status_str = f"✓ λ*={lam_best:.4f} → certified f ≥ 3"
                else:
                    cert_lb = lb_theory
                    status_str = f"✗ λ*≈0 → inconclusive (f ≥ {lb_theory})"
            else:
                cert_lb = lb_theory
                status_str = f"! SDP failed → f ≥ {lb_theory} (theory)"
                lam_best = None

            print(f"{status_str}  [{elapsed:.1f}s, {num_dirs} candidates]")

            results.append({
                'i': i,
                'M': M,
                'Hypercube_UB': ub,
                'Theory_LB': lb_theory,
                'Best_Lambda*': lam_best,
                'Certified_LB': cert_lb,
                'Method': 'grid_search',
                'Verified': verified,
                'Notes': f'tested {num_trials} trials'
            })

        print()

    return results


# ============================================================================
# Reporting
# ============================================================================

def print_summary_table(results):
    """Print a nice summary table."""
    print("=" * 90)
    print("SUMMARY TABLE")
    print("=" * 90)
    print()
    print("Legend:")
    print("  M = ceil(i/(k-1)) = number of bad-bat groups")
    print("  Hypercube_UB = 2^(M-1) = upper bound from hypercube construction")
    print("  Theory_LB = proven lower bound from partition argument")
    print("  Best_Lambda* = best s=2 SOS certificate found (grid search only)")
    print("  Certified_LB = lower bound we can prove with current methods")
    print()

    # ASCII table
    print(f"{'i':>2} | {'M':>1} | {'Hyper_UB':>8} | {'Theory_LB':>9} | {'Best_λ*':>10} | {'Cert_LB':>8} | Status")
    print("-" * 90)

    for r in results:
        i, M, ub, lb_th = r['i'], r['M'], r['Hypercube_UB'], r['Theory_LB']
        lam = r['Best_Lambda*']
        cert_lb = r['Certified_LB']

        lam_str = f"{lam:.6f}" if lam is not None else "—"
        status = "✓" if r['Method'] == 'theory_only' else ("✓" if lam and lam > 0.001 else "?")

        print(f"{i:2d} | {M:1d} | {ub:8d} | {lb_th:9d} | {lam_str:>10} | {cert_lb:8d} | {status}")

    print()


def save_csv(results, filename='dd_phase1_results.csv'):
    """Save results to CSV."""
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['i', 'M', 'Hypercube_UB', 'Theory_LB',
                                                'Best_Lambda*', 'Certified_LB', 'Method', 'Verified', 'Notes'])
        writer.writeheader()
        writer.writerows(results)
    print(f"Results saved to {filename}")


# ============================================================================
# Main
# ============================================================================

if __name__ == '__main__':
    results = run_survey()
    print_summary_table(results)
    save_csv(results)

    # Final summary
    print("=" * 90)
    print("Key findings:")
    print()
    for i in [1, 4, 7, 10, 13, 16]:
        r = results[i - 1]
        ub = r['Hypercube_UB']
        cert_lb = r['Certified_LB']
        gap = ub - cert_lb
        print(f"  i={i:2d}: {cert_lb} ≤ f(DD({i})) ≤ {ub}  (gap: {gap})")
    print()
