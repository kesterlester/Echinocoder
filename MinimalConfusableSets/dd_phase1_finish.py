"""
dd_phase1_finish.py

Completes the Phase 1 survey for i=14, 15, 16 (k=4).
Appends results to dd_phase1_results.txt and writes a fresh CSV for all i.

Usage:
    python dd_phase1_finish.py
"""

import sys
import os
import time
import csv

# Ensure imports work from this directory
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from sympy import Matrix
import confusable_dd_search as cds

# ============================================================================
# Configuration (must match dd_phase1_survey.py)
# ============================================================================

K = 4
MAX_NORM = 12

NUM_TRIALS = {
    14: 700,
    15: 700,
    16: 800,
}

# ============================================================================
# Helpers (copied from dd_phase1_survey.py)
# ============================================================================

def theory_lower_bound(i, k=4):
    if i < k:
        return 1
    elif i <= 2*k - 2:
        return 2
    else:
        return 3

def hypercube_upper_bound(i, k=4):
    M = (i + k - 2) // (k - 1)
    return 2 ** (M - 1)

def search_best_directions(i, k=4, num_trials=500, max_norm=12, verbose=True):
    if verbose:
        print(f"  Generating candidate directions (max_norm={max_norm})...", end='', flush=True)

    # Enumerate as (float_vec, canonical_int_vec) pairs for exact verification
    directions_both = list(cds.enumerate_directions(k, max_norm, use_float=True,
                                                    lazy=False, return_int_vecs=True))
    directions_float = [d for d, _ in directions_both]

    if verbose:
        print(f" {len(directions_float)} candidates")
        print(f"  Searching {num_trials} random {i}-tuples...", end='', flush=True)

    best_lambda   = -np.inf
    best_verified = False
    best_Q        = None
    best_int_vecs = None

    import random
    random.seed(42 + i)  # Same seed as original survey

    for trial in range(num_trials):
        D_indices = [random.randint(0, len(directions_float) - 1) for _ in range(i)]
        D     = [directions_float[idx] for idx in D_indices]
        D_sym = [Matrix(d) for d in D]

        lam, Q, status = cds.sos_sdp_for_directions(D_sym, k=k, s=2, verbose=False)

        if lam is None:
            continue

        if lam > best_lambda:
            best_lambda   = lam
            best_Q        = Q
            best_int_vecs = [directions_both[idx][1] for idx in D_indices]

        num_tried = trial + 1
        if verbose and num_tried % max(50, num_trials // 10) == 0:
            print(f"\r  Searching {num_trials} random {i}-tuples... {num_tried}/{num_trials}",
                  end='', flush=True)

    if verbose:
        print()

    # Verify the best candidate using exact rational arithmetic
    if best_lambda > 0.001 and best_int_vecs is not None and best_Q is not None:
        best_verified, _ = cds.verify_sos_certificate_exact(
            best_int_vecs, best_lambda, best_Q, k=k, verbose=False)

    return best_lambda if best_lambda > -np.inf else None, best_verified, len(directions_float)

# ============================================================================
# Previously completed results (i=1..13) for the final CSV
# ============================================================================

PREVIOUS_RESULTS = [
    # i, M, Hypercube_UB, Theory_LB, Best_Lambda*, Certified_LB, Method, Verified, Notes
    {'i': 1,  'M': 1, 'Hypercube_UB': 1,   'Theory_LB': 1, 'Best_Lambda*': None,   'Certified_LB': 1, 'Method': 'theory_only', 'Verified': 'N/A', 'Notes': 'No search needed'},
    {'i': 2,  'M': 1, 'Hypercube_UB': 1,   'Theory_LB': 1, 'Best_Lambda*': None,   'Certified_LB': 1, 'Method': 'theory_only', 'Verified': 'N/A', 'Notes': 'No search needed'},
    {'i': 3,  'M': 1, 'Hypercube_UB': 1,   'Theory_LB': 1, 'Best_Lambda*': None,   'Certified_LB': 1, 'Method': 'theory_only', 'Verified': 'N/A', 'Notes': 'No search needed'},
    {'i': 4,  'M': 2, 'Hypercube_UB': 2,   'Theory_LB': 2, 'Best_Lambda*': None,   'Certified_LB': 2, 'Method': 'theory_only', 'Verified': 'N/A', 'Notes': 'No search needed'},
    {'i': 5,  'M': 2, 'Hypercube_UB': 2,   'Theory_LB': 2, 'Best_Lambda*': None,   'Certified_LB': 2, 'Method': 'theory_only', 'Verified': 'N/A', 'Notes': 'No search needed'},
    {'i': 6,  'M': 2, 'Hypercube_UB': 2,   'Theory_LB': 2, 'Best_Lambda*': None,   'Certified_LB': 2, 'Method': 'theory_only', 'Verified': 'N/A', 'Notes': 'No search needed'},
    {'i': 7,  'M': 3, 'Hypercube_UB': 4,   'Theory_LB': 3, 'Best_Lambda*': 0.0,    'Certified_LB': 3, 'Method': 'grid_search', 'Verified': False, 'Notes': 'tested 500 trials; lambda*≈0 inconclusive'},
    {'i': 8,  'M': 3, 'Hypercube_UB': 4,   'Theory_LB': 3, 'Best_Lambda*': 0.0188, 'Certified_LB': 3, 'Method': 'grid_search', 'Verified': True,  'Notes': 'tested 500 trials'},
    {'i': 9,  'M': 3, 'Hypercube_UB': 4,   'Theory_LB': 3, 'Best_Lambda*': 0.0219, 'Certified_LB': 3, 'Method': 'grid_search', 'Verified': True,  'Notes': 'tested 500 trials'},
    {'i': 10, 'M': 4, 'Hypercube_UB': 8,   'Theory_LB': 3, 'Best_Lambda*': 0.0575, 'Certified_LB': 3, 'Method': 'grid_search', 'Verified': True,  'Notes': 'tested 600 trials'},
    {'i': 11, 'M': 4, 'Hypercube_UB': 8,   'Theory_LB': 3, 'Best_Lambda*': 0.0737, 'Certified_LB': 3, 'Method': 'grid_search', 'Verified': True,  'Notes': 'tested 600 trials'},
    {'i': 12, 'M': 4, 'Hypercube_UB': 8,   'Theory_LB': 3, 'Best_Lambda*': 0.1233, 'Certified_LB': 3, 'Method': 'grid_search', 'Verified': True,  'Notes': 'tested 600 trials'},
    {'i': 13, 'M': 5, 'Hypercube_UB': 16,  'Theory_LB': 3, 'Best_Lambda*': 0.1525, 'Certified_LB': 3, 'Method': 'grid_search', 'Verified': True,  'Notes': 'tested 700 trials'},
]

# ============================================================================
# Main
# ============================================================================

def run_finish():
    new_results = []

    with open('dd_phase1_results.txt', 'a') as log:
        log.write('\n--- Resuming: running i=14,15,16 ---\n\n')
        log.flush()

        for i in [14, 15, 16]:
            M = (i + K - 2) // (K - 1)
            ub = hypercube_upper_bound(i, K)
            lb_theory = theory_lower_bound(i, K)
            num_trials = NUM_TRIALS[i]

            line_prefix = f"i={i:2d} (M={M}): "
            print(line_prefix, end='', flush=True)
            log.write(line_prefix)
            log.flush()

            t0 = time.time()
            lam_best, verified, num_dirs = search_best_directions(i, K, num_trials, MAX_NORM, verbose=True)
            elapsed = time.time() - t0

            if lam_best is not None and lam_best > 0.001:
                cert_lb = 3
                status_str = f"✓ λ*={lam_best:.4f} → certified f ≥ 3"
            else:
                cert_lb = lb_theory
                lam_val = lam_best if lam_best is not None else 0.0
                status_str = f"✗ λ*≈0 → inconclusive (f ≥ {lb_theory})"
                lam_best = lam_val

            full_line = f"{status_str}  [{elapsed:.1f}s, {num_dirs} candidates]\n"
            print(full_line, end='', flush=True)
            log.write(full_line)
            log.write('\n')
            log.flush()

            new_results.append({
                'i': i,
                'M': M,
                'Hypercube_UB': ub,
                'Theory_LB': lb_theory,
                'Best_Lambda*': lam_best,
                'Certified_LB': cert_lb,
                'Method': 'grid_search',
                'Verified': verified,
                'Notes': f'tested {num_trials} trials',
            })

    # Write complete CSV (all i=1..16)
    all_results = PREVIOUS_RESULTS + new_results
    csv_path = 'dd_phase1_results.csv'
    fieldnames = ['i', 'M', 'Hypercube_UB', 'Theory_LB', 'Best_Lambda*', 'Certified_LB', 'Method', 'Verified', 'Notes']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nComplete results written to {csv_path}")

    # Print summary table
    print()
    print("=" * 90)
    print("COMPLETE SUMMARY TABLE (i=1..16)")
    print("=" * 90)
    print(f"{'i':>2} | {'M':>1} | {'Hyper_UB':>8} | {'Theory_LB':>9} | {'Best_λ*':>10} | {'Cert_LB':>8} | Status")
    print("-" * 90)
    for r in all_results:
        lam = r['Best_Lambda*']
        lam_str = f"{lam:.6f}" if lam is not None else "—"
        status = "✓" if r['Method'] == 'theory_only' else ("✓" if lam and lam > 0.001 else "?")
        print(f"{r['i']:2d} | {r['M']:1d} | {r['Hypercube_UB']:8d} | {r['Theory_LB']:9d} | {lam_str:>10} | {r['Certified_LB']:8d} | {status}")
    print()

    return all_results


if __name__ == '__main__':
    run_finish()
