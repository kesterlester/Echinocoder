"""
symcoder.encoders._decoder_helpers
====================================
Shared inverse operations used by decode() methods on row-pair encoders.

These are pure computational helpers — they do not touch any encoder state
and impose no constraints on the encoding path.
"""
from __future__ import annotations

import numpy as np


def _reals_to_complex(r: np.ndarray) -> np.ndarray:
    """Pack 2n reals [re0, im0, re1, im1, …] → n complex values.

    Inverse of _complex_to_reals in sign_correlation.py.
    """
    return r[0::2] + 1j * r[1::2]


def _zip_decode(coeffs: np.ndarray) -> np.ndarray:
    """Recover the original n complex values from zip_embed's coefficient output.

    zip_embed(z) stores the n non-leading coefficients of the monic degree-n
    polynomial whose roots are {-z_k}, ordered highest-to-lowest degree.
    Inverting: find the roots of that polynomial, then negate.

    coeffs : complex array of length n (as returned by zip_embed, before any
             further packing by _complex_to_reals or symmetry extraction).
    Returns: complex array of length n = the original z values (as a multiset;
             order is arbitrary).
    """
    poly = np.concatenate([[1.0 + 0j], coeffs])   # prepend monic leading coeff
    roots = np.roots(poly)                          # roots are {-z_k}
    return -roots                                   # recover {z_k}


def _polish_component(decoded_vals: list, phase1_vals: list) -> list:
    """Replace decoded_vals with optimally-matched phase1_vals (Hungarian assignment).

    Phase 1 values are more accurate (sort-based) than Phase 2 values
    (polynomial root-finding).  This replaces each noisy decoded value with the
    closest Phase 1 value, using optimal bijective assignment to avoid the
    collapse problem that greedy nearest-neighbour suffers.

    If len(decoded_vals) == k * len(phase1_vals) for integer k > 1, phase1_vals
    are replicated k times before matching.  This handles encoder types where
    each u- (or v-) Phase 1 atom appears k times in the pair orbit (e.g. TYPE_12
    where u is always positive-signed but v varies over ±, doubling the pair count
    relative to the u-orbit size).

    Returns decoded_vals unchanged if sizes are incompatible (no silent crash).
    """
    from scipy.optimize import linear_sum_assignment
    n, m = len(decoded_vals), len(phase1_vals)
    if n == 0 or m == 0 or n % m != 0:
        return list(decoded_vals)
    k = n // m
    p1 = np.array(phase1_vals * k, dtype=float)
    d  = np.array(decoded_vals,    dtype=float)
    cost = (d[:, None] - p1[None, :]) ** 2
    _, col_ind = linear_sum_assignment(cost)
    return [float(p1[j]) for j in col_ind]


def _apply_polish(pairs: list, u_phase1_vals: list, v_phase1_vals: list) -> list:
    """Polish (u, v) pairs by matching u and v components independently to Phase 1 values.

    Independent matching means the permutation used for u need not equal the
    permutation used for v.  This is correct: the output is a multiset of pairs,
    so no per-element u↔v correspondence is claimed.
    """
    if not pairs:
        return pairs
    u_vals = [u for u, v in pairs]
    v_vals = [v for u, v in pairs]
    u_polished = _polish_component(u_vals, u_phase1_vals)
    v_polished = _polish_component(v_vals, v_phase1_vals)
    return list(zip(u_polished, v_polished))
