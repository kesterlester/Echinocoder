"""
Safe (Level 1) reductions of the Bursarial outputs: dropping outputs that are exact linear combinations (over Q) of
outputs that are kept, as polynomials in the slots of w. See Sections 4.2 and 4.3 of the notes.

The E_{c,j} are expanded exactly, as integer polynomials in w, but only for c <= c_max, because the number of
monomials of degree c grows quickly. Outputs with c > c_max are never dropped, which is always safe.

Candidate dependencies are found modulo a prime, reconstructed as rationals, and then verified exactly in integer
arithmetic. A candidate that fails reconstruction or verification is discarded and its output is kept, so a
mistake can only ever cost compression, never injectivity.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb, gcd

import numpy as np

PRIME = 2_147_483_629          # < 2^31, so products of residues fit in int64
_INT64_SAFE = 2 ** 62


def monomials(m: int, c: int) -> list[tuple]:
    """Monomials of degree c in m slot variables, as non-decreasing tuples of slot indices."""
    return list(combinations_with_replacement(range(m), c))


def expansion_cost(m: int, c: int) -> int:
    """Number of integer coefficients in the exact expansion of E_c: (#monomials) x (#powers of x)."""
    return comb(m + c - 1, c) * (c * (m - 1) + 1)


def choose_c_max(N: int, m: int, budget: int) -> int:
    """The largest c_max <= N whose exact expansions E_1..E_{c_max} fit within the budget."""
    total, c_max = 0, 0
    for c in range(1, N + 1):
        total += expansion_cost(m, c)
        if total > budget:
            break
        c_max = c
    return c_max


def exact_expansion(index: list, sign: list, m: int, c_max: int) -> list:
    """
    The exact expansion of E_1..E_{c_max} as integer polynomials in w.

    Returns a list, for c = 1..c_max, of (monos, coeffs), where ``coeffs[j, t]`` is the coefficient of monomial
    ``monos[t]`` in E_{c,j}. Uses the recurrence E_c += q_g * E_{c-1} over g, with q_g(x) = sum_j' sign*w_index x^j'.
    """
    N = len(index)
    # |coefficient| <= C(N, c) m^c; use exact Python ints if that might overflow int64.
    big = any(comb(N, c) * m ** c >= _INT64_SAFE for c in range(1, c_max + 1))
    dtype = object if big else np.int64

    monos = [monomials(m, c) for c in range(c_max + 1)]
    where = [{mono: t for t, mono in enumerate(ms)} for ms in monos]
    # times[c][t, s] = index (in degree c) of monos[c-1][t] * w_s
    times = [None] + [np.array([[where[c][tuple(sorted(mono + (s,)))] for s in range(m)] for mono in monos[c - 1]],
                               dtype=np.int64).reshape(len(monos[c - 1]), m)
                      for c in range(1, c_max + 1)]
    E = [np.ones((1, 1), dtype=dtype)] + [np.zeros((len(monos[c]), c * (m - 1) + 1), dtype=dtype)
                                          for c in range(1, c_max + 1)]
    for g in range(N):
        for c in range(min(g + 1, c_max), 0, -1):          # downwards, so E[c-1] is still the old value
            prev = E[c - 1]
            Lp = prev.shape[1]
            for jp in range(m):
                s, sgn = index[g][jp], sign[g][jp]
                E[c][times[c][:, s], jp:jp + Lp] += sgn * prev
    return [(monos[c], E[c].T.copy()) for c in range(1, c_max + 1)]


def _rational_reconstruct(a: int, p: int = PRIME) -> Fraction | None:
    """The fraction r/s with |r|, |s| < sqrt(p/2) congruent to a mod p, if one exists."""
    a %= p
    bound = int((p // 2) ** 0.5)
    r0, r1, s0, s1 = p, a, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or gcd(r1, abs(s1)) != 1:
        return None
    return Fraction(r1, s1)


@dataclass(frozen=True)
class LinearRelation:
    """Output ``dropped`` (a row index) equals sum of coefficient * row, over the kept rows given."""
    dropped: int
    terms: tuple            # tuple of (row index, Fraction)


def find_linear_relations(coeffs: np.ndarray) -> tuple[list[int], list[LinearRelation]]:
    """
    Greedy, in row order: keep a row unless it is an exact rational linear combination of rows already kept.

    Parameters
    ----------
    coeffs
        Integer matrix (rows = outputs, columns = monomials).

    Returns
    -------
    (kept row indices, relations for the dropped rows). Every relation has been verified exactly.
    """
    rows, cols = coeffs.shape
    if rows == 0:
        return [], []
    if cols == 0:
        return [], [LinearRelation(i, ()) for i in range(rows)]
    A = (np.array(coeffs, dtype=object) % PRIME).astype(np.int64)
    # Augment with the identity to track which combination of original rows each reduced row is.
    aug = np.concatenate([A, np.eye(rows, dtype=np.int64)], axis=1)
    basis: list[tuple[int, np.ndarray]] = []              # (pivot column, row with pivot entry 1)
    kept: list[int] = []
    relations: list[LinearRelation] = []
    for i in range(rows):
        r = aug[i].copy()
        for pivot, b in basis:
            f = int(r[pivot])
            if f:
                r = (r - f * b) % PRIME
        nz = np.nonzero(r[:cols])[0]
        if nz.size:
            pivot = int(nz[0])
            r = (r * pow(int(r[pivot]), PRIME - 2, PRIME)) % PRIME
            basis.append((pivot, r))
            kept.append(i)
            continue
        # sum_t r[cols+t] * row_t == 0 (mod p), with r[cols+i] == 1, so row_i == -sum_{t != i} r[cols+t] row_t.
        terms = []
        for t in kept:
            residue = int(r[cols + t])
            if residue:
                frac = _rational_reconstruct(-residue)
                if frac is None:
                    break
                terms.append((t, frac))
        else:
            if _verify(coeffs, i, terms):
                relations.append(LinearRelation(i, tuple(terms)))
                continue
        # Reconstruction or verification failed (an unlucky prime, or coefficients too large to reconstruct): keep
        # the row, which is always safe. It cannot join the modular basis (it reduced to zero mod p), so later rows
        # that depend on it will also fail verification and be kept.
        kept.append(i)
    return kept, relations


def _verify(coeffs: np.ndarray, i: int, terms: list) -> bool:
    """Exact check, in Python integers, that row i == sum of fraction * row t."""
    denominator = 1
    for _, frac in terms:
        denominator = denominator * frac.denominator // gcd(denominator, frac.denominator)
    lhs = np.array(coeffs[i], dtype=object) * denominator
    rhs = np.zeros(coeffs.shape[1], dtype=object)
    for t, frac in terms:
        rhs = rhs + np.array(coeffs[t], dtype=object) * (frac.numerator * (denominator // frac.denominator))
    return bool(np.all(lhs == rhs))
