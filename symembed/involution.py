"""
Involution packing (notes, Section 3.3): turning invariance of a multiset under a linear involution into real
coefficients of a complex Bursarial embedding.

Let T be a real linear map with T^2 = I, and let a multiset of real m-vectors be invariant (as a multiset) under T.
Write each vector in coordinates adapted to T: the +1 coordinates as real parts and the -1 coordinates as imaginary
parts of complex numbers (padding the shorter list with zeros). In these coordinates T acts as complex conjugation, so
the multiset of packed vectors is closed under conjugation, and the complex Bursarial polynomial
prod (y + q(x)) has real coefficients. Only those need storing.

The core (``packing_from_matrix``) is generic in T. For a multiset {g.w : g in G}, any g0 in G whose action on w is an
involution provides T, because g -> g0 g permutes G; ``packing_from_signed_involution`` gives the corresponding
integer coordinates explicitly.
"""
from __future__ import annotations

from dataclasses import dataclass

import sympy


@dataclass(frozen=True)
class Packing:
    """
    Packed coordinate j is z_j = sum_s plus[j][s] w_s + i sum_s minus[j][s] w_s (integer coefficients).
    Rows beyond the number of +1 (resp. -1) coordinates are zero (padding).
    """
    plus: tuple            # r rows of m integers
    minus: tuple           # r rows of m integers

    @property
    def r(self) -> int:
        return len(self.plus)

    @property
    def m(self) -> int:
        return len(self.plus[0]) if self.plus else 0

    def describe(self, names: list[str]) -> list[tuple[str, str]]:
        """Each packed coordinate as (real part, imaginary part) in terms of the given slot names."""
        def form(row):
            terms = []
            for coeff, name in zip(row, names):
                if coeff:
                    sign = "-" if coeff < 0 else "+"
                    mag = "" if abs(coeff) == 1 else f"{abs(coeff)}*"
                    terms.append(f"{sign} {mag}{name}")
            text = " ".join(terms)
            return "0" if not text else (text[2:] if text.startswith("+ ") else "-" + text[2:])
        return [(form(p), form(q)) for p, q in zip(self.plus, self.minus)]


def _pad(rows: list, r: int, m: int) -> tuple:
    return tuple(tuple(row) for row in rows) + tuple((0,) * m for _ in range(r - len(rows)))


def packing_from_matrix(T) -> Packing:
    """
    Generic core: integer bases of the +1 and -1 eigenspaces of an integer (or rational) involution T, acting on
    vectors as w -> T w. The coordinate functionals f satisfy f(T w) = +f(w) (real parts) and -f(w) (imaginary parts).
    """
    T = sympy.Matrix(T)
    m = T.shape[0]
    if T.shape != (m, m) or T * T != sympy.eye(m):
        raise ValueError("T must be a square matrix with T^2 = I")

    def integer_rows(space):
        rows = []
        for v in space:
            v = v.T                                              # a functional f with f T = +-f
            denominator = sympy.ilcm(*[sympy.fraction(x)[1] for x in v])
            v = v * denominator
            g = sympy.igcd(*[int(x) for x in v]) or 1
            rows.append([int(x) // g for x in v])
        return rows

    # Functionals f (row vectors) with f T = f, i.e. (T^T - I) f^T = 0, and similarly for -1.
    plus = integer_rows((T.T - sympy.eye(m)).nullspace())
    minus = integer_rows((T.T + sympy.eye(m)).nullspace())
    r = max(len(plus), len(minus))
    return Packing(_pad(plus, r, m), _pad(minus, r, m))


def packing_from_signed_involution(index: list, sign: list) -> Packing:
    """
    Integer packing for a signed-permutation involution, (T w)[s] = sign[s] * w[index[s]], in slot form:
    a fixed slot s contributes w_s to the +1 or -1 part according to sign[s]; a swapped pair (s, t) contributes
    w_s + sign[s] w_t (+1 part) and w_s - sign[s] w_t (-1 part).
    """
    m = len(index)
    for s in range(m):
        if index[index[s]] != s or sign[s] * sign[index[s]] != 1:
            raise ValueError("not an involution")
    plus, minus, done = [], [], set()
    for s in range(m):
        if s in done:
            continue
        t = index[s]
        done.update({s, t})
        if t == s:
            row = [0] * m
            row[s] = 1
            (plus if sign[s] > 0 else minus).append(row)
        else:
            p, q = [0] * m, [0] * m
            p[s], p[t] = 1, sign[s]
            q[s], q[t] = 1, -sign[s]
            plus.append(p)
            minus.append(q)
    r = max(len(plus), len(minus))
    return Packing(_pad(plus, r, m), _pad(minus, r, m))


def is_involution(index: list, sign: list) -> bool:
    return all(index[index[s]] == s and sign[s] * sign[index[s]] == 1 for s in range(len(index)))


def trace(index: list, sign: list) -> int:
    return sum(sign[s] for s in range(len(index)) if index[s] == s)


def choose_involution(index: list, sign: list) -> int | None:
    """
    The group element (row of the signed-permutation table) whose action on w is a non-trivial involution with the
    smallest |trace|, which minimises the number of packed coordinates (m + |trace|)/2. Ties go to the first in
    enumeration order. None if there is no such element.
    """
    best = None
    for g, (row_index, row_sign) in enumerate(zip(index, sign)):
        if all(row_index[s] == s and row_sign[s] == 1 for s in range(len(row_index))):
            continue                                            # acts trivially on w
        if is_involution(row_index, row_sign):
            tr = abs(trace(row_index, row_sign))
            if best is None or tr < best[0]:
                best = (tr, g)
    return None if best is None else best[1]


def packed_forms(index: list, sign: list, packing: Packing) -> list:
    """
    The packed polynomials as linear forms in the slots of w, for ``reductions.exact_expansion_forms``:
    forms[g][j] lists (slot t, (a, b)) such that coefficient j of q_g is sum (a + i b) w_t.
    """
    forms = []
    for row_index, row_sign in zip(index, sign):
        row = []
        for p_row, q_row in zip(packing.plus, packing.minus):
            terms = {}
            for s, (pc, qc) in enumerate(zip(p_row, q_row)):
                if pc or qc:
                    t, sg = row_index[s], row_sign[s]
                    a, b = terms.get(t, (0, 0))
                    terms[t] = (a + pc * sg, b + qc * sg)
            row.append([(t, ab) for t, ab in sorted(terms.items()) if ab != (0, 0)])
        forms.append(row)
    return forms
