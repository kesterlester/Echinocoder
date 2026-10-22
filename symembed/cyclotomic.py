"""
Exact arithmetic in the cyclotomic integers Z[zeta_L], zeta_L = exp(2 pi i / L).

An element a_0 + a_1 zeta + ... + a_{f-1} zeta^{f-1} (f = phi(L), the degree of the cyclotomic polynomial Phi_L) is
stored as an integer vector of length f. Multiplication by a fixed element is an f x f integer matrix acting on such
vectors; complex conjugation (zeta -> zeta^-1) is another. For L divisible by 4 the field contains i = zeta^(L/4).

This is all that graded packing needs to be handled exactly: its coordinates have coefficients in Z[zeta_r] and i,
both inside Z[zeta_L] for L = lcm(r, 4). (L = 4 gives the Gaussian integers.)
"""
from __future__ import annotations

from functools import cached_property
from math import gcd

import numpy as np
import sympy


class Cyclotomic:
    """The ring Z[zeta_L], with elements as integer vectors in the power basis 1, zeta, ..., zeta^(f-1)."""

    def __init__(self, L: int):
        if L < 1:
            raise ValueError("L must be positive")
        self.L = L
        x = sympy.Symbol("x")
        phi = sympy.Poly(sympy.cyclotomic_poly(L, x), x)
        self.f = phi.degree()
        self.phi = [int(c) for c in reversed(phi.all_coeffs())]    # phi[t] = coefficient of x^t; phi[f] = 1

    def zero(self) -> np.ndarray:
        return np.zeros(self.f, dtype=object)

    def integer(self, n: int) -> np.ndarray:
        v = self.zero()
        v[0] = int(n)
        return v

    def zeta_power(self, e: int) -> np.ndarray:
        """zeta^e, reduced modulo Phi_L."""
        v = self.zero()
        v[0] = 1
        for _ in range(e % self.L):
            v = self.times_zeta(v)
        return v

    def times_zeta(self, v: np.ndarray) -> np.ndarray:
        """v * zeta: shift up, then reduce the zeta^f term using Phi_L(zeta) = 0."""
        top = v[-1]
        out = np.concatenate([[0], v[:-1]]).astype(object)
        if top:
            out = out - top * np.array(self.phi[:-1], dtype=object)
        return out

    @cached_property
    def i(self) -> np.ndarray:
        if self.L % 4:
            raise ValueError(f"i is not in Z[zeta_{self.L}]")
        return self.zeta_power(self.L // 4)

    def mul_matrix(self, a: np.ndarray) -> np.ndarray:
        """The f x f integer matrix M with M @ v = a * v (as vectors)."""
        columns, basis = [], self.integer(1)
        for _ in range(self.f):
            columns.append(self.mul(a, basis))
            basis = self.times_zeta(basis)
        return np.array(columns, dtype=object).T

    def mul(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        out, power = self.zero(), np.array(b, dtype=object)
        for coeff in a:
            if coeff:
                out = out + coeff * power
            power = self.times_zeta(power)
        return out

    @cached_property
    def conj_matrix(self) -> np.ndarray:
        """The matrix of complex conjugation zeta -> zeta^-1 (a ring automorphism)."""
        return np.array([self.zeta_power(-t) for t in range(self.f)], dtype=object).T

    def conj(self, v: np.ndarray) -> np.ndarray:
        return self.conj_matrix @ np.array(v, dtype=object)

    @cached_property
    def _powers(self) -> np.ndarray:
        return np.exp(2j * np.pi * np.arange(self.f) / self.L)

    def to_complex(self, v) -> complex | np.ndarray:
        """Numerical value(s); v may have extra leading axes, the last axis being the power basis."""
        return np.asarray(v, dtype=float) @ self._powers

    def to_sympy(self, v) -> sympy.Expr:
        """Exact value as a sympy expression (radicals where sympy can find them)."""
        z = sympy.exp(2 * sympy.pi * sympy.I / self.L)
        return sympy.nsimplify(sympy.expand_complex(sum(int(c) * z ** t for t, c in enumerate(v))))


def field_for(r: int) -> Cyclotomic:
    """The smallest Z[zeta_L] containing zeta_r and i: L = lcm(r, 4)."""
    return Cyclotomic(r * 4 // gcd(r, 4))
