# Unit-aware Bursarial embedder ("UnitAwareBursarialEncoder").
#
# Writeup: DOCS/dimensional_multiset_embedding/main.tex
#
# Encodes a length-n multiset of k-vectors whose k components carry FIXED physical
# units, by the "encode by factorise" engine made unit-honest.  Each component c is
# homogenised by a monomial in B base-unit bookkeeping variables xi_0..xi_{B-1}, so
# that the per-vector polynomial
#
#     Psi_j  =  sum_c  v_{j,c} * t^c * prod_b xi_b^{ e_{c,b} }
#
# adds only commensurable terms (all of one common unit -- see the exponent note
# below), never two incommensurable quantities.  The whole multiset becomes
#
#     Phi    =  prod_j ( w - Psi_j ),
#
# and the embedding is the list of coefficients of Phi (a polynomial in the
# B+2 bookkeeping variables w, t, xi_0..xi_{B-1}), dropping only the trivial monic
# leading coefficient (the w^n term, which is always 1).  Injectivity is the
# uniqueness of the factorisation of Phi back into the linear factors (w - Psi_j);
# permutation-invariance is the symmetry of the product; continuity is polynomial.
#
# The component units live in k, not n, so this Embedder is specialised to ONE k,
# fixed at construction by the list of dimension exponent vectors.  size_from_n_k(.)
# returns -1 for any k that does not match.  (n is free.)
#
# Non-negative-exponent convention.  The writeup uses xi^{-d_c}.  Here we use
# xi^{ D - d_c } with D = elementwise max of the d_c, i.e. we multiply every Psi by
# the single fixed monomial xi^{D}.  That is one legal (unit-honest) multiplication:
# it merely gives every term the common unit base^D in place of "dimensionless".  It
# keeps all xi exponents >= 0 so we can use ordinary (non-Laurent) polynomials, and
# it changes neither injectivity nor permutation-invariance nor continuity.

import numpy as np
from math import prod
import sympy
from sympy import symbols, Integer, Poly
from MultisetEmbedder import MultisetEmbedder
from typing import Any


class UnitAwareBursarialEncoder(MultisetEmbedder):

    def __init__(self, dimensions):
        """
        dimensions : a length-k iterable of length-B integer exponent vectors, one
                     per component, in a fixed base-unit order.  Examples:

                       (mass, length, length) over base order (M, L):
                           [(1, 0), (0, 1), (0, 1)]
                       (velocity, energy) over base order (M, L, T):
                           [(0, 1, -1), (1, 2, -2)]

                     Fixes k = len(dimensions) and B = len(dimensions[0]).
        """
        dims = [tuple(int(e) for e in d) for d in dimensions]
        k = len(dims)
        if k == 0:
            raise ValueError("UnitAwareBursarialEncoder needs at least one component (k>=1).")
        B = len(dims[0])
        if any(len(d) != B for d in dims):
            raise ValueError("All dimension exponent vectors must share one length B.")

        self._dims = dims
        self._k = k
        self._B = B
        # Elementwise max of the exponent vectors, used to shift all xi-exponents >= 0.
        self._D = tuple(max(d[b] for d in dims) for b in range(B))

        # Bookkeeping variables, in a fixed generator order:  w, t, xi_0 .. xi_{B-1}.
        self._w = symbols('w')
        self._t = symbols('t')
        self._xis = tuple(symbols(f'xi0:{B}')) if B > 0 else tuple()
        self._gens = (self._w, self._t) + self._xis

        self._support_cache = {}  # n -> sorted list of monomial exponent-tuples

    # ------------------------------------------------------------------ #
    #  Introspection                                                     #
    # ------------------------------------------------------------------ #
    @property
    def k(self):
        return self._k

    @property
    def B(self):
        return self._B

    @property
    def dimensions(self):
        return list(self._dims)

    # ------------------------------------------------------------------ #
    #  Core construction                                                 #
    # ------------------------------------------------------------------ #
    def _address(self, c):
        """The (coefficient-free) monomial  t^c * prod_b xi_b^{D_b - d_{c,b}}  for component c."""
        mono = self._t ** c
        for b in range(self._B):
            e = self._D[b] - self._dims[c][b]
            if e:
                mono *= self._xis[b] ** e
        return mono

    def _psi(self, vec):
        """Psi for one vector: sum_c vec[c] * address(c)."""
        return sum(vec[c] * self._address(c) for c in range(self._k))

    def _phi_poly(self, rows):
        """Phi = prod_j (w - Psi_j) as an expanded sympy Poly in self._gens."""
        expr = prod(self._w - self._psi(rows[j]) for j in range(len(rows)))  # n>=2 so non-empty
        return Poly(sympy.expand(expr), *self._gens)

    def _canonical_support(self, n):
        """
        The data-independent list of monomials Phi can carry, EXCLUDING the trivial
        monic leading term w^n (always 1).  Computed once per n from a symbolic
        template (so it is exactly the generically-present support) and cached.
        """
        if n in self._support_cache:
            return self._support_cache[n]
        a = [[symbols(f'a_{j}_{c}') for c in range(self._k)] for j in range(n)]
        poly = self._phi_poly(a)
        support = sorted(m for m in poly.monoms() if m[0] != n)  # m[0] is the w-power
        self._support_cache[n] = support
        return support

    # ------------------------------------------------------------------ #
    #  Framework interface                                               #
    # ------------------------------------------------------------------ #
    def size_from_n_k(self, n: int, k: int) -> int:
        # Enforce the promised k for EVERY shape (including the n==1 / k==1 edges that
        # the base class would otherwise answer without consulting us).
        if k != self._k:
            return -1
        return super().size_from_n_k(n, k)

    def size_from_n_k_generic(self, n: int, k: int) -> int:
        if k != self._k:
            return -1
        return len(self._canonical_support(n))

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        # k==1: a multiset of reals, all sharing the single component's unit, so the
        # ordinary (Cinf) polynomial embedder is already unit-honest.
        return MultisetEmbedder.embed_kOne_polynomial(data), None

    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        n = len(data)
        k = len(data[0])
        if k != self._k:
            raise ValueError(
                f"This UnitAwareBursarialEncoder was built for k={self._k}, got k={k}.")

        support = self._canonical_support(n)
        rows = [[float(x) for x in row] for row in data]   # plain floats for sympy
        poly = self._phi_poly(rows)
        d = poly.as_dict()  # {exponent_tuple: coefficient}
        coeffs = [float(d.get(mono, Integer(0))) for mono in support]

        expected = self.size_from_n_k(n, k)
        if len(coeffs) != expected:
            raise Exception(
                f"Bug: produced {len(coeffs)} coeffs but size_from_n_k said {expected}.")
        return np.asarray(coeffs, dtype=np.float64), None


# Keep the repo-wide convention that a module exposes `Embedder`.
Embedder = UnitAwareBursarialEncoder
