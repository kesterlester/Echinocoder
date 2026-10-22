"""Exact cyclotomic integers Z[zeta_L], checked against complex floating point."""
import numpy as np
import pytest
import sympy

from symembed.cyclotomic import Cyclotomic, field_for


@pytest.mark.parametrize("L", [1, 2, 3, 4, 5, 6, 8, 12, 20])
def test_zeta_powers(L):
    F = Cyclotomic(L)
    assert F.f == sympy.totient(L)
    assert np.array_equal(F.zeta_power(L), F.integer(1))
    for e in range(-L, 2 * L):
        assert np.isclose(F.to_complex(F.zeta_power(e)), np.exp(2j * np.pi * e / L))


@pytest.mark.parametrize("L", [3, 5, 12, 20])
def test_multiplication_and_conjugation(L):
    F = Cyclotomic(L)
    rng = np.random.default_rng(L)
    for _ in range(20):
        a = np.array(rng.integers(-5, 6, F.f), dtype=object)
        b = np.array(rng.integers(-5, 6, F.f), dtype=object)
        assert np.isclose(F.to_complex(F.mul(a, b)), F.to_complex(a) * F.to_complex(b))
        assert np.array_equal(F.mul_matrix(a) @ b, F.mul(a, b))
        assert np.isclose(F.to_complex(F.conj(a)), np.conj(F.to_complex(a)))
        assert np.array_equal(F.conj(F.mul(a, b)), F.mul(F.conj(a), F.conj(b)))     # an automorphism


def test_i_and_field_for():
    for r, L in [(1, 4), (2, 4), (3, 12), (4, 4), (5, 20), (6, 12)]:
        F = field_for(r)
        assert F.L == L
        assert np.isclose(F.to_complex(F.i), 1j)
        assert np.isclose(F.to_complex(F.zeta_power(L // r)), np.exp(2j * np.pi / r))


def test_to_sympy():
    F = Cyclotomic(12)
    assert sympy.simplify(F.to_sympy(F.zeta_power(1)) - (sympy.sqrt(3) / 2 + sympy.I / 2)) == 0
