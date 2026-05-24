"""
Tests for symcoder.operations — the standard operations library.

Checks three things for each module:
  - symbolic properties (name, rank, odd_parity, argument_symmetry, tex set)
  - numerical correctness of eval_fn on known inputs
  - singleton identity (repeated imports give the same object)
"""
import pytest
import numpy as np
from symatom import ArgumentSymmetry


# ─────────────────────────────────────────────────────────────────────────────
# Euclidean 3-D
# ─────────────────────────────────────────────────────────────────────────────

class TestEuclidean3:

    # -- symbolic properties --------------------------------------------------

    def test_mag_properties(self):
        from symcoder.operations.euclidean3 import mag
        assert mag.name == "mag"
        assert mag.rank == 1
        assert mag.odd_parity is False
        assert mag.argument_symmetry == ArgumentSymmetry.SYMMETRIC
        assert mag.tex is not None

    def test_dot_properties(self):
        from symcoder.operations.euclidean3 import dot
        assert dot.name == "dot"
        assert dot.rank == 2
        assert dot.odd_parity is False
        assert dot.argument_symmetry == ArgumentSymmetry.SYMMETRIC
        assert dot.tex is not None

    def test_eps_properties(self):
        from symcoder.operations.euclidean3 import eps
        assert eps.name == "eps"
        assert eps.rank == 3
        assert eps.odd_parity is True
        assert eps.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        assert eps.tex is not None

    # -- numerical correctness ------------------------------------------------

    def test_mag_eval(self):
        from symcoder.operations.euclidean3 import mag
        assert mag.eval_fn([np.array([3.0, 4.0, 0.0])]) == pytest.approx(5.0)

    def test_dot_eval_orthogonal(self):
        from symcoder.operations.euclidean3 import dot
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 1.0, 0.0])
        assert dot.eval_fn([a, b]) == pytest.approx(0.0)

    def test_dot_eval_parallel(self):
        from symcoder.operations.euclidean3 import dot
        a = np.array([2.0, 0.0, 0.0])
        b = np.array([3.0, 0.0, 0.0])
        assert dot.eval_fn([a, b]) == pytest.approx(6.0)

    def test_eps_eval_right_handed_frame(self):
        from symcoder.operations.euclidean3 import eps
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 1.0, 0.0])
        c = np.array([0.0, 0.0, 1.0])
        assert eps.eval_fn([a, b, c]) == pytest.approx(+1.0)

    def test_eps_eval_left_handed_frame(self):
        from symcoder.operations.euclidean3 import eps
        # swap two args → sign flips
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 1.0, 0.0])
        c = np.array([0.0, 0.0, 1.0])
        assert eps.eval_fn([b, a, c]) == pytest.approx(-1.0)

    # -- singleton identity ---------------------------------------------------

    def test_mag_singleton(self):
        from symcoder.operations.euclidean3 import mag
        from symcoder.operations.euclidean3 import mag as mag2
        assert mag is mag2

    def test_dot_singleton(self):
        from symcoder.operations.euclidean3 import dot
        from symcoder.operations.euclidean3 import dot as dot2
        assert dot is dot2

    def test_eps_singleton(self):
        from symcoder.operations.euclidean3 import eps
        from symcoder.operations.euclidean3 import eps as eps2
        assert eps is eps2


# ─────────────────────────────────────────────────────────────────────────────
# Euclidean 2-D
# ─────────────────────────────────────────────────────────────────────────────

class TestEuclidean2:

    # -- symbolic properties --------------------------------------------------

    def test_mag_properties(self):
        from symcoder.operations.euclidean2 import mag
        assert mag.rank == 1
        assert mag.odd_parity is False
        assert mag.argument_symmetry == ArgumentSymmetry.SYMMETRIC

    def test_dot_properties(self):
        from symcoder.operations.euclidean2 import dot
        assert dot.rank == 2
        assert dot.odd_parity is False
        assert dot.argument_symmetry == ArgumentSymmetry.SYMMETRIC

    def test_eps_properties(self):
        from symcoder.operations.euclidean2 import eps
        assert eps.rank == 2
        assert eps.odd_parity is True
        assert eps.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

    # -- numerical correctness ------------------------------------------------

    def test_mag_eval(self):
        from symcoder.operations.euclidean2 import mag
        assert mag.eval_fn([np.array([3.0, 4.0])]) == pytest.approx(5.0)

    def test_dot_eval(self):
        from symcoder.operations.euclidean2 import dot
        a = np.array([1.0, 2.0])
        b = np.array([3.0, 4.0])
        assert dot.eval_fn([a, b]) == pytest.approx(11.0)

    def test_eps_eval_standard_basis(self):
        from symcoder.operations.euclidean2 import eps
        a = np.array([1.0, 0.0])
        b = np.array([0.0, 1.0])
        assert eps.eval_fn([a, b]) == pytest.approx(+1.0)
        assert eps.eval_fn([b, a]) == pytest.approx(-1.0)

    # -- singleton identity ---------------------------------------------------

    def test_singletons(self):
        from symcoder.operations.euclidean2 import mag, dot, eps
        from symcoder.operations.euclidean2 import mag as mag2, dot as dot2, eps as eps2
        assert mag is mag2
        assert dot is dot2
        assert eps is eps2

    # -- cross-module distinctness --------------------------------------------

    def test_mag_distinct_from_euclidean3(self):
        from symcoder.operations.euclidean2 import mag as mag2
        from symcoder.operations.euclidean3 import mag as mag3
        assert mag2 is not mag3

    def test_dot_distinct_from_euclidean3(self):
        from symcoder.operations.euclidean2 import dot as dot2
        from symcoder.operations.euclidean3 import dot as dot3
        assert dot2 is not dot3

    def test_eps_distinct_from_euclidean3(self):
        from symcoder.operations.euclidean2 import eps as eps2
        from symcoder.operations.euclidean3 import eps as eps3
        assert eps2 is not eps3


# ─────────────────────────────────────────────────────────────────────────────
# Minkowski 3+1-D
# ─────────────────────────────────────────────────────────────────────────────

class TestMinkowski:

    # -- symbolic properties --------------------------------------------------

    def test_dot_properties(self):
        from symcoder.operations.minkowski import dot
        assert dot.rank == 2
        assert dot.odd_parity is False
        assert dot.argument_symmetry == ArgumentSymmetry.SYMMETRIC

    def test_mag_sq_properties(self):
        from symcoder.operations.minkowski import mag_sq
        assert mag_sq.rank == 1
        assert mag_sq.odd_parity is False
        assert mag_sq.argument_symmetry == ArgumentSymmetry.SYMMETRIC

    def test_eps_properties(self):
        from symcoder.operations.minkowski import eps
        assert eps.rank == 4
        assert eps.odd_parity is True
        assert eps.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

    # -- numerical correctness ------------------------------------------------

    def test_dot_timelike(self):
        # two identical timelike unit vectors: p·p = +1 in (+−−−)
        from symcoder.operations.minkowski import dot
        p = np.array([1.0, 0.0, 0.0, 0.0])
        assert dot.eval_fn([p, p]) == pytest.approx(+1.0)

    def test_dot_spacelike(self):
        # two identical spacelike unit vectors: p·p = −1 in (+−−−)
        from symcoder.operations.minkowski import dot
        p = np.array([0.0, 1.0, 0.0, 0.0])
        assert dot.eval_fn([p, p]) == pytest.approx(-1.0)

    def test_dot_orthogonal(self):
        from symcoder.operations.minkowski import dot
        t = np.array([1.0, 0.0, 0.0, 0.0])
        x = np.array([0.0, 1.0, 0.0, 0.0])
        assert dot.eval_fn([t, x]) == pytest.approx(0.0)

    def test_mag_sq_timelike(self):
        from symcoder.operations.minkowski import mag_sq
        p = np.array([2.0, 1.0, 0.0, 0.0])   # p² = 4 − 1 = 3
        assert mag_sq.eval_fn([p]) == pytest.approx(3.0)

    def test_mag_sq_spacelike_negative(self):
        from symcoder.operations.minkowski import mag_sq
        p = np.array([0.0, 0.0, 0.0, 1.0])   # pure spatial → p² = −1
        assert mag_sq.eval_fn([p]) == pytest.approx(-1.0)

    def test_eps_standard_basis(self):
        # det(I_4) = 1
        from symcoder.operations.minkowski import eps
        t = np.array([1.0, 0.0, 0.0, 0.0])
        x = np.array([0.0, 1.0, 0.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 0.0])
        z = np.array([0.0, 0.0, 0.0, 1.0])
        assert eps.eval_fn([t, x, y, z]) == pytest.approx(+1.0)

    def test_eps_swap_flips_sign(self):
        from symcoder.operations.minkowski import eps
        t = np.array([1.0, 0.0, 0.0, 0.0])
        x = np.array([0.0, 1.0, 0.0, 0.0])
        y = np.array([0.0, 0.0, 1.0, 0.0])
        z = np.array([0.0, 0.0, 0.0, 1.0])
        assert eps.eval_fn([x, t, y, z]) == pytest.approx(-1.0)

    # -- singleton identity ---------------------------------------------------

    def test_singletons(self):
        from symcoder.operations.minkowski import dot, mag_sq, eps
        from symcoder.operations.minkowski import dot as dot2, mag_sq as ms2, eps as eps2
        assert dot is dot2
        assert mag_sq is ms2
        assert eps is eps2

    def test_dot_distinct_from_euclidean3(self):
        from symcoder.operations.minkowski import dot as dot_m
        from symcoder.operations.euclidean3 import dot as dot_e
        assert dot_m is not dot_e
