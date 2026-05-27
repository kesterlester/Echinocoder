"""
symcoder.encoders.phase3_vandermonde
=====================================
Phase 3 encoder: applies a **Vandermonde-power-sums** embedding (or its
sorted-projection cousin) to the full alignment table T of shape
(|G|, |rep|).  Like the simplicial Phase 3 encoder, this is a one-shot
whole-table encoder that produces a permutation-invariant embedding of T
viewed as a multiset of |G| tuples in R^{|rep|}.

See ``DOCS/phase3_vandermonde.pdf`` for a worked-through derivation of the
mathematics; the summary below is just enough to read the code.

Mathematics (summary)
---------------------
For each Vandermonde node ``t_c`` (one per ``c = 0, 1, ..., m_proj - 1``),
define the linear projection

    L_c(x)  =  sum_{j=1}^{k} t_c^{j-1} · x_j

Apply ``L_c`` to every row of T:

    y_{c, i}  =  L_c(T[i, :])           i = 1..n,  c = 0..m_proj-1

so that ``y[:, c]`` is a multiset of n reals (or complex numbers, when t_c
is complex).

The encoding stores, for each c, an n-real encoding of the multiset
``y[:, c]``:

* **mode = "Cinf"** — power-sum mode, **DFT nodes** (complex roots of unity).

  Use ``m_proj = (k-1)·n + 1`` and choose

      t_c  =  exp(2π i c / m_proj),          c = 0, …, m_proj - 1.

  These are the m_proj-th roots of unity.  All distinct (so the Vandermonde
  non-singularity argument for separation still holds), all of modulus 1
  (so no over/underflow at any power), and the implicit Vandermonde matrix
  is exactly the DFT matrix (condition number = 1, perfectly conditioned).

  For each c and p ∈ {1, …, n}, the (complex) power sum

      M_{c, p}  =  sum_i (y_{c, i})^p

  is a polynomial in T with real coefficients, so it is C-infinity in T.
  By construction M_{c, p}(conjugate(t)) = conjugate(M_{c, p}(t)).  We
  exploit this conjugate symmetry to emit a real-valued vector of length
  m_proj · n:

    - At the real-axis node t_0 = +1 we emit ``Re(M_{c=0, p})`` for each p
      (n reals).
    - When m_proj is even, t_{m_proj/2} = -1 is also on the real axis;
      emit ``Re(M_{m_proj/2, p})`` for each p (n more reals).
    - For each upper-half conjugate pair (t_c, conjugate(t_c)) — i.e.
      c = 1, …, ⌊(m_proj-1)/2⌋ — emit ``Re(M_{c, p})`` and ``Im(M_{c, p})``
      for each p (2n reals per pair).

  Total reals: ``m_proj · n``.  Same size as the naive real-node variant
  would be, but no over/underflow and perfect numerical conditioning.

* **mode = "C0"** — sort-encoding, **real-valued evenly spaced nodes**.

  Use ``m_proj`` distinct real nodes evenly spaced in (-1, 1):

      t_c  =  -1 + 2(c + 0.5) / m_proj,    c = 0, …, m_proj - 1.

  All distinct, all bounded.  For each c, emit the ascending sort of
  ``y[:, c]`` (n reals per projection).  Total reals: ``m_proj · n``.

  Continuous but not differentiable (corners at projection-value ties).
  Useful when the dimensionful blow-up of Cinf power sums (which have
  units of ``input^p`` for p = 1..n) is unacceptable downstream: sort-mode
  outputs share units with the inputs.

Both modes are deterministically faithful, by Newton-Girard plus
Vandermonde non-singularity.

Optional scale parameter
------------------------
The Cinf mode's outputs span units ``input^1`` through ``input^n``.  When
the input vectors come from a known energy scale (e.g. typical
momentum ~1 GeV), the magnitudes of M_{c,p} span 12 orders of magnitude
for n=12, which is numerically awkward.  Passing ``scale=λ`` to the encoder
divides each column of T by ``λ ** mass_dimension[j]`` before the
projections are computed, where ``mass_dimension[j]`` is the homogeneity
degree of the j-th rep atom's operation.  This produces a *dimensionless*
table whose power sums are all of order unity when λ matches the event
scale.  The scale cannot be derived from the data (which can in principle
be zero); it must be supplied as a fixed encoder parameter.

Output size
-----------
For T of shape (n, k), both modes:

    m_proj = (k - 1) · n + 1     (when n > 1 and k > 1)
    output_dim = m_proj · n

For the parity-fixture plan (n = |G| = 12, k = |rep| = 25):

    m_proj = 24 · 12 + 1 = 289
    output_dim = 289 · 12 = 3468 reals

Larger than the simplicial encoder's 576 reals.  Trade-off: smoothness
(Cinf) or unit-homogeneity (C0), and a one-line proof of correctness.
"""
from __future__ import annotations

import numpy as np

from symcoder.encoders._base import EncodingResult
from symcoder.encoders.phase3_common import (
    apply_scale, build_alignment_table, column_mass_dimensions, rep_atoms_for,
)


# Lazy access to symcoder.describe types to avoid import cycles.
def _Phase3Tree_class():
    from symcoder.describe import Phase3Tree
    return Phase3Tree


def _SegmentInfo_class():
    from symcoder.describe import SegmentInfo
    return SegmentInfo


_VALID_MODES = ("Cinf", "C0")


class Phase3VandermondeEncoder:
    """Vandermonde-projection Phase 3 encoder.

    Parameters
    ----------
    plan
        The plan whose ``context.the_group`` and operations determine the
        alignment table shape.
    mode
        Either ``"Cinf"`` (DFT-node polynomial power sums; smooth
        everywhere) or ``"C0"`` (real-node sorted projections; continuous
        but not differentiable at ties).  See the module docstring for the
        mathematical justification.
    scale : float | None
        Optional energy / length scale.  When not None, each column j of
        the alignment table T is divided by ``scale ** mass_dimension_j``
        before the Vandermonde projections are computed, where
        ``mass_dimension_j`` is the homogeneity degree of the j-th rep
        atom's operation.  This makes the per-power-sum magnitudes
        approximately scale-free in Cinf mode.  Has no effect on the
        *faithfulness* of the encoding (the scaling is invertible) — only
        on its numerical conditioning.  ``None`` (the default) skips the
        scaling step entirely.

        Setting ``scale`` requires every operation in the plan to declare
        a ``mass_dimension`` on its ``Operation`` constructor; an
        informative ``ValueError`` is raised at ``encode()`` time if not.
    """

    def __init__(self, plan, mode: str = "Cinf", scale: float | None = None) -> None:
        if mode not in _VALID_MODES:
            raise ValueError(
                f"Phase3VandermondeEncoder: mode must be one of {_VALID_MODES}, "
                f"got {mode!r}"
            )
        if scale is not None and not (isinstance(scale, (int, float)) and scale > 0):
            raise ValueError(
                f"Phase3VandermondeEncoder: scale must be a positive real or None, "
                f"got {scale!r}"
            )
        self._plan       = plan
        self._mode       = mode
        self._scale      = None if scale is None else float(scale)
        self._the_group  = plan.context.the_group
        self._rep_atoms  = rep_atoms_for(plan)
        self._n          = self._the_group.order()
        self._k          = len(self._rep_atoms)

        # If the user wants dimensional rescaling, every operation in the plan
        # must have declared a mass_dimension.  Validate eagerly so the error
        # surfaces at factory build / encoder construction rather than at the
        # first encode() call.  Note: assessment in the factory (assess())
        # already filters out plans that cannot satisfy the scale request, so
        # this assertion-style check is the "defensive layer" — callers that
        # construct Phase3VandermondeEncoder directly with a scale and an
        # incompatible plan get a clear immediate diagnostic.
        if self._scale is not None and self._k > 0:
            column_mass_dimensions(self._rep_atoms)   # raises ValueError on missing

        if self._n == 0 or self._k == 0:
            self._m_proj     = 0
            self._output_dim = 0
        elif self._k == 1:
            # Degenerate k=1 case: the alignment table is already a multiset
            # of n reals.  One projection (L_c(x) = x_0) suffices; that yields
            # n power sums (Cinf) or the sort (C0).
            self._m_proj     = 1
            self._output_dim = self._n
        else:
            # General case.  (k-1)·n + 1 Vandermonde nodes recover all
            # multi-symmetric power sums of total degree ≤ n.
            self._m_proj     = (self._k - 1) * self._n + 1
            self._output_dim = self._m_proj * self._n

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @property
    def output_dim(self) -> int:
        return self._output_dim

    @property
    def mode(self) -> str:
        return self._mode

    @property
    def scale(self) -> float | None:
        return self._scale

    def encode(self, event: dict) -> EncodingResult:
        n, k = self._n, self._k
        if n == 0 or k == 0:
            return EncodingResult(values=np.array([], dtype=np.float64))

        T = build_alignment_table(self._plan, event, self._the_group, self._rep_atoms)
        # Optional dimensional rescaling (no-op when self._scale is None).
        T = apply_scale(T, self._rep_atoms, self._scale)

        if self._mode == "Cinf":
            values = self._encode_dft_power_sums(T)
        else:  # "C0"
            values = self._encode_sorted(T)

        assert values.shape == (self._output_dim,), (
            f"Phase3VandermondeEncoder({self._mode}): produced {values.shape[0]} "
            f"reals; expected {self._output_dim}."
        )
        return EncodingResult(values=values)

    def describe(self, start_offset: int = 0):
        SegmentInfo = _SegmentInfo_class()
        Phase3Tree  = _Phase3Tree_class()
        if self._output_dim == 0:
            return Phase3Tree(segments=[], n=self._n, k=self._k)
        method_name = (
            "vandermonde_cinf_dft" if self._mode == "Cinf" else "vandermonde_c0_real"
        )
        scale_str = "none" if self._scale is None else f"{self._scale:g}"
        seg = SegmentInfo(
            kind            = "VANDERMONDE",
            start           = start_offset,
            length          = self._output_dim,
            op_u            = None,
            flavour_u       = None,
            op_v            = None,
            flavour_v       = None,
            overlap         = None,
            symmetry_class  = None,
            sign_compressed = None,
            method_name     = method_name,
            example         = (
                f"Vandermonde ({self._mode}) embedding of full rep orbit table "
                f"(n=|G|={self._n}, k=|rep|={self._k}, m_proj={self._m_proj}, "
                f"scale={scale_str})"
            ),
        )
        return Phase3Tree(segments=[seg], n=self._n, k=self._k)

    # ------------------------------------------------------------------
    # Cinf internals: DFT-node power sums
    # ------------------------------------------------------------------

    def _encode_dft_power_sums(self, T: np.ndarray) -> np.ndarray:
        """Cinf: per-projection power sums M_{c, p} for p=1..n at DFT nodes,
        emitted as a real-valued vector via conjugate-pair compression.

        Layout of the output (length m_proj · n):

          - First block of n reals: ``Re(M_{0, p})`` for p = 1..n.
            (c = 0 corresponds to t_0 = +1, the real-axis node — and
            the imaginary part is identically zero.)
          - For each c = 1, …, ⌊(m_proj-1)/2⌋ (upper-half nodes), 2n reals:
            ``Re(M_{c, 1}), Im(M_{c, 1}), Re(M_{c, 2}), Im(M_{c, 2}), …``
          - If m_proj is even, a final block of n reals:
            ``Re(M_{m_proj/2, p})`` for p = 1..n.  (c = m_proj/2 is the
            second real-axis node t = -1; its imaginary part is also zero.)

        The lower-half nodes (c > m_proj/2) are NOT stored: their values
        are the complex conjugates of the upper-half ones (because
        M_p has real coefficients) and therefore carry no new information.
        """
        n, k = self._n, self._k

        # m_proj-th roots of unity: t_c = exp(2π i c / m_proj)
        m_proj = self._m_proj
        c_idx  = np.arange(m_proj, dtype=np.int64)
        ts     = np.exp(2j * np.pi * c_idx / m_proj)             # (m_proj,)

        # Projection matrix V of shape (m_proj, k): V[c, j] = t_c^j  (j = 0..k-1)
        exponents = np.arange(k, dtype=np.int64)
        V = ts[:, None] ** exponents[None, :]                    # complex128

        # y[i, c] = sum_j V[c, j] · T[i, j]
        y = T.astype(np.complex128) @ V.T                        # (n, m_proj)

        # Power sums M[c, p] = sum_i y[i, c]^p for p = 1..n
        ps    = np.arange(1, n + 1, dtype=np.int64)              # 1..n
        y_pow = y[:, :, None] ** ps[None, None, :]               # (n, m_proj, n)
        M     = y_pow.sum(axis=0)                                # (m_proj, n) complex

        # Conjugate-pair compression to real-valued output.
        return _conjugate_compress_to_real(M, m_proj, n)

    # ------------------------------------------------------------------
    # C0 internals: real-node sort encoding
    # ------------------------------------------------------------------

    def _encode_sorted(self, T: np.ndarray) -> np.ndarray:
        """C0: per-projection ascending sort of y[:, c] at evenly-spaced
        real nodes in (-1, 1).

        Output is m_proj · n reals, laid out as
        ``[sorted(y[:,0]), sorted(y[:,1]), …]`` (all n entries for c=0 first,
        then c=1, …).  Real-valued nodes are used here (rather than DFT)
        because sorting is undefined on complex numbers and unnecessary for
        bounding the output magnitudes (sort doesn't amplify scale).
        """
        n, k = self._n, self._k
        m_proj = self._m_proj

        # Evenly spaced real nodes in (-1, 1), distinct by construction.
        ts = (-1.0 + 2.0 * (np.arange(m_proj) + 0.5) / m_proj).astype(np.float64)

        exponents = np.arange(k, dtype=np.int64)
        V = ts[:, None] ** exponents[None, :]                    # (m_proj, k)
        y = T @ V.T                                              # (n, m_proj)
        y_sorted = np.sort(y, axis=0)                            # ascending per col
        return y_sorted.T.reshape(-1)                            # c-major flatten


def _conjugate_compress_to_real(M: np.ndarray, m_proj: int, n: int) -> np.ndarray:
    """Compress a complex (m_proj, n) array M to a real vector of length
    ``m_proj * n``, exploiting conjugate symmetry: rows c and m_proj-c
    are complex conjugates of each other (because M_p has real coefficients).

    Layout (matches the docstring on ``_encode_dft_power_sums``):

      - n reals for c = 0 (the real-axis node t = +1):  Re(M[0, :])
      - For each c = 1, …, ⌊(m_proj-1)/2⌋ (upper-half nodes):
            interleaved Re and Im of M[c, :], 2n reals total per c
      - If m_proj is even, n more reals for c = m_proj/2 (t = -1):
            Re(M[m_proj/2, :])

    Sanity-checked here: the imaginary parts at the real-axis nodes are
    asserted to be zero modulo float noise.
    """
    assert M.shape == (m_proj, n)
    parts = []

    # c = 0: real-axis node t = +1
    re0 = M[0].real
    im0 = M[0].imag
    assert np.max(np.abs(im0)) < 1e-9, (
        f"DFT power-sum at t = +1 has non-trivial imaginary part "
        f"(max |im| = {np.max(np.abs(im0)):.3e}); this should be zero "
        f"by symmetry of M_p over R."
    )
    parts.append(re0)                                            # n reals

    # Upper-half conjugate-pair nodes (c = 1 .. floor((m_proj-1)/2))
    n_pairs = (m_proj - 1) // 2
    for c in range(1, n_pairs + 1):
        row = M[c]                                               # complex (n,)
        # Interleave Re and Im so that ``M_{c, p}`` lives at indices
        # 2(p-1), 2(p-1)+1 within the block.
        interleaved = np.empty(2 * n, dtype=np.float64)
        interleaved[0::2] = row.real
        interleaved[1::2] = row.imag
        parts.append(interleaved)                                # 2n reals each

    # If m_proj is even, the final real-axis node is t = -1 at c = m_proj/2.
    if m_proj % 2 == 0:
        c_neg1 = m_proj // 2
        re_neg = M[c_neg1].real
        im_neg = M[c_neg1].imag
        assert np.max(np.abs(im_neg)) < 1e-9, (
            f"DFT power-sum at t = -1 has non-trivial imaginary part "
            f"(max |im| = {np.max(np.abs(im_neg)):.3e}); this should be "
            f"zero by symmetry of M_p over R."
        )
        parts.append(re_neg)                                     # n reals

    out = np.concatenate(parts)
    # Sanity: total = n (c=0) + 2n * n_pairs + (n if m_proj even else 0)
    #               = n * (1 + 2·n_pairs + [m_proj even ? 1 : 0])
    #               = n * m_proj
    assert out.shape == (m_proj * n,), (
        f"_conjugate_compress_to_real: got {out.shape[0]} reals, "
        f"expected {m_proj * n}."
    )
    return out


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


class Phase3VandermondeEncoderFactory:
    """Factory for ``Phase3VandermondeEncoder``.

    The ``mode`` and (optional) ``scale`` arguments are fixed at factory
    construction time.  To materialise multiple variants, instantiate
    multiple factories::

        Phase3VandermondeEncoderFactory(mode="Cinf", scale=1.0)
        Phase3VandermondeEncoderFactory(mode="C0")

    and feed them all into a stacking ``Phase3EncoderFactory``.

    Plan compatibility (``assess``)
    -------------------------------
    When ``scale`` is set, every operation in the plan must declare a
    ``mass_dimension``.  Plans containing an operation without one cannot
    be dimensionally rescaled — for such a plan, ``assess(plan)`` returns
    the empty list (i.e. this factory refuses to offer a bound encoder),
    forcing the caller to either drop the scale request or extend their
    operations with mass dimensions.  When ``scale`` is ``None`` no
    rescaling is performed and ``assess(plan)`` always returns a singleton
    bound encoder.
    """

    def __init__(self, mode: str = "Cinf", scale: float | None = None) -> None:
        if mode not in _VALID_MODES:
            raise ValueError(
                f"Phase3VandermondeEncoderFactory: mode must be one of "
                f"{_VALID_MODES}, got {mode!r}"
            )
        if scale is not None and not (isinstance(scale, (int, float)) and scale > 0):
            raise ValueError(
                f"Phase3VandermondeEncoderFactory: scale must be a positive "
                f"real or None, got {scale!r}"
            )
        self._mode  = mode
        self._scale = None if scale is None else float(scale)

    def assess(self, plan) -> list:
        """Return ``[bound_encoder]`` or ``[]`` depending on plan compatibility.

        Empty when ``scale`` is set and the plan contains any operation
        lacking ``mass_dimension`` — i.e. when there is no defined power of
        the scale to divide out of that operation's eval, and so the
        encoder cannot honour the scale request.  Non-empty otherwise.
        """
        if self._scale is not None:
            # Check up-front whether all rep atoms have mass_dimension.  If
            # not, refuse to offer a bound encoder.
            atoms = rep_atoms_for(plan)
            try:
                column_mass_dimensions(atoms)
            except ValueError:
                return []
        return [Phase3VandermondeEncoder(plan, mode=self._mode, scale=self._scale)]

    def build(self, plan) -> Phase3VandermondeEncoder:
        """Build a single bound encoder, or raise if the plan is incompatible.

        Defined in terms of ``assess`` so the two APIs cannot drift apart.
        """
        offered = self.assess(plan)
        if not offered:
            raise RuntimeError(
                f"Phase3VandermondeEncoderFactory(mode={self._mode!r}, "
                f"scale={self._scale!r}).build: plan contains operation(s) "
                f"without mass_dimension; cannot honour the scale request.  "
                f"Either drop scale (pass scale=None) or declare a "
                f"mass_dimension on every operation in the plan."
            )
        return offered[0]
