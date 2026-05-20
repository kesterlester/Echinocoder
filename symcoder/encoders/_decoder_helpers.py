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
