"""
symcoder.sign_correlation
==========================
Packing helper shared by row-pair encoders.
"""
from __future__ import annotations

import numpy as np


def _complex_to_reals(c: np.ndarray) -> np.ndarray:
    """Unpack n complex values to 2n reals: [re0, im0, re1, im1, ...]."""
    out = np.empty(2 * len(c), dtype=float)
    out[0::2] = c.real
    out[1::2] = c.imag
    return out
