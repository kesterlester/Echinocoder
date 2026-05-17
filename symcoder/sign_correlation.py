"""
symcoder.sign_correlation
==========================
SignCorrelationType enum and the canonical algorithm for computing it from a
PairFlavour, plus the _complex_to_reals packing helper.

These live here (not in encode.py) so that both encode.py (legacy path) and
the row_pair_encoders (registry path) can share them without circular imports.
"""
from __future__ import annotations

from enum import Enum

import numpy as np
from symatom.atoms import ArgumentSymmetry


class SignCorrelationType(Enum):
    """
    Describes how the G-orbit of an atom-pair (u, v) relates the sign degrees
    of freedom of u and v.

    The achievable set is the subgroup of Z_2 × Z_2 consisting of all
    (sign_u, sign_v) pairs that appear in the orbit.  There are exactly five
    subgroups of Z_2 × Z_2, corresponding to the five types below.

    TYPE_11   Achievable = {(+1,+1)}.  Neither sign changes.
    TYPE_NEG  Achievable = {(+1,+1),(-1,-1)}.  Only correlated negation.
    TYPE_12   Achievable = {(+1,+1),(+1,-1)}.  Only v's sign flips freely.
    TYPE_21   Achievable = {(+1,+1),(-1,+1)}.  Only u's sign flips freely.
    TYPE_22   Achievable = {(+1,+1),(-1,+1),(+1,-1),(-1,-1)}.  Both flip freely.

    CRITICAL: the type is NOT determined by ArgumentSymmetry alone.
    Always compute via _sign_correlation_type_from_pf(), never by inspecting
    op.argument_symmetry directly.
    """
    TYPE_11  = "11"
    TYPE_NEG = "NEG"
    TYPE_12  = "12"
    TYPE_21  = "21"
    TYPE_22  = "22"


def _sign_correlation_type_from_pf(pf) -> SignCorrelationType:
    """
    Compute the SignCorrelationType for a PairFlavour from its structural
    parameters alone (no specific atom instances required).

    Algorithm: per-group achievable-set product.
    """
    antisym_u = pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
    antisym_v = pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

    achievable: set[tuple[int, int]] = {(1, 1)}

    for ku_g, kv_g, s_g in zip(
        pf.flavour_u.counts, pf.flavour_v.counts, pf.overlap
    ):
        pu_can_flip = antisym_u and ku_g >= 2
        pv_can_flip = antisym_v and kv_g >= 2
        u_g_eq_v_g  = (s_g == ku_g == kv_g)

        if not pu_can_flip and not pv_can_flip:
            g_ach: set[tuple[int, int]] = {(1, 1)}
        elif pu_can_flip and not pv_can_flip:
            g_ach = {(1, 1), (-1, 1)}
        elif not pu_can_flip and pv_can_flip:
            g_ach = {(1, 1), (1, -1)}
        elif u_g_eq_v_g:
            g_ach = {(1, 1), (-1, -1)}
        else:
            g_ach = {(1, 1), (1, -1), (-1, 1), (-1, -1)}

        achievable = {
            (pu * qu, pv * qv)
            for (pu, pv) in achievable
            for (qu, qv) in g_ach
        }

    u_flips   = (-1,  1) in achievable
    v_flips   = ( 1, -1) in achievable
    both_flip = (-1, -1) in achievable

    if u_flips and v_flips:
        return SignCorrelationType.TYPE_22
    elif u_flips:
        return SignCorrelationType.TYPE_21
    elif v_flips:
        return SignCorrelationType.TYPE_12
    elif both_flip:
        return SignCorrelationType.TYPE_NEG
    else:
        return SignCorrelationType.TYPE_11


def _complex_to_reals(c: np.ndarray) -> np.ndarray:
    """Unpack n complex values to 2n reals: [re0, im0, re1, im1, ...]."""
    out = np.empty(2 * len(c), dtype=float)
    out[0::2] = c.real
    out[1::2] = c.imag
    return out
