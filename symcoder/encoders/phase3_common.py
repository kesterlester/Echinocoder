"""
symcoder.encoders.phase3_common
================================
Shared helpers for the Phase 3 one-shot whole-table encoders
(``phase3_simplicial.py``, ``phase3_vandermonde.py``).

A "Phase 3 encoder" takes a plan + event and produces a permutation-invariant
embedding of the *full* alignment table

    T(event) ∈ R^{|G| x |rep|},     T[g, i] = eval(g · rep_atoms[i], event)

viewed as a multiset of |G| tuples in R^{|rep|}.  Different concrete Phase 3
encoders (simplicial, Vandermonde C∞, Vandermonde C⁰) all consume the same
``T``; they differ only in how the multiset is rendered as reals.

This module provides:

* ``build_alignment_table(plan, event, the_group, rep_atoms)`` — the table
  computation, in one place rather than duplicated across encoders.
* ``rep_atoms_for(plan)`` — the canonical list of repS atoms used by every
  Phase 3 encoder, so they all agree on column ordering.

Why the rep atom list lives here
--------------------------------
The pruning from ``rep`` → ``repS`` is unsafe for whole-table encoders (see
``DOCS/parity_blindness_and_repS_concern.pdf``).  Every Phase 3 encoder must
therefore use ``rep`` (the unpruned set).  Centralising the choice here makes
it impossible for a future Phase 3 implementation to accidentally reach for
``repS`` instead.
"""
from __future__ import annotations

import numpy as np

from symatom import rep as _rep_fn
from symcoder.eval import evaluate


def rep_atoms_for(plan) -> list:
    """Return the unpruned list of atoms over which Phase 3 encoders work.

    Always returns ``rep`` (not ``repS``); see module docstring.
    """
    return list(_rep_fn(plan.context, plan.operations))


def build_alignment_table(plan, event: dict, the_group, rep_atoms) -> np.ndarray:
    """Return T of shape (|G|, |rep|), dtype float64.

    T[i, j] = eval(g_i · rep_atoms[j], event)

    The simplicial and Vandermonde encoders both consume this same table; only
    the way they render it as a vector of reals differs.
    """
    group_elements = list(the_group.all_group_elements())
    n = len(group_elements)
    k = len(rep_atoms)
    T = np.empty((n, k), dtype=np.float64)
    for i, g in enumerate(group_elements):
        for j, atom in enumerate(rep_atoms):
            T[i, j] = float(evaluate(g.apply(atom), event))
    return T


def column_mass_dimensions(rep_atoms) -> list:
    """Return ``[atom.operation.mass_dimension for atom in rep_atoms]``.

    Raises ``ValueError`` if any operation lacks a declared mass_dimension.
    Used by encoders that wish to pre-divide each column of the alignment
    table by a user-supplied length scale raised to that operation's
    homogeneity degree.
    """
    out = []
    for j, atom in enumerate(rep_atoms):
        md = getattr(atom.operation, "mass_dimension", None)
        if md is None:
            raise ValueError(
                f"column_mass_dimensions: rep atom {j} (operation "
                f"{atom.operation.name!r}) has no declared mass_dimension. "
                f"Encoders that wish to use scale-based pre-normalisation "
                f"require every operation in the plan to declare a "
                f"mass_dimension on its Operation constructor.  "
                f"E.g. ``Operation(..., mass_dimension=2)`` for a rank-2 "
                f"symmetric dot product."
            )
        out.append(int(md))
    return out


def apply_scale(T: np.ndarray, rep_atoms, scale) -> np.ndarray:
    """Divide each column of T by ``scale ** rep_atoms[j].operation.mass_dimension``.

    The result is a copy with the same shape and dtype as T.  When ``scale``
    is ``None`` the input is returned unchanged (no copy).

    The intent is to bring rows of T to a common dimensional footing before
    they are mixed in a non-linear way (e.g. raised to high powers in the
    Cinf Vandermonde encoder).  Without rescaling, a row carrying values
    with units of length-cubed and another carrying values with units of
    length-squared will produce wildly different magnitudes once raised to
    the n-th power, making the float64 output numerically unbalanced.

    Mathematical effect
    -------------------
    If the original event vectors were in some unit (say ``MeV``), then
    ``mag(v)`` has units ``MeV``, ``dot(u, v)`` has units ``MeV^2``, and
    ``eps3(u, v, w)`` has units ``MeV^3``.  Dividing the ``mag`` column by
    ``scale``, the ``dot`` column by ``scale^2``, and the ``eps3`` column
    by ``scale^3`` produces a *dimensionless* table whose entries are
    all of order unity when ``scale`` matches the typical event scale.

    Parameters
    ----------
    T : np.ndarray, shape (n, k)
        Alignment table as produced by ``build_alignment_table``.
    rep_atoms : sequence of Atom
        One entry per column of T; each must have
        ``atom.operation.mass_dimension`` set.
    scale : float | None
        If None, return T unchanged.  Otherwise must be a positive real;
        the j-th column is divided by ``scale ** mass_dim[j]``.

    Returns
    -------
    np.ndarray
        Scaled copy of T.  Same shape and dtype.
    """
    if scale is None:
        return T
    if scale <= 0:
        raise ValueError(
            f"apply_scale: scale must be a positive real, got {scale!r}.  "
            f"(A length scale cannot be derived from the event itself "
            f"because the event values can be arbitrarily close to zero; "
            f"the scale must be supplied externally as an encoder parameter.)"
        )
    mds = column_mass_dimensions(rep_atoms)
    # Build the divisor row in float64 once, then broadcast across rows.
    divisor = np.array([float(scale) ** md for md in mds], dtype=np.float64)
    return (T / divisor[None, :]).astype(np.float64, copy=False)
