"""
symcoder.encoders.phase3_simplicial
====================================
Phase 3 encoder: applies the C0-Hom-degree-1 simplicial-complex multiset
embedder (lives at the Echinocoder repo root) to the full alignment table
T of shape (|G|, |rep|) — i.e. the multiset of repS-evaluated column
vectors over the discrete group orbit at the event.

Motivation
----------
Phase 1 (per-row multisets) and Phase 2 (per-row-pair joint multisets)
are not in general sufficient to distinguish two events that lie in
distinct SO(3) × S_n orbits but which happen to have palindromic
chirality structure (see symcoder/tests/test_parity_blindness_shadows.py).

Phase 3 closes that gap by embedding the *entire* alignment table as a
single multiset of |G| tuples in R^|repS|.  The simplicial-complex
embedder is deterministically faithful for any multiset of (n, k) reals
and outputs ``2*n*k + 1 - k`` reals (its ``size_from_n_k_generic``).

Faithfulness vs smoothness
--------------------------
The simplicial embedder is piecewise linear (C^0 but not C^1).  It is
intentionally appended to the end of the encoding, with a Phase3Tree
descriptor exposing its start offset, so downstream consumers that care
only about smooth structure can slice it off and use Phase 1 + Phase 2
alone.

Architecture
------------
Mirrors the OrbitEncoder / Phase2Encoder pattern:

* ``Phase3SimplicialEncoder``         — built object, .encode(event) → EncodingResult,
                                        .describe(start_offset) → Phase3Tree.
* ``Phase3SimplicialEncoderFactory``  — .build(plan) → Phase3SimplicialEncoder.

Use via ``symcoder.encode`` and ``symcoder.encode_and_describe``; pass a
factory instance for ``phase3_factory`` (or ``None`` to skip Phase 3).
The top-level entry points default to a Phase3SimplicialEncoderFactory()
so that Phase 3 is on by default.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np

from symatom import rep as _rep_fn
from symcoder.encoders._base import EncodingResult
from symcoder.eval import evaluate

# ---------------------------------------------------------------------------
# Lazy-import the top-level simplicial embedder (it lives outside any package).
# ---------------------------------------------------------------------------
_REPO_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# Two interchangeable implementations exist at the repo root; embedder 2 has
# a tighter output dim and is otherwise indistinguishable from embedder 1.
from C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset import (  # type: ignore
    Embedder as _SimplicialEmbedder,
)


# ---------------------------------------------------------------------------
# Phase3Tree (lightweight descriptor; parallels Phase1Tree / Phase2Tree).
# Imported lazily to avoid circular dependency with symcoder.describe at
# module-load time.
# ---------------------------------------------------------------------------

def _Phase3Tree_class():
    from symcoder.describe import Phase3Tree
    return Phase3Tree

def _SegmentInfo_class():
    from symcoder.describe import SegmentInfo
    return SegmentInfo


# ---------------------------------------------------------------------------
# Phase 3 encoder
# ---------------------------------------------------------------------------

class Phase3SimplicialEncoder:
    """
    Encode the full alignment table T(event) of shape (|G|, |rep|) as a
    permutation-invariant multiset via the simplicial-complex embedder.

    The alignment table is built canonically: rows indexed by group elements
    in TheGroup.all_group_elements() order, columns indexed by the rep atoms
    in canonical order.  T_{g, i} = eval(g · atom_i, event).  Although the
    row index ordering is canonical here, the simplicial embedder is
    multiset-invariant over rows, so the encoding is in fact independent of
    any permutation of group-element ordering — exactly what label-symmetry
    invariance requires.

    Why ``rep`` and not ``repS``
    ----------------------------
    ``repS`` is a pre-pruned version of ``rep`` whose pruning was justified
    only under the (now-known-false) assumption that Phase 1 + Phase 2
    marginals are sufficient for faithful encoding.  Using ``repS`` here
    would repeat the same mistake: rows that look ``duplicative'' at the
    pair-marginal level may carry distinct information once we look at the
    full alignment table.  See
    ``DOCS/parity_blindness_and_repS_concern.pdf`` for the full history.
    """

    def __init__(self, plan) -> None:
        self._plan       = plan
        self._the_group  = plan.context.the_group
        self._rep_atoms  = _rep_fn(plan.context, plan.operations)
        self._embedder   = _SimplicialEmbedder()

        n = self._the_group.order()
        k = len(self._rep_atoms)
        self._n = n
        self._k = k
        self._output_dim = int(self._embedder.size_from_n_k(n, k))

    @property
    def output_dim(self) -> int:
        return self._output_dim

    def _build_table(self, event: dict) -> np.ndarray:
        """Return T of shape (|G|, |rep|), dtype float64."""
        group_elements = list(self._the_group.all_group_elements())
        n = len(group_elements)
        k = len(self._rep_atoms)
        T = np.empty((n, k), dtype=np.float64)
        for i, g in enumerate(group_elements):
            for j, atom in enumerate(self._rep_atoms):
                T[i, j] = float(evaluate(g.apply(atom), event))
        return T

    def encode(self, event: dict) -> EncodingResult:
        T = self._build_table(event)
        embedding, _shape, _meta = self._embedder.embed(T)
        values = np.asarray(embedding, dtype=np.float64)
        # Defensive: size must match algebraic prediction.
        assert values.shape == (self._output_dim,), (
            f"Phase 3 simplicial embedding produced {values.shape[0]} reals; "
            f"expected {self._output_dim} (=size_from_n_k({self._n},{self._k}))."
        )
        return EncodingResult(values=values)

    def describe(self, start_offset: int = 0):
        """Return a Phase3Tree describing this Phase-3 block.

        When ``output_dim`` is zero (degenerate plan — no operations or an
        empty group) we emit an empty Phase3Tree so that the overall
        EncodingTree is empty for trivial plans, preserving the original
        invariant that ``describe_encoding`` returns an empty list for an
        empty plan.
        """
        SegmentInfo = _SegmentInfo_class()
        Phase3Tree  = _Phase3Tree_class()
        if self._output_dim == 0:
            return Phase3Tree(segments=[], n=self._n, k=self._k)
        seg = SegmentInfo(
            kind         = "SIMPLICIAL",
            start        = start_offset,
            length       = self._output_dim,
            op_u         = None,          # spans all of repS — not per-operation
            flavour_u    = None,
            op_v         = None,
            flavour_v    = None,
            overlap      = None,
            symmetry_class  = None,
            sign_compressed = None,
            method_name  = "simplicial_c0_deg1",
            example      = (f"simplicial embedding of full rep orbit table "
                            f"(n=|G|={self._n}, k=|rep|={self._k})"),
        )
        return Phase3Tree(segments=[seg], n=self._n, k=self._k)


# ---------------------------------------------------------------------------
# Phase 3 factory
# ---------------------------------------------------------------------------

class Phase3SimplicialEncoderFactory:
    """
    Factory for Phase3SimplicialEncoder.  Sole entry point for callers that
    want to opt in/out of Phase 3 at the top-level ``encode`` API.

    Passing ``phase3_factory=None`` to ``encode`` / ``encode_and_describe``
    suppresses Phase 3 entirely; passing an instance of this factory turns
    it on.  Defaults at the call sites are ``Phase3SimplicialEncoderFactory()``
    — i.e. Phase 3 is on by default.
    """

    def build(self, plan) -> Phase3SimplicialEncoder:
        return Phase3SimplicialEncoder(plan)
