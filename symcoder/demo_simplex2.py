#!/usr/bin/env python3
"""
symcoder/demo_simplex2.py
=========================
Standalone Phase-3 simplicial encoder demonstration.

Runs only the Phase-3 simplicial-complex encoder
(``Phase3SimplicialEncoderFactory``) on a chosen event.  No Phase 1, no
Phase 2, no alignment decoder — the simplicial encoder is faithful on its
own.  This script exists to show what the Phase-3 simplicial encoder
sees, computes, and outputs in isolation.

Outputs simultaneously to:
  stdout              — plain text
  demo_simplex2.tex   — LaTeX source

Run from repo root:
    venv/bin/python symcoder/demo_simplex2.py

This script is a *thin driver*; the LaTeX preamble, the rep atom and
alignment-table narration, and the Phase-3 segment descriptor display all
live in ``symcoder/_demo_helper.py`` and are shared with the other demos.
"""
from __future__ import annotations

import os
import sys

import numpy as np

# ── locate repo root and activate imports ──────────────────────────────────
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan
import symcoder.operations.euclidean3
from symcoder._demo_helper import run_demo
from symcoder.encoders.phase3_simplicial import Phase3SimplicialEncoderFactory
from symcoder.tests._parity_shadow_events import SHADOW_EVENTS


_HERE          = os.path.dirname(os.path.abspath(__file__))
TEX_FILE       = os.path.join(_HERE, "demo_simplex2.tex")
SHADOW_INDEX   = 0   # use shadow event #0 so the demo also illustrates the
                     # parity-blindness fix (Phase 3 distinguishes E from P(E)
                     # on this event, where Phase 1 + Phase 2 cannot).


def _build_plan():
    """3-electron / 2-muon plan with Mymag + dot3 + eps3 (no euclidean2.dot)."""
    mymag = Operation(
        "Mymag", rank=1, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))),
        tex=r"|\vc{#1}|",
        mass_dimension=1,
    )
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    plan      = Plan(context=ctx, operations=(
        mymag,
        symcoder.operations.euclidean3.dot,
        symcoder.operations.euclidean3.eps,
    ))
    return plan, ctx


def _intro_cb(out):
    out.section("About this document")
    out.line(
        "This document demonstrates the Phase-3 simplicial-complex one-shot "
        "encoder running on its own, with Phase 1 and Phase 2 disabled.  "
        "The simplicial encoder is deterministically faithful on the full "
        "rep alignment table, so it is a valid standalone encoder; this "
        "document shows exactly what it consumes (the rep atoms and the "
        "alignment table T) and what it emits.",
        tex=(
            r"This document demonstrates the Phase-3 simplicial-complex "
            r"one-shot encoder running on its own, with Phase 1 and Phase 2 "
            r"disabled.  The simplicial encoder is deterministically "
            r"faithful on the full rep alignment table, so it is a valid "
            r"standalone encoder; this document shows exactly what it "
            r"consumes (the rep atoms and the alignment table $T$) and what "
            r"it emits.\\[1ex]"
        ),
    )


def main():
    raw       = SHADOW_EVENTS[SHADOW_INDEX]
    event     = {lbl: np.array(v, dtype=float) for lbl, v in raw.items()}
    plan, ctx = _build_plan()

    print(f"\n### Phase-3 simplicial demo, shadow event #{SHADOW_INDEX} ###")
    run_demo(
        plan, ctx, event, TEX_FILE,
        title_main="symcoder Phase-3 simplicial demo",
        subtitle=f"standalone simplicial encoder, shadow event #{SHADOW_INDEX}",
        intro_callback=_intro_cb,
        # Phase 1 and Phase 2 disabled; Phase 3 set to simplicial only.
        orbit_factory=None,
        phase2_factory=None,
        phase3_factory=Phase3SimplicialEncoderFactory(),
    )


if __name__ == "__main__":
    main()
