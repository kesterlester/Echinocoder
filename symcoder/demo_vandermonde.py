#!/usr/bin/env python3
"""
symcoder/demo_vandermonde.py
============================
Standalone Phase-3 Vandermonde encoder demonstration.

Runs only the Phase-3 Vandermonde encoder
(``Phase3VandermondeEncoderFactory``) on a chosen event.  No Phase 1, no
Phase 2, no alignment decoder — the Vandermonde encoder is faithful on
its own.  This script exists to show what the Vandermonde encoder sees,
computes, and outputs in isolation.

Switch ``USE_C0`` below to toggle between the two modes:

* ``USE_C0 = False`` (default): the C-infinity DFT-node power-sum
  variant.  Smooth everywhere, internally complex, output dimensionally
  inhomogeneous (units span ``input^1`` … ``input^n``).
* ``USE_C0 = True``: the C-zero sort variant.  Continuous (Lipschitz)
  but not differentiable at value-ties, internally real, output shares
  units with the input.

Outputs simultaneously to:
  stdout                — plain text
  demo_vandermonde.tex  — LaTeX source

Run from repo root:
    venv/bin/python symcoder/demo_vandermonde.py

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
from symcoder.encoders.phase3_vandermonde import Phase3VandermondeEncoderFactory
from symcoder.tests._parity_shadow_events import SHADOW_EVENTS


# ─── User-configurable switches ─────────────────────────────────────────────
USE_C0       = False  # True → C0 sort-mode; False → Cinf DFT-power-sum mode.
SCALE        = 4      # Optional positive real to pre-normalise alignment-table
                      # columns by the operation's mass_dimension power.  Only
                      # affects the Cinf mode in any visible way.  None = off.
SHADOW_INDEX = 0      # which shadow event to demonstrate on.
# ────────────────────────────────────────────────────────────────────────────


_HERE     = os.path.dirname(os.path.abspath(__file__))
TEX_FILE  = os.path.join(_HERE, "demo_vandermonde.tex")


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


def _intro_cb(out, mode: str, scale):
    scale_str = "none" if scale is None else f"{scale:g}"
    out.section("About this document")
    out.line(
        f"This document demonstrates the Phase-3 Vandermonde one-shot "
        f"encoder running on its own, with Phase 1 and Phase 2 disabled.  "
        f"The Vandermonde encoder is deterministically faithful on the "
        f"full rep alignment table, so it is a valid standalone encoder; "
        f"this document shows exactly what it consumes (the rep atoms and "
        f"the alignment table T) and what it emits.",
        tex=(
            r"This document demonstrates the Phase-3 Vandermonde one-shot "
            r"encoder running on its own, with Phase 1 and Phase 2 "
            r"disabled.  The Vandermonde encoder is deterministically "
            r"faithful on the full rep alignment table, so it is a valid "
            r"standalone encoder; this document shows exactly what it "
            r"consumes (the rep atoms and the alignment table $T$) and "
            r"what it emits.\\[1ex]"
        ),
    )
    # Note: in stdout we use the literal lambda character; in TeX we use the
    # math-mode macro \lambda since LaTeX's default font cannot render U+03BB.
    out.kv("Mode", mode,
           tex=rf"\textbf{{Mode:}} \texttt{{{mode}}}\\")
    out.kv("Scale (lambda)", scale_str,
           tex=rf"\textbf{{Scale ($\lambda$):}} \texttt{{{scale_str}}}\\")
    out.line(
        "See DOCS/phase3_vandermonde.pdf for the mathematical derivation "
        "(DFT-node power sums, conjugate-pair compression, optional "
        "dimensional rescaling).",
        tex=(r"See \texttt{DOCS/phase3\_vandermonde.pdf} for the "
             r"mathematical derivation (DFT-node power sums, conjugate-pair "
             r"compression, optional dimensional rescaling).\\[1ex]"),
    )


def main():
    mode = "C0" if USE_C0 else "Cinf"
    raw       = SHADOW_EVENTS[SHADOW_INDEX]
    event     = {lbl: np.array(v, dtype=float) for lbl, v in raw.items()}
    plan, ctx = _build_plan()

    print(f"\n### Phase-3 Vandermonde demo, mode={mode}, scale={SCALE}, "
          f"shadow event #{SHADOW_INDEX} ###")

    run_demo(
        plan, ctx, event, TEX_FILE,
        title_main=f"symcoder Phase-3 Vandermonde demo ({mode})",
        subtitle=(f"standalone Vandermonde encoder ({mode}, "
                  f"scale={SCALE!r}), shadow event #{SHADOW_INDEX}"),
        intro_callback=lambda out: _intro_cb(out, mode=mode, scale=SCALE),
        # Phase 1 and Phase 2 disabled; Phase 3 set to Vandermonde only.
        orbit_factory=None,
        phase2_factory=None,
        phase3_factory=Phase3VandermondeEncoderFactory(mode=mode, scale=SCALE),
    )


if __name__ == "__main__":
    main()
