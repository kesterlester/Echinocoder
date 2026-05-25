#!/usr/bin/env python3
"""
symcoder/demo_roundtrip.py
===========================
Round-trip encode → decode demonstration for a small plan.

Context   : electrons = (a, b, c),  muons = (p, q)
Operations: Mymag (user-defined), euclidean3.dot, euclidean2.dot, euclidean3.eps
Group     : G = S_3 × S_2  (order 12)

Outputs simultaneously to:
  stdout               — plain text
  demo_roundtrip.tex   — LaTeX source; compile with  pdflatex demo_roundtrip.tex

Run from repo root:
  venv/bin/python symcoder/demo_roundtrip.py

This script is a *thin driver*.  Everything from the LaTeX preamble down to
the alignment-decoder narration lives in ``symcoder/_demo_helper.py`` and is
shared with ``demo_shadow_event.py``; this file is just a place to choose a
plan and a fixed event.
"""
from __future__ import annotations

import os
import sys

import numpy as np

# ── locate repo root and activate imports ──────────────────────────────────
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan
import symcoder.operations.euclidean2
import symcoder.operations.euclidean3
from symcoder._demo_helper import run_demo


_TEX_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "demo_roundtrip.tex")


def main():
    # Mymag is defined here to show that user-defined operations work
    # just as well as library ones.  The others come from the standard library.
    mymag = Operation(
        "Mymag", rank=1, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))),
        tex=r"|\vc{#1}|",
    )

    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    plan      = Plan(context=ctx, operations=(
        mymag,
        symcoder.operations.euclidean3.dot,
        symcoder.operations.euclidean2.dot,
        symcoder.operations.euclidean3.eps,
    ))

    # Fixed event (3-D vectors; eps3 needs ≥3 dimensions).
    event = {
        "a": np.array([ 1.2,  0.3,  0.5]),
        "b": np.array([ 0.4,  1.1,  0.2]),
        "c": np.array([-0.5,  1.2, -0.7]),
        "p": np.array([ 0.0,  0.0, +1.0]),
        "q": np.array([ 0.0,  0.0, -1.0]),
    }

    subtitle = ",  ".join(
        f"{t.name}=({','.join(t.labels)})" for t in ctx.types
    )

    run_demo(
        plan, ctx, event, _TEX_FILE,
        title_main="symcoder round-trip demo",
        subtitle=subtitle,
        use_comp_drop=False,
    )


if __name__ == "__main__":
    main()
