#!/usr/bin/env python3
"""
symcoder/demo_shadow_event.py
=============================
Parity-blindness demonstration.

This script runs the full encode → decode pipeline **twice** with the same
plan, once on a chosen "shadow event" E and once on its parity flip
P(E) = {label: -v for label, v in E.items()}:

  1.  E      →  demo_shadow_event_parity_+1.tex
  2.  P(E)   →  demo_shadow_event_parity_-1.tex

Side-by-side, the two PDFs make the parity-blindness of the current
Echinocoder Phase 1 + Phase 2 encoder visible: the encoded numerical values
are identical in the two documents even though E and P(E) lie in distinct
SO(3) orbits (their signed scalar-triple products eps3(...) carry opposite
sign).  See ``symcoder/tests/test_parity_blindness_shadows.py`` for a
formal pytest-driven regression that documents this property.

The plan used here is deliberately limited to ``Mymag`` + ``euclidean3.dot``
+ ``euclidean3.eps`` (no ``euclidean2.dot``), because the parity-blindness
shadows were discovered against exactly this plan and the matching
SO(3) × S_n equivalence oracle.

Run from repo root:
    venv/bin/python symcoder/demo_shadow_event.py

This script is a *thin driver*; the heavy lifting (LaTeX preamble, Phase 1
and Phase 2 narration, alignment-decoder section, summary) lives in
``symcoder/_demo_helper.py`` and is shared with ``demo_roundtrip.py``.
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
from symcoder.tests._parity_shadow_events import SHADOW_EVENTS


_HERE          = os.path.dirname(os.path.abspath(__file__))
TEX_FILE_POS   = os.path.join(_HERE, "demo_shadow_event_parity_+1.tex")
TEX_FILE_NEG   = os.path.join(_HERE, "demo_shadow_event_parity_-1.tex")
SHADOW_INDEX   = 0   # which entry of SHADOW_EVENTS to demonstrate


def _build_plan():
    """Pure 3D plan: Mymag + euclidean3.dot + euclidean3.eps."""
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


def _make_intro_callback(parity_sign: int, shadow_index: int):
    """Return a callable(out) that emits the 'About this document' banner.

    The banner explains what the document is for and identifies the
    companion document so reviewers can find both halves.
    """
    parity_str = "+1" if parity_sign > 0 else "-1"
    other_str  = "-1" if parity_sign > 0 else "+1"
    companion  = f"demo_shadow_event_parity_{other_str}.tex"

    def cb(out):
        out.section("About this document")
        out.line(
            "This document is one half of a parity-blindness demonstration.  "
            "A companion document, identical in structure but built from the "
            "parity-flipped event, exists alongside it.  Their Phase 1 and "
            "Phase 2 encoded values should match exactly, despite the events "
            "lying in distinct SO(3) orbits.",
            tex=(
                r"This document is one half of a parity-blindness demonstration.  "
                r"A companion document, identical in structure but built from the "
                r"parity-flipped event, exists alongside it.  Their Phase 1 and "
                r"Phase 2 encoded values should match exactly, despite the events "
                r"lying in distinct $SO(3)$ orbits."
                r"\\[1ex]"
            ),
        )
        out.kv(
            "This run uses event with parity sign", parity_str,
            tex=(rf"\textbf{{This run uses event with parity sign:}} "
                 rf"${parity_str}$\\")
        )
        out.kv(
            "Shadow event index", str(shadow_index),
            tex=(rf"\textbf{{Shadow event index:}} ${shadow_index}$ "
                 rf"(from \texttt{{symcoder.tests.\_parity\_shadow\_events.SHADOW\_EVENTS}})\\")
        )
        companion_tex = companion.replace("_", r"\_")
        out.kv(
            "Companion document", companion,
            tex=(rf"\textbf{{Companion document:}} "
                 rf"\texttt{{{companion_tex}}}\\")
        )
    return cb


def main():
    raw       = SHADOW_EVENTS[SHADOW_INDEX]
    event_pos = {lbl: np.array(v, dtype=float) for lbl, v in raw.items()}
    event_neg = {lbl: -v for lbl, v in event_pos.items()}
    plan, ctx = _build_plan()

    runs = [
        (+1, event_pos, TEX_FILE_POS),
        (-1, event_neg, TEX_FILE_NEG),
    ]
    for sign, event, tex_path in runs:
        parity_str = "+1" if sign > 0 else "-1"
        print(f"\n### shadow event #{SHADOW_INDEX}, parity sign = {parity_str} ###")
        run_demo(
            plan, ctx, event, tex_path,
            title_main="symcoder parity-blindness demo",
            subtitle=f"shadow event #{SHADOW_INDEX},  parity = {parity_str}",
            intro_callback=_make_intro_callback(sign, SHADOW_INDEX),
            use_comp_drop=False,
        )


if __name__ == "__main__":
    main()
