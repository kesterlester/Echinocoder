"""
cross_validate.py
-----------------
Run this script before and after a refactoring to verify that no numerical
output has changed.

Usage
-----
Before refactoring — save a baseline:
    python scripts/cross_validate.py --save

After refactoring — check outputs match:
    python scripts/cross_validate.py --check

The baseline is stored in scripts/cross_validate_baseline.npy.

The test cases
--------------
  - Context: 3 electrons (A,B,C), 2 muons (P,Q), 5D vectors.
  - Operations: dot (SYMMETRIC rank-2) and eps (ANTISYMMETRIC rank-2).
  - Events: 5 deterministic pseudo-random events (rng seed 42).

If the refactoring changes the NAMES of things but not the algorithm, the
numerical output must be bit-for-bit identical.  Any difference is a bug.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# ── locate repo root so we can import the packages regardless of CWD ─────────
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO))

# ── imports — update these names if the API was renamed ──────────────────────
# Rename VectorType → VectorType here once step 1 is committed.
from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan   # noqa: E402
from symcoder import EvaluableOperation, encode                                # noqa: E402

_BASELINE = Path(__file__).parent / "cross_validate_baseline.npy"
_VEC_DIM  = 5       # dimensionality of each vector
_N_EVENTS = 5
_SEED     = 42


# ── fixed physics setup ───────────────────────────────────────────────────────

def _make_plan() -> Plan:
    electrons = VectorType("electrons", ("A", "B", "C"))
    muons     = VectorType("muons",     ("P", "Q"))
    ctx = Context((electrons, muons))

    dot = EvaluableOperation(
        "dot", rank=2, parity=+1,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda vecs: float(np.dot(vecs[0], vecs[1])),
    )
    eps = EvaluableOperation(
        "eps", rank=2, parity=-1,
        argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
        # 2D antisymmetric invariant embedded in higher-dim vectors
        eval_fn=lambda vecs: float(vecs[0][0] * vecs[1][1] - vecs[0][1] * vecs[1][0]),
    )
    return Plan(context=ctx, operations=(dot, eps))


def _make_events(plan: Plan) -> list[dict]:
    rng    = np.random.default_rng(_SEED)
    labels = plan.context.all_labels
    return [
        {lbl: rng.standard_normal(_VEC_DIM) for lbl in labels}
        for _ in range(_N_EVENTS)
    ]


def _run() -> np.ndarray:
    plan   = _make_plan()
    events = _make_events(plan)
    return np.stack([encode(plan, ev) for ev in events])


# ── CLI ───────────────────────────────────────────────────────────────────────

def save_baseline() -> None:
    data = _run()
    np.save(str(_BASELINE), data)
    print(f"Baseline saved to {_BASELINE}")
    print(f"  shape: {data.shape}  dtype: {data.dtype}")
    print(f"  first row (first 8 values): {data[0, :8]}")


def check_against_baseline() -> None:
    if not _BASELINE.exists():
        print(f"ERROR: baseline file not found: {_BASELINE}")
        print("Run with --save first.")
        sys.exit(1)

    baseline = np.load(str(_BASELINE))
    current  = _run()

    shape_ok   = baseline.shape == current.shape
    values_ok  = np.allclose(baseline, current, rtol=0, atol=0)   # bit-exact
    print(f"Baseline shape : {baseline.shape}")
    print(f"Current  shape : {current.shape}  {'OK' if shape_ok else 'MISMATCH'}")
    if not shape_ok:
        sys.exit(1)
    print(f"Max abs diff   : {np.max(np.abs(baseline - current))}")
    print(f"Values match   : {'YES (bit-exact)' if values_ok else 'NO — MISMATCH'}")
    if not values_ok:
        sys.exit(1)
    print("Cross-validation PASSED.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cross-validate encoding outputs.")
    group  = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--save",  action="store_true", help="Save baseline.")
    group.add_argument("--check", action="store_true", help="Check against baseline.")
    args = parser.parse_args()

    if args.save:
        save_baseline()
    else:
        check_against_baseline()
