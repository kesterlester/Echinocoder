"""
symcoder.tests.test_parity_blindness_shadows
=============================================
Regression test for parity-blindness in the current Echinocoder Phase 1 +
Phase 2 encoder.

Background
----------
For an event E (a labelled tuple of vectors in R^3), define its parity flip
as P(E) = {label: -v for label, v in E.items()}.  P is a determinant-(-1)
orthogonal transformation; events E and P(E) lie in *different* SO(3)
orbits of R^3, distinguished by the sign of every signed triple product
(every value of the eps3 operation).

Because the plan in this test uses ``euclidean3.eps`` (a pseudoscalar), the
encoder is *intended* to be SO(3)-faithful: encode(E) and encode(P(E))
should differ whenever E is not itself a fixed point of parity.

It turns out (empirically — see ``symcoder.experiments.encoding_fidelity_search``)
that for certain structured events the encoder silently collapses E and P(E)
into the same output.  The 31 events in ``_parity_shadow_events.SHADOW_EVENTS``
are such "shadow" events: each one satisfies

    encode(E) == encode(P(E))   exactly, even though Gram + signed-det
                                 invariants distinguish them.

The root cause is that Phase 1 stores per-row *multisets* of operation
values, and the multiset of any antisymmetric operation's orbit values is
necessarily palindromic under negation (every value comes with its negation
elsewhere in the orbit).  Phase 2's joint multisets break this symmetry on
generic inputs but, on these structured inputs, happen to be palindromic
too.

Status
------
Each parametrised case is marked ``xfail(strict=True)``: it is expected to
fail today.  When the encoder is extended (e.g. with a Phase 3 step that
captures chirality, or with a tweak to Phase 2 that preserves orientation
on degenerate inputs), all 31 cases should pass and the xfail marker should
be removed.

To regenerate the fixture (e.g. after adding more event families to the
search), use the recipe in ``_parity_shadow_events.py``'s module docstring.
"""
from __future__ import annotations

import numpy as np
import pytest

from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan
import symcoder.operations.euclidean3
from symcoder.encoders import (
    OrbitEncoderFactory, SortEncoderFactory, HalfSortEncoderFactory,
    standard_row_pair_factories, OverlapBlockEncoderFactory, Phase2EncoderFactory,
)
from symcoder.tests._parity_shadow_events import SHADOW_EVENTS


ATOL = 1e-9
"""Absolute tolerance below which two encoding values are considered equal.

The encoder is deterministic and uses sorted multisets, so identical inputs
produce *exactly* identical outputs in IEEE 754.  Differences at or below
ATOL therefore indicate "encoder cannot distinguish", not floating-point
noise.  ATOL is set 9 decimal places to leave plenty of headroom against
spurious noise without masking genuine sub-1e-6 differences.
"""


def _make_plan():
    """Pure 3D plan with eps3, matching the search that produced the shadows."""
    mymag = Operation(
        "Mymag", rank=1, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))),
    )
    electrons = VectorType("electrons", ("a", "b", "c"))
    muons     = VectorType("muons",     ("p", "q"))
    ctx       = Context((electrons, muons))
    plan      = Plan(context=ctx, operations=(
        mymag,
        symcoder.operations.euclidean3.dot,
        symcoder.operations.euclidean3.eps,
    ))
    return plan


def _build_encoders(plan):
    orbit_fac  = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
    phase2_fac = Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories(),
                                   use_complementarity_drop=False)
    ])
    return orbit_fac.build(plan), phase2_fac.build(plan)


def _encode(orbit_enc, phase2_enc, event):
    ph1 = np.asarray(orbit_enc.encode(event).values,  dtype=float)
    ph2 = np.asarray(phase2_enc.encode(event).values, dtype=float)
    return np.concatenate([ph1, ph2])


def _to_event(event_dict):
    return {k: np.array(v, dtype=float) for k, v in event_dict.items()}


def _parity_flip(event):
    return {k: -v for k, v in event.items()}


@pytest.fixture(scope="module")
def encoders():
    plan = _make_plan()
    return _build_encoders(plan)


@pytest.mark.parametrize("shadow_idx,raw_event",
                         list(enumerate(SHADOW_EVENTS)),
                         ids=[f"shadow_{i:02d}" for i in range(len(SHADOW_EVENTS))])
@pytest.mark.xfail(
    strict=True,
    reason=(
        "Known parity-blindness of the current Phase 1 + Phase 2 encoder on "
        "structured events.  encode(E) == encode(P(E)) for this shadow event "
        "even though E and its parity flip are in different SO(3) orbits.  "
        "See module docstring."
    ),
)
def test_parity_flip_distinguishable(encoders, shadow_idx, raw_event):
    """encode(E) and encode(P(E)) must differ by more than ATOL.

    Currently fails (xfail) for every shadow event: the encoder produces
    *exactly* the same output for E and P(E).
    """
    orbit_enc, phase2_enc = encoders
    event = _to_event(raw_event)
    flipped = _parity_flip(event)

    enc_pos = _encode(orbit_enc, phase2_enc, event)
    enc_neg = _encode(orbit_enc, phase2_enc, flipped)

    max_diff = float(np.max(np.abs(enc_pos - enc_neg)))
    assert max_diff > ATOL, (
        f"shadow #{shadow_idx}: encoder is parity-blind on this event.\n"
        f"  encode(E)  - encode(P(E)) max |diff| = {max_diff:.3e}\n"
        f"  required:                            >  {ATOL:.3e}\n"
        f"  E = {raw_event}"
    )
