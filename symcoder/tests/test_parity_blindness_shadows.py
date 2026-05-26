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
from symcoder import encode as symcoder_encode
from symcoder.encoders import (
    OrbitEncoderFactory, SortEncoderFactory, HalfSortEncoderFactory,
    standard_row_pair_factories, OverlapBlockEncoderFactory, Phase2EncoderFactory,
)
from symcoder.experiments.encoding_fidelity_search import canonical_invariant
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


def _orbit_factory():
    return OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])


def _phase2_factory():
    return Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories(),
                                   use_complementarity_drop=False)
    ])


def _encode_default(plan, event):
    """Top-level encode with default-on Phase 3 (simplicial-complex embedding).

    Phase 3 is enabled by ``symcoder.encode``'s default, so it is included
    here without an explicit factory argument.
    """
    return symcoder_encode(plan, event,
                           orbit_factory=_orbit_factory(),
                           phase2_factory=_phase2_factory())


def _encode_no_phase3(plan, event):
    """Top-level encode with Phase 3 explicitly disabled — the historical
    Phase 1 + Phase 2 only output, used by the no-Phase-3 regression case."""
    return symcoder_encode(plan, event,
                           orbit_factory=_orbit_factory(),
                           phase2_factory=_phase2_factory(),
                           phase3_factory=None)


def _to_event(event_dict):
    return {k: np.array(v, dtype=float) for k, v in event_dict.items()}


def _parity_flip(event):
    return {k: -v for k, v in event.items()}


@pytest.fixture(scope="module")
def plan():
    return _make_plan()


@pytest.mark.parametrize("shadow_idx,raw_event",
                         list(enumerate(SHADOW_EVENTS)),
                         ids=[f"shadow_{i:02d}" for i in range(len(SHADOW_EVENTS))])
def test_parity_flip_distinguishable(plan, shadow_idx, raw_event):
    """encode(E) and encode(P(E)) must differ by more than ATOL.

    Two-part claim:

      (i)  E is chiral: E and P(E) are NOT related by any rotation +
           label permutation, hence lie in distinct SO(3) x S_n orbits.
           Verified up front using the classical Gram + signed-determinants
           invariant (Weyl's first fundamental theorem for SO(n)).  If this
           ever fails for a fixture event, the test is vacuous on that event
           and is *skipped* (not failed) — Skip rather than assert because a
           failed assertion would be indistinguishable from a real encoder
           bug; ``pytest.skip`` is the honest signal that the fixture has
           gone stale and needs regeneration.

      (ii) encode(E) ≠ encode(P(E)).  The encoder, with the default-on
           Phase 3 simplicial-complex embedding of the full alignment table,
           must produce different outputs for chiral E and P(E).

    History: before Phase 3 was added, this test failed (xfail) on every
    shadow event because Phase 1 + Phase 2 alone are parity-blind on these
    structured inputs.  ``test_parity_flip_no_phase3_still_blind`` below
    pins the historical Phase 1 + Phase 2 behaviour as an xfail to make the
    regression explicit.
    """
    event   = _to_event(raw_event)
    flipped = _parity_flip(event)

    # (i) Chirality precondition: confirm E and P(E) are in DIFFERENT
    # SO(3) x S_n orbits using Weyl's complete invariant for SO(n).  If they
    # happen to share an orbit (event is achiral), the parity test is
    # vacuous on this event and we skip rather than fail.
    inv_pos = canonical_invariant(event,   group="SO3")
    inv_neg = canonical_invariant(flipped, group="SO3")
    if inv_pos == inv_neg:
        pytest.skip(
            f"shadow #{shadow_idx} is achiral (E is SO(3) x S_n-equivalent to "
            f"P(E)); parity-blindness test is vacuous here.  "
            f"Fixture should be regenerated if this happens — see "
            f"symcoder/tests/_parity_shadow_events.py."
        )

    # (ii) Encoder distinguishability: for a chiral E, an SO(3)-faithful
    # encoder must produce different outputs for E and P(E).  Phase 3 is
    # on by default in symcoder.encode.
    enc_pos = _encode_default(plan, event)
    enc_neg = _encode_default(plan, flipped)

    max_diff = float(np.max(np.abs(enc_pos - enc_neg)))
    assert max_diff > ATOL, (
        f"shadow #{shadow_idx}: event is chiral (SO(3) x S_n distinguishes "
        f"E from P(E)), yet the encoder is parity-blind on it.\n"
        f"  encode(E)  - encode(P(E)) max |diff| = {max_diff:.3e}\n"
        f"  required:                            >  {ATOL:.3e}\n"
        f"  E = {raw_event}"
    )


@pytest.mark.parametrize("shadow_idx,raw_event",
                         list(enumerate(SHADOW_EVENTS)),
                         ids=[f"shadow_{i:02d}" for i in range(len(SHADOW_EVENTS))])
@pytest.mark.xfail(
    strict=True,
    reason=(
        "Historical parity-blindness of Phase 1 + Phase 2 alone (no Phase 3) "
        "on structured events.  encode(E, phase3=None) == encode(P(E), phase3=None) "
        "for this shadow event even though E and its parity flip are in "
        "different SO(3) orbits.  This test pins the historical bug as an "
        "xfail to make explicit which encoder weakness Phase 3 was added to "
        "address.  See module docstring."
    ),
)
def test_parity_flip_no_phase3_still_blind(plan, shadow_idx, raw_event):
    """With Phase 3 explicitly disabled the encoder must still be parity-blind.

    This is the *companion* to test_parity_flip_distinguishable: the latter
    confirms that Phase 3 (default-on) restores parity-faithfulness; this
    one pins the historical Phase-1+2-only behaviour so the regression
    remains visible.  Both are needed to document the fix.
    """
    event   = _to_event(raw_event)
    flipped = _parity_flip(event)

    # Same chirality precondition as the main test.
    inv_pos = canonical_invariant(event,   group="SO3")
    inv_neg = canonical_invariant(flipped, group="SO3")
    if inv_pos == inv_neg:
        pytest.skip(
            f"shadow #{shadow_idx} is achiral; no-Phase-3 regression vacuous."
        )

    enc_pos = _encode_no_phase3(plan, event)
    enc_neg = _encode_no_phase3(plan, flipped)

    max_diff = float(np.max(np.abs(enc_pos - enc_neg)))
    assert max_diff > ATOL, (
        f"shadow #{shadow_idx}: Phase 1 + Phase 2 alone is parity-blind here "
        f"(expected — this test is xfail-pinned).  max |diff| = {max_diff:.3e}"
    )
