"""
symcoder.encoders.embed_compressed_encoder
==========================================
Concrete Phase 2 encoder wrapping the legacy _embed_compressed path.

This bridges the existing algebraic pair-orbit embedding machinery into the
PairOrbitEncoder / PairOrbitEncoderFactory registry framework, making Phase 2
encoding pluggable in the same way that Phase 1 encoding is.

EmbedCompressedEncoder
    Ready-to-use encoder for one PairFlavour.  Pre-binds pf, plan, and
    type_sizes at assess time so that encode(event) needs only event data.
    output_dim = 2 * pf.count(type_sizes) reals (real + imaginary parts of n
    complex polynomial coefficients).

EmbedCompressedEncoderFactory
    Factory: accepts any non-self PairOrbitSpec (self-pairs produce no output
    and are handled by the NULL_SELF drop logic in the orchestrator).
    Priority 0.5 (same as SortEncoderFactory — the only current option).
"""
from __future__ import annotations

import numpy as np

from symatom.context import Plan
from symatom.rep import PairFlavour

from .pair_base import (
    PairOrbitEncoder,
    PairOrbitEncoderFactory,
    PairOrbitSpec,
    EncodingResult,
)

# Import the two core Phase 2 computation functions from encode.py.
# We use a lazy import to avoid a circular dependency: encode.py imports from
# encoders/, and embed_compressed_encoder.py would otherwise import from encode.py.
# The lazy import fires only the first time encode() is called on an instance,
# by which point all modules are fully initialised.
def _get_phase2_fns():
    """Lazily import the two Phase 2 functions from symcoder.encode."""
    from symcoder.encode import _embed_compressed, _complex_to_reals  # noqa: PLC0415
    from symcoder.pairs import eval_pair_orbit_positive               # noqa: PLC0415
    return eval_pair_orbit_positive, _embed_compressed, _complex_to_reals


_PRIORITY    = 0.5
_METHOD_NAME = "embed_compressed"


# ---------------------------------------------------------------------------
# EmbedCompressedEncoder — the ready-to-use encoder produced by the factory
# ---------------------------------------------------------------------------

class EmbedCompressedEncoder(PairOrbitEncoder):
    """
    Embeds a pair orbit using the algebraic _embed_compressed pipeline.

    Constructed by EmbedCompressedEncoderFactory.assess() with pf and plan
    pre-bound.  encode(event) requires only the event data.

    Algorithm (delegated to encode.py)
    ------------------------------------
    1. Evaluate z_pos = eval_pair_orbit_positive(pf, plan, event)  — the
       complex z-values for the sign=(+1,+1) subset of the pair orbit.
    2. Embed via _embed_compressed(z_pos, pf)  — produces n complex numbers
       using the appropriate SignCorrelationType branch.
    3. Unpack to 2n reals via _complex_to_reals.

    output_dim = 2 * pf.count(type_sizes)
    """

    def __init__(self, pf: PairFlavour, plan: Plan) -> None:
        self._pf         = pf
        self._plan       = plan
        self._type_sizes = tuple(g.size for g in plan.context.types)
        self._output_dim = 2 * pf.count(self._type_sizes)

    @property
    def output_dim(self) -> int:
        return self._output_dim

    @property
    def priority(self) -> float:
        return _PRIORITY

    @property
    def method_name(self) -> str:
        return _METHOD_NAME

    def encode(self, event: dict) -> EncodingResult:
        eval_pair_orbit_positive, _embed_compressed, _complex_to_reals = _get_phase2_fns()

        z_pos = np.array(
            eval_pair_orbit_positive(self._pf, self._plan, event),
            dtype=complex,
        )
        assert len(z_pos) == self._pf.count(self._type_sizes), (
            f"eval_pair_orbit_positive returned {len(z_pos)} values "
            f"but expected pf.count={self._pf.count(self._type_sizes)} for {self._pf!r}"
        )
        coeffs = _embed_compressed(z_pos, self._pf)
        values = _complex_to_reals(coeffs)
        assert len(values) == self._output_dim, (
            f"_complex_to_reals produced {len(values)} values "
            f"but output_dim={self._output_dim} for {self._pf!r}"
        )
        return EncodingResult(
            values=values,
            metadata={
                "method":     _METHOD_NAME,
                "pair_count": self._pf.count(self._type_sizes),
            },
        )


# ---------------------------------------------------------------------------
# EmbedCompressedEncoderFactory — creates EmbedCompressedEncoder instances
# ---------------------------------------------------------------------------

class EmbedCompressedEncoderFactory(PairOrbitEncoderFactory):
    """
    Factory that produces EmbedCompressedEncoder instances for non-self pair orbits.

    Self-pairs (where pf.op_u == pf.op_v and flavour_u == flavour_v with full
    overlap) are handled by the NULL_SELF drop logic in the orchestrator and
    should never reach the registry; this factory returns [] for them as a
    safety net.

    Returns a single EmbedCompressedEncoder for all other pairs.
    """

    def assess(self, spec: PairOrbitSpec, plan: Plan) -> list[PairOrbitEncoder]:
        from symcoder.pairs import _is_self_pair  # lazy, avoids circular at import time
        pf = spec.pf
        if _is_self_pair(pf):
            # Self-pairs are dropped before the registry is consulted; return []
            # here as a safety net so a misconfigured caller gets no encoder
            # rather than a wrong one.
            return []
        return [EmbedCompressedEncoder(pf, plan)]
