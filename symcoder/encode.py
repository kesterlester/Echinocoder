from __future__ import annotations
import numpy as np

from .encoders.orbit_encoder import OrbitEncoderFactory, OrbitEncoder
from .encoders.phase2_encoder import Phase2EncoderFactory
from .describe import SegmentInfo


def encode_and_describe(
    plan,
    event: dict | None,
    orbit_factory: OrbitEncoderFactory | None = None,
    phase2_factory: Phase2EncoderFactory | None = None,
) -> tuple[np.ndarray, list[SegmentInfo]]:
    """
    Core encoding function.  Runs the full two-phase encoding and simultaneously
    collects a SegmentInfo descriptor for every segment produced.

    Returns (values, segments) where:
      values   — 1-D float64 numpy array (the permutation-invariant embedding).
                 When event is None, values is a zero array of the correct length
                 (useful for describe-only callers that discard it).
      segments — list[SegmentInfo] mirroring the structure of values exactly.

    Passing event=None is the describe-only mode: segment lengths and metadata
    are computed identically (they are purely algebraic / factory-derived), but
    no atom evaluations are performed and the values array contains zeros.
    describe_encoding() uses this mode so that both functions share exactly one
    code path.

    Phase 1 (orbit_factory): an OrbitEncoderFactory built from sub-factories
      (e.g. SortEncoderFactory).  When None, Phase 1 is skipped entirely.
    Phase 2 (phase2_factory): a Phase2EncoderFactory built from the hierarchical
      OverlapBlockEncoderFactory → row-pair factories chain.  When None, Phase 2
      is skipped entirely.
    """
    parts    = []
    segments = []
    cursor   = 0

    # ------------------------------------------------------------------
    # Phase 1: encode each FlavouredOperator orbit.
    # ------------------------------------------------------------------
    if orbit_factory is not None:
        orbit_enc = orbit_factory.build(plan)
        if event is not None:
            parts.append(orbit_enc.encode(event).values)
        else:
            parts.append(np.zeros(orbit_enc.output_dim, dtype=np.float64))
        segments.extend(orbit_enc.describe(start_offset=0))
        cursor += orbit_enc.output_dim

    # ------------------------------------------------------------------
    # Phase 2: compressed pair encoding.
    # ------------------------------------------------------------------
    if phase2_factory is not None:
        phase2_enc = phase2_factory.build(plan)
        if event is not None:
            parts.append(phase2_enc.encode(event).values)
        else:
            parts.append(np.zeros(phase2_enc.output_dim, dtype=np.float64))
        segments.extend(phase2_enc.describe(start_offset=cursor))
        cursor += phase2_enc.output_dim

    values = np.concatenate(parts) if parts else np.array([], dtype=float)
    return values, segments


def encode(
    plan,
    event: dict,
    orbit_factory: OrbitEncoderFactory | None = None,
    phase2_factory: Phase2EncoderFactory | None = None,
) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector.
    See encode_and_describe() for full documentation of the two-phase algorithm.
    Returns a 1-D float64 numpy array.
    """
    values, _segments = encode_and_describe(plan, event, orbit_factory, phase2_factory)
    return values


def describe_encoding(
    plan,
    orbit_factory: OrbitEncoderFactory | None = None,
    phase2_factory: Phase2EncoderFactory | None = None,
) -> list[SegmentInfo]:
    """
    Return a list of SegmentInfo objects describing the full structure of the
    vector produced by encode(plan, event, orbit_factory, phase2_factory).

    Delegates to encode_and_describe(plan, event=None, ...) — the single
    authoritative code path — and discards the (zero) values array.
    """
    _, segments = encode_and_describe(plan, event=None, orbit_factory=orbit_factory, phase2_factory=phase2_factory)
    return segments
