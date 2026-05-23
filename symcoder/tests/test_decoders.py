"""
Roundtrip decode tests for all leaf-level encoders (Step A, Phase 1 and Phase 2).

Strategy: encode a known event, then call decode() on each encoder's slice of
the output and verify that the decoded value multiset matches the ground-truth
values obtained by evaluating the atoms directly.

Each test runs in two modes — float and integer — via the parametrised
``event_mode`` fixture in conftest.py.  Tolerances are supplied by the
``atol`` and ``atol_exact`` fixtures, which adjust automatically:

float mode  (event_mode="float"):
    N(0,1) random vectors; atol = 1e-10.

integer mode  (event_mode="integer"):
    Exact-integer vectors from [-4, 4], including zeros (frequent in
    real-world inputs).  atol = 1e-6 for unpolished outputs (double-root
    splitting from Phase 1 collisions reaches ~2e-7).  atol_exact = 0.0
    for polished outputs (bitwise exact after polishing replaces noisy
    roots with the exact Phase 1 float).
"""
from __future__ import annotations

import numpy as np
import pytest
from symatom import Plan
from symcoder.encoders.row_pair_encoders import (
    Type11PairEncoder, Type12PairEncoder, Type21PairEncoder, Type22PairEncoder,
    NegPairEncoder, SelfPairType11Encoder, SelfPairNegEncoder,
)
from symcoder.decoded_types import AnnotatedMultisetOfReals, AnnotatedMultisetOfRealPairs
from symcoder.eval import evaluate


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sorted_vals(vals):
    return sorted(vals)


def _multisets_close(pairs_a, pairs_b, atol):
    """Order-insensitive check: every pair in pairs_a has a close match in pairs_b.

    Uses Chebyshev distance for matching so that both components must be within
    atol.  Runs a greedy O(n²) bipartite match — fine for small n.
    """
    if len(pairs_a) != len(pairs_b):
        return False
    a = np.array(pairs_a, dtype=float).reshape(-1, 2)
    b = np.array(pairs_b, dtype=float).reshape(-1, 2)
    used = np.zeros(len(b), dtype=bool)
    for row in a:
        dists = np.max(np.abs(b - row), axis=1)
        dists[used] = np.inf
        best = int(np.argmin(dists))
        if dists[best] > atol:
            return False
        used[best] = True
    return True


def _ground_truth_orbit(enc, event):
    """Evaluate all atoms in a Phase 1 encoder and return sorted values."""
    # SortEncoder has _orbit; HalfSortEncoder has _representatives + negatives
    from symcoder.encoders.sort_encoder import SortEncoder, HalfSortEncoder
    if isinstance(enc, SortEncoder):
        return sorted(evaluate(a, event) for a in enc._orbit)
    if isinstance(enc, HalfSortEncoder):
        pos = [evaluate(a, event) for a in enc._representatives]
        return sorted(pos + [-v for v in pos])
    raise TypeError(f"Unknown encoder {type(enc)}")


def _ground_truth_pairs(enc, event):
    """Evaluate all atom-pairs in a Phase 2 row-pair encoder and return sorted pairs."""
    from symcoder.encoders.row_pair_encoders import (
        Type11PairEncoder, Type12PairEncoder, Type21PairEncoder,
        Type22PairEncoder, NegPairEncoder,
        SelfPairType11Encoder, SelfPairNegEncoder,
    )
    from symatom.atoms import Atom

    def neg(a):
        return Atom(a.operation, a.labels, sign=-a.sign)

    reps = enc._reps

    if isinstance(enc, (Type11PairEncoder, SelfPairType11Encoder)):
        atom_pairs = list(reps)
    elif isinstance(enc, Type21PairEncoder):
        atom_pairs = list(reps) + [(neg(u), v) for u, v in reps]
    elif isinstance(enc, Type12PairEncoder):
        atom_pairs = list(reps) + [(u, neg(v)) for u, v in reps]
    elif isinstance(enc, Type22PairEncoder):
        atom_pairs = [(u, v) for u, v in reps]
        atom_pairs += [(neg(u), v) for u, v in reps]
        atom_pairs += [(u, neg(v)) for u, v in reps]
        atom_pairs += [(neg(u), neg(v)) for u, v in reps]
    elif isinstance(enc, (NegPairEncoder, SelfPairNegEncoder)):
        atom_pairs = list(reps) + [(neg(u), neg(v)) for u, v in reps]
    else:
        raise TypeError(f"Unknown encoder {type(enc)}")

    return [(evaluate(u, event), evaluate(v, event)) for u, v in atom_pairs]


# ---------------------------------------------------------------------------
# Phase 1 roundtrip tests
# ---------------------------------------------------------------------------

def test_sort_decoder_roundtrip(ops, ctx, event, orbit_factory, atol):
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc = orbit_factory.build(plan)
    encoded = orbit_enc.encode(event).values

    cursor = 0
    for fo, enc in orbit_enc._selections:
        from symcoder.encoders.sort_encoder import SortEncoder
        if not isinstance(enc, SortEncoder):
            cursor += enc.output_dim
            continue
        chunk = encoded[cursor:cursor + enc.output_dim]
        decoded = enc.decode(chunk)
        assert isinstance(decoded, AnnotatedMultisetOfReals)
        assert np.allclose(_sorted_vals(decoded.values),
                           _ground_truth_orbit(enc, event), atol=atol)
        cursor += enc.output_dim


def test_halfsort_decoder_roundtrip(ops, ctx, event, orbit_factory, atol):
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc = orbit_factory.build(plan)
    encoded = orbit_enc.encode(event).values

    cursor = 0
    found = False
    for fo, enc in orbit_enc._selections:
        from symcoder.encoders.sort_encoder import HalfSortEncoder
        chunk = encoded[cursor:cursor + enc.output_dim]
        if isinstance(enc, HalfSortEncoder):
            found = True
            decoded = enc.decode(chunk)
            assert isinstance(decoded, AnnotatedMultisetOfReals)
            assert len(decoded.values) == 2 * enc.output_dim
            assert np.allclose(_sorted_vals(decoded.values),
                               _ground_truth_orbit(enc, event), atol=atol)
        cursor += enc.output_dim
    assert found, "No HalfSortEncoder found in plan — check fixtures"


# ---------------------------------------------------------------------------
# Phase 2 row-pair roundtrip tests
# ---------------------------------------------------------------------------

def _roundtrip_all_pair_encoders(enc_type, ops, ctx, event, orbit_factory, phase2_factory, atol):
    """Generic roundtrip for a given PairOrbitEncoder subclass."""
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc = orbit_factory.build(plan)
    phase2_enc = phase2_factory.build(plan)

    phase1_vals = orbit_enc.encode(event).values
    phase2_vals = phase2_enc.encode(event).values
    full_encoded = np.concatenate([phase1_vals, phase2_vals])
    phase1_dim = orbit_enc.output_dim

    found = False
    cursor = phase1_dim
    for _, block_enc in phase2_enc._block_encoders:
        for sel in block_enc._selections:
            enc = sel.encoder
            length = enc.output_dim
            if sel.is_comp_drop or length == 0:
                continue
            if isinstance(enc, enc_type):
                found = True
                chunk = full_encoded[cursor:cursor + length]
                decoded = enc.decode(chunk)
                assert isinstance(decoded, AnnotatedMultisetOfRealPairs)
                truth = _ground_truth_pairs(enc, event)
                assert _multisets_close(decoded.pairs, truth, atol=atol), (
                    f"{enc_type.__name__}: decoded {sorted(decoded.pairs)} "
                    f"!= truth {sorted(truth)}"
                )
            cursor += length

    return found


def test_type11_decoder_roundtrip(ops, ctx, event, orbit_factory, phase2_factory, atol):
    found = _roundtrip_all_pair_encoders(
        Type11PairEncoder, ops, ctx, event, orbit_factory, phase2_factory, atol)
    assert found, "No Type11PairEncoder ASSOC found — check fixtures"


def test_type12_decoder_roundtrip(ops, ctx, event, orbit_factory, phase2_factory, atol):
    found = _roundtrip_all_pair_encoders(
        Type12PairEncoder, ops, ctx, event, orbit_factory, phase2_factory, atol)
    assert found, "No Type12PairEncoder ASSOC found — check fixtures"


def test_type21_decoder_roundtrip(ops, ctx, event, orbit_factory, phase2_factory, atol):
    found = _roundtrip_all_pair_encoders(
        Type21PairEncoder, ops, ctx, event, orbit_factory, phase2_factory, atol)
    assert found, "No Type21PairEncoder ASSOC found — check fixtures"


def test_type22_decoder_roundtrip(ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol):
    found = _roundtrip_all_pair_encoders(
        Type22PairEncoder, ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol)
    assert found, "No Type22PairEncoder found — check fixtures"


def test_neg_decoder_roundtrip(ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol):
    found = _roundtrip_all_pair_encoders(
        NegPairEncoder, ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol)
    assert found, "No NegPairEncoder found — check fixtures"


def test_selfpair_type11_sigma_decoder_roundtrip(ops, ctx, event, orbit_factory, phase2_factory, atol):
    found = _roundtrip_all_pair_encoders(
        SelfPairType11Encoder, ops, ctx, event, orbit_factory, phase2_factory, atol)
    assert found, "No SelfPairType11Encoder ASSOC found — check fixtures"


def test_selfpair_neg_sigma_decoder_roundtrip(ops, ctx, event, orbit_factory, phase2_factory, atol):
    found = _roundtrip_all_pair_encoders(
        SelfPairNegEncoder, ops, ctx, event, orbit_factory, phase2_factory, atol)
    assert found, "No SelfPairNegEncoder ASSOC found — check fixtures"


# ---------------------------------------------------------------------------
# OverlapBlock-level decoder test (complementarity drop disabled)
# ---------------------------------------------------------------------------

def _build_phase1_results(orbit_enc, phase1_vals):
    """Decode all Phase 1 orbits and return a dict keyed by (op_name, flavour_counts)."""
    results = {}
    cursor = 0
    for fo, enc in orbit_enc._selections:
        chunk = phase1_vals[cursor:cursor + enc.output_dim]
        results[(fo.operation.name, tuple(fo.flavour.counts))] = enc.decode(chunk)
        cursor += enc.output_dim
    return results


def test_overlap_block_decoder_no_drop(ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol):
    """OverlapBlockEncoder.decode() roundtrip with complementarity drop disabled.

    With use_complementarity_drop=False every block contains only ASSOC and
    NULL_SELF selections (no NULL_COMP).  The block decoder must:
      - slice the encoded array correctly for each ASSOC encoder and call its decode()
      - reconstruct NULL_SELF pairs from Phase 1 decoded values
    """
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc = orbit_factory.build(plan)
    phase2_enc = phase2_factory_no_drop.build(plan)

    phase1_vals = orbit_enc.encode(event).values
    phase2_vals = phase2_enc.encode(event).values

    phase1_results = _build_phase1_results(orbit_enc, phase1_vals)

    cursor = 0
    for _, block_enc in phase2_enc._block_encoders:
        block_slice = phase2_vals[cursor:cursor + block_enc.output_dim]
        decoded_list = block_enc.decode(block_slice, phase1_results)

        active_sels = [s for s in block_enc._selections if not s.is_comp_drop]
        assert len(decoded_list) == len(active_sels)

        for sel, decoded in zip(active_sels, decoded_list):
            assert isinstance(decoded, AnnotatedMultisetOfRealPairs)
            enc = sel.encoder

            if enc.output_dim == 0:
                # NULL_SELF: truth is (eval(atom), eval(atom)) for each Phase 1 atom
                key = (sel.pf.op_u.name, tuple(sel.pf.flavour_u.counts))
                truth = [(evaluate(a, event), evaluate(a, event))
                         for a in phase1_results[key].atoms]
            else:
                truth = _ground_truth_pairs(enc, event)

            assert _multisets_close(decoded.pairs, truth, atol=atol), (
                f"Block decode mismatch for {type(enc).__name__}: "
                f"decoded {sorted(decoded.pairs)} != truth {sorted(truth)}"
            )

        cursor += block_enc.output_dim


def test_overlap_block_decoder_with_drop(ops, ctx, event, orbit_factory, phase2_factory, atol):
    """OverlapBlockEncoder.decode() roundtrip with complementarity drop enabled (default).

    With use_complementarity_drop=True each block has exactly one NULL_COMP selection
    (the largest non-self association is dropped at encode time).  The block decoder must:
      - decode all non-dropped selections (ASSOC + NULL_SELF) as before
      - reconstruct the NULL_COMP pairs via multiset complement:
          Z_full(U × V) ∖ NULL_SELF ∖ all-ASSOC = NULL_COMP
      - return one result per selection (len == len(_selections), not len(active_sels))

    The complement works cleanly because _plan_owned_atom_pairs claims all seven
    sign/swap partner orbits for each selection (not just the base orbit), and
    decode() subtracts their value-space transforms from Z_full accordingly.
    """
    mag, dot, eps3 = ops
    plan = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc  = orbit_factory.build(plan)
    phase2_enc = phase2_factory.build(plan)

    phase1_vals = orbit_enc.encode(event).values
    phase2_vals = phase2_enc.encode(event).values

    phase1_results = _build_phase1_results(orbit_enc, phase1_vals)

    found_comp_drop = False
    cursor = 0
    for _, block_enc in phase2_enc._block_encoders:
        block_slice = phase2_vals[cursor:cursor + block_enc.output_dim]
        decoded_list = block_enc.decode(block_slice, phase1_results)

        # Result list must be one entry per selection (including the comp-dropped one).
        assert len(decoded_list) == len(block_enc._selections)

        for sel, decoded in zip(block_enc._selections, decoded_list):
            assert isinstance(decoded, AnnotatedMultisetOfRealPairs)
            enc = sel.encoder

            if enc.output_dim == 0:
                # NULL_SELF: truth = (eval(atom), eval(atom)) for each Phase 1 atom
                key = (sel.pf.op_u.name, tuple(sel.pf.flavour_u.counts))
                truth = [(evaluate(a, event), evaluate(a, event))
                         for a in phase1_results[key].atoms]
            elif sel.is_comp_drop:
                # NULL_COMP: ground truth = all owned atom-pairs evaluated.
                # owned_atom_pairs covers the base G-orbit plus all sign/swap
                # partner orbits claimed by this selection during planning.
                truth = [(evaluate(u, event), evaluate(v, event))
                         for u, v in sel.owned_atom_pairs]
                found_comp_drop = True
            else:
                # ASSOC: ground truth from encoder reps only (base orbit).
                truth = _ground_truth_pairs(enc, event)

            assert _multisets_close(decoded.pairs, truth, atol=atol), (
                f"Block decode mismatch for {type(enc).__name__} "
                f"({'NULL_COMP' if sel.is_comp_drop else 'ASSOC/NULL_SELF'}): "
                f"decoded {sorted(decoded.pairs)} != truth {sorted(truth)}"
            )

        cursor += block_enc.output_dim

    assert found_comp_drop, "No comp-dropped selection found — check fixtures or factory"
