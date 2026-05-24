"""
symcoder.tests.test_alignment_decoder
======================================
Roundtrip tests for the alignment decoder.

The alignment decoder (decode_alignment) recovers the G-orbit of repS
evaluations from Phase 1 and Phase 2 decoded outputs.

Ground truth
------------
For each group element g ∈ G (enumerated via TheGroup.all_group_elements()),
we evaluate eval(g·repS[r], E) for all r = 0...|repS|-1.  This gives
n_cols = |G| column vectors.  The decoded multiset must match this ground truth.

Phase 2 factory: complementarity drop DISABLED (phase2_factory_no_drop)
------------------------------------------------------------------------
All tests here use phase2_factory_no_drop rather than the standard
phase2_factory.  The reason is a fundamental constraint on ANTISYMMETRIC
cross-FO pairs:

  When complementarity drop is enabled, a cross-FO block that contains a
  single PairFlavour (e.g. (eps3, 2e+1m) × (eps3, 3e) with forced overlap=2)
  has that sole PairFlavour comp-dropped, giving the block zero output
  dimensions.  The NULL_COMP reconstruction then returns U × V (the full
  Cartesian product of the two row-value lists) — 24 pairs for the context
  here — instead of the 12-pair joint-orbit required by the alignment decoder's
  pair-relation filter (len(full_pairs) == n_cols).

  Without the direct ASSOC pair (eps3(a,b,p), eps3(a,b,c)) being encoded, no
  available pair constraint can uniquely determine the sign of eps3(a,b,c) for
  each column.  This is because the stabiliser of every other repS atom always
  contains an odd electron permutation that flips eps3(a,b,c)'s sign, making
  all such constraints degenerate.

  With complementarity drop disabled the ASSOC pair is encoded directly (12
  pairs, stab=1, full_pairs==n_cols), resolving the sign unambiguously.

Modes
-----
Each test runs in both float and integer event modes via the parametrised
``event_mode`` fixture in conftest.py (→ both modes on a single pytest run).

Float mode  (event_mode="float"):
    Continuous N(0,1) vectors; polished decoded values match ground truth to
    within atol_exact = 1e-10.

Integer mode  (event_mode="integer"):
    Exact-integer vectors from [-4, 4]; polished decoded values are bitwise
    equal to ground truth (atol_exact = 0.0).
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pytest
from symatom import Plan, repS as _repS_fn
from symcoder.alignment_decoder import decode_alignment
from symcoder.decoded_types import AnnotatedMultisetOfRepSEvalVectors
from symcoder.eval import evaluate


# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------

def _build_phase1_results(orbit_enc, phase1_vals):
    """Decode all Phase 1 orbits; return dict keyed by (operation, flavour_counts)."""
    results = {}
    cursor  = 0
    for fo, enc in orbit_enc._selections:
        chunk = phase1_vals[cursor:cursor + enc.output_dim]
        results[(fo.operation, tuple(fo.flavour.counts))] = enc.decode(chunk)
        cursor += enc.output_dim
    return results


def _build_all_pair_decoded(phase2_enc, phase2_vals, phase1_results):
    """Run block.decode() on every block; return flat list of AnnotatedMultisetOfRealPairs."""
    all_decoded = []
    cursor      = 0
    for _, block_enc in phase2_enc._block_encoders:
        block_slice  = phase2_vals[cursor:cursor + block_enc.output_dim]
        decoded_list = block_enc.decode(block_slice, phase1_results)
        all_decoded.extend(decoded_list)
        cursor += block_enc.output_dim
    return all_decoded


def _ground_truth_alignment(plan, repS_atoms, event):
    """
    Compute the n_cols = |G| ground-truth column vectors.

    For each g ∈ G (via TheGroup.all_group_elements()), column c = tuple of
    eval(g·repS[r], E) for r = 0...|repS|-1.  Duplicates are kept: for a
    non-trivial stabiliser some group elements produce identical columns.
    """
    the_group = plan.context.the_group
    vectors   = []
    for g in the_group.all_group_elements():
        col = tuple(
            evaluate(g.apply(atom), event)
            for atom in repS_atoms
        )
        vectors.append(col)
    return vectors


def _multisets_close(vecs_a, vecs_b, atol):
    """Greedy order-insensitive comparison: each element of vecs_a must be
    close (max-norm ≤ atol) to some unused element of vecs_b.

    Returns True iff the two lists are equal as approximate multisets.
    """
    if len(vecs_a) != len(vecs_b):
        return False
    a    = [np.array(v, dtype=float) for v in vecs_a]
    b    = [np.array(v, dtype=float) for v in vecs_b]
    used = [False] * len(b)
    for va in a:
        matched = False
        for i, vb in enumerate(b):
            if not used[i] and np.max(np.abs(va - vb)) <= atol:
                used[i]  = True
                matched  = True
                break
        if not matched:
            return False
    return True


# ---------------------------------------------------------------------------
# Main roundtrip test
# ---------------------------------------------------------------------------

def test_alignment_decoder_roundtrip(
    ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol_exact
):
    """Full encode → decode Phase1 + Phase2 → alignment decoder → compare multiset.

    Verifies that the decoded multiset of column vectors equals the ground-truth
    multiset produced by evaluating all n_cols = |G| group elements on repS.
    """
    mag, dot, eps3 = ops
    plan       = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc  = orbit_factory.build(plan)
    phase2_enc = phase2_factory_no_drop.build(plan)

    phase1_vals = orbit_enc.encode(event).values
    phase2_vals = phase2_enc.encode(event).values

    phase1_results  = _build_phase1_results(orbit_enc, phase1_vals)
    all_pair_decoded = _build_all_pair_decoded(phase2_enc, phase2_vals, phase1_results)

    decoded = decode_alignment(plan, phase1_results, all_pair_decoded, atol=atol_exact)

    assert isinstance(decoded, AnnotatedMultisetOfRepSEvalVectors)

    # Output should contain exactly n_cols vectors.
    n_cols = plan.context.the_group.order()
    assert len(decoded.vectors) == n_cols, (
        f"Expected {n_cols} column vectors, got {len(decoded.vectors)}"
    )

    # atoms field should list the repS canonical representatives.
    fo_list    = _repS_fn(plan.context, plan.operations)
    repS_atoms = [fo.canonical_representative() for fo in fo_list]
    assert decoded.atoms == repS_atoms

    # Ground truth: evaluate all group elements on repS.
    truth = _ground_truth_alignment(plan, repS_atoms, event)
    assert len(truth) == n_cols

    # Multiset comparison — order-insensitive.
    if atol_exact == 0.0:
        # Integer mode: polished values are bitwise exact → use Counter equality.
        decoded_counter = Counter(decoded.vectors)
        truth_counter   = Counter(tuple(float(v) for v in vec) for vec in truth)
        assert decoded_counter == truth_counter, (
            f"Alignment decoder exact-multiset mismatch.\n"
            f"Decoded unique vectors ({len(decoded_counter)} distinct):\n"
            + "\n".join(f"  {k}: ×{v}" for k, v in sorted(decoded_counter.items()))
            + f"\nTruth unique vectors ({len(truth_counter)} distinct):\n"
            + "\n".join(f"  {k}: ×{v}" for k, v in sorted(truth_counter.items()))
        )
    else:
        # Float mode: approximate comparison (atol = 1e-10 after polishing).
        assert _multisets_close(decoded.vectors, truth, atol=atol_exact), (
            f"Alignment decoder approximate-multiset mismatch (atol={atol_exact}).\n"
            f"Decoded ({len(decoded.vectors)} vectors): "
            f"{sorted(decoded.vectors)[:5]}{'...' if len(decoded.vectors) > 5 else ''}\n"
            f"Truth   ({len(truth)} vectors): "
            f"{sorted(truth)[:5]}{'...' if len(truth) > 5 else ''}"
        )


# ---------------------------------------------------------------------------
# Structural / metadata tests
# ---------------------------------------------------------------------------

def test_alignment_decoder_vector_length(
    ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol_exact
):
    """Every decoded column vector has length |repS|."""
    mag, dot, eps3 = ops
    plan       = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc  = orbit_factory.build(plan)
    phase2_enc = phase2_factory_no_drop.build(plan)

    phase1_vals = orbit_enc.encode(event).values
    phase2_vals = phase2_enc.encode(event).values

    phase1_results   = _build_phase1_results(orbit_enc, phase1_vals)
    all_pair_decoded = _build_all_pair_decoded(phase2_enc, phase2_vals, phase1_results)

    decoded = decode_alignment(plan, phase1_results, all_pair_decoded, atol=atol_exact)

    n_reps = len(decoded.atoms)
    for vec in decoded.vectors:
        assert len(vec) == n_reps, (
            f"Vector length {len(vec)} ≠ |repS| = {n_reps}"
        )


def test_alignment_decoder_row_multisets_correct(
    ops, ctx, event, orbit_factory, phase2_factory_no_drop, atol_exact
):
    """For each repS atom, the multiset of its row values in the decoded table
    matches the Phase 1 decoded values (replicated by stabiliser factor).

    This checks that the row-level structure is preserved regardless of how
    pair constraints assigned values to individual columns.
    """
    mag, dot, eps3 = ops
    plan       = Plan(context=ctx, operations=(mag, dot, eps3))
    orbit_enc  = orbit_factory.build(plan)
    phase2_enc = phase2_factory_no_drop.build(plan)

    phase1_vals = orbit_enc.encode(event).values
    phase2_vals = phase2_enc.encode(event).values

    phase1_results   = _build_phase1_results(orbit_enc, phase1_vals)
    all_pair_decoded = _build_all_pair_decoded(phase2_enc, phase2_vals, phase1_results)

    decoded = decode_alignment(plan, phase1_results, all_pair_decoded, atol=atol_exact)

    n_cols = plan.context.the_group.order()
    fo_list = _repS_fn(plan.context, plan.operations)

    for r, (fo, atom) in enumerate(zip(fo_list, decoded.atoms)):
        key        = (fo.operation, tuple(fo.flavour.counts))
        phase1     = phase1_results[key]
        orbit_size = len(phase1.values)
        stab_u     = n_cols // orbit_size
        expected   = sorted(phase1.values * stab_u)

        # Extract row r from the decoded table.
        row_vals = sorted(vec[r] for vec in decoded.vectors)

        if atol_exact == 0.0:
            assert row_vals == expected, (
                f"Row {r} ({fo!r}): row multiset {row_vals} ≠ Phase1 expected {expected}"
            )
        else:
            assert np.allclose(row_vals, expected, atol=atol_exact), (
                f"Row {r} ({fo!r}): row multiset {row_vals} "
                f"not close to Phase1 expected {expected} (atol={atol_exact})"
            )
