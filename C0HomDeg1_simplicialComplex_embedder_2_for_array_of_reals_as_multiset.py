"""Algorithm 2 — wrapper that selects between the _2a_ and _2b_ implementations.

Two interchangeable implementations exist:

  _2a_  Original direct implementation (Eji / Maximal_Simplex_Vertex classes).
  _2b_  EncDec-backed implementation (uses EncDec.py as a service library).

Both produce numerically identical output under the default configuration.
Switch the active implementation by changing ACTIVE_IMPLEMENTATION below.

Current default: _2a_  (the original; known-good and well-tested).
Switch to _2b_ once the parity tests in test_simplex2_ab_parity.py pass and
you are satisfied with performance.
"""

# ── Choose the active implementation ─────────────────────────────────────────
#
# ACTIVE_IMPLEMENTATION = "2a"  ← original direct implementation (default)
# ACTIVE_IMPLEMENTATION = "2b"  ← EncDec-backed implementation
#
ACTIVE_IMPLEMENTATION = "2a"
# ─────────────────────────────────────────────────────────────────────────────

if ACTIVE_IMPLEMENTATION == "2a":
    from C0HomDeg1_simplicialComplex_embedder_2a_for_array_of_reals_as_multiset import (
        Embedder,
        eji_set_to_np_array,
        tost,
        run_unit_tests,
    )
elif ACTIVE_IMPLEMENTATION == "2b":
    from C0HomDeg1_simplicialComplex_embedder_2b_for_array_of_reals_as_multiset import (
        Embedder,
        tost,
        run_unit_tests,
    )
else:
    raise ValueError(f"Unknown ACTIVE_IMPLEMENTATION={ACTIVE_IMPLEMENTATION!r}; "
                     f"expected '2a' or '2b'.")
