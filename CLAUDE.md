# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

The venv lives one level up from this repo. Activate it or use it directly:

```bash
source ../venv/bin/activate
```

Run all tests:
```bash
pytest
# or with stdout visible:
pytest -s
```

Run a single test file:
```bash
pytest MinimalConfusableSets/test_bistate_vertex_matches.py
```

Run a single test by name:
```bash
pytest test_EncDec.py::test_name
```

Run the example script:
```bash
python example.py
```

**Note on `MinimalConfusableSets/` tests:** 18 tests in that subdirectory currently fail with `AssertionError: not (cfg.use_tristate and cfg.use_bistate)` — this is a pre-existing bug where both flags in `MinimalConfusableSets/config.py` are being set `True` simultaneously, which the code forbids. The 62 top-level tests all pass.

## Architecture

### Core abstractions

**`MultisetEncoder`** (`MultisetEncoder.py`) — base class for encoders. An encoder maps a multiset of n k-vectors (an `(n,k)` numpy array) to a 1D float array. Not necessarily injective.

**`MultisetEmbedder`** (`MultisetEmbedder.py`) — subclass of `MultisetEncoder` that is additionally injective (a true embedding). Handles edge cases (n=0, k=0, n=1, k=1) and dispatches to `embed_kOne` or `embed_generic` for the interesting cases. All concrete embedders subclass this.

Each concrete embedder lives in its own file named after its mathematical properties:
- `C0*` = piecewise linear (C⁰), `Cinf*` = infinitely differentiable
- `*array*` = accepts 2D arrays of k-vectors, `*list*` = accepts 1D lists (k=1 only)
- `*multiset*` = embeds symmetric product space, `*projectivespace*` = embeds ℝᵏ/Z₂

### `EncDec.py` — the Simplex preprocessing pipeline

Contains `simplex_1_preprocess_steps` and `simplex_2_preprocess_steps`, which implement the core mathematics behind the two simplicial complex embedders. The pipeline is:

1. Represent the input array as a `LinComb` of `MonoLinComb` basis elements (Eji vectors)
2. Barycentric-subdivide once (first differences) via `barycentric_subdivide()`
3. Barycentric-subdivide again (second differences)
4. Optionally canonicalise by sorting rows lexicographically — this makes the output permutation-invariant

`LinComb` / `MonoLinComb` are lightweight algebraic types (linear combination over numpy array basis vectors) defined in this file and used throughout the Simplex encoders.

### `MinimalConfusableSets/` — separate sub-project

A research sub-project studying which multisets are "confusable" by the dotting encoder. It has its own `config.py` with global flags (`use_bistate`, `use_tristate`) that must not both be `True`. The key entry points are `vertex_matches.py` (generating vertex match signatures) and `confusable_multisets.py`. Tests in this subdirectory should be run via `MinimalConfusableSets/run_tests.sh` rather than bare `pytest` due to global config state.

### Key utility modules

- `tools.py` — numpy helpers (e.g. `sort_np_array_rows_lexicographically`, `expand_complex_to_real_pairs`)
- `distinct_permutations_with_leftovers.py`, `distinct_partitions_with_start.py` — combinatorial generators used by the encoders and the `MinimalConfusableSets` sub-project
- `injection.py` — injectivity testing utilities
- `brute_force_array_decoder.py` — brute-force decoder for validation

### Embedding sizes (order)

All embedders implement `size_from_n_k(n, k) -> int`. The base class handles n≤1 or k≤1 (always returns nk); the generic case is delegated to `size_from_n_k_generic`. Embedding sizes determine output array length and are checked via assertions inside `embed()`.
