# symembed

Continuous, algebraic embeddings of `n` vectors in `k` dimensions that are invariant under a labelling group `G` and under
rotations, and that separate orbits on **every** input.

- `G` is a direct product of symmetric groups, one per species, optionally times parity `P`. Parity reflects every
  vector, so it negates every determinant.
- The specification and the mathematics are in [`docs/invariant_multisets.tex`](docs/invariant_multisets.tex). Section
  numbers below refer to it.

## Status

Two levels are implemented (certification levels are defined in Section 4.3):

- **`reductions="none"`**: the unpruned Bursarial baseline (Level 0).
  - Layer 1 (Section 1): `w` holds all the dots `d` and determinants `e`, and `G` acts on it by signed permutations.
  - Layer 2 (Section 3): the Bursarial embedding of the multiset `{g·w}`.
- **`reductions="safe"`**: the reductions that need no judgement (Section 4.2).
  - With parity in `G`, the determinants are omitted (Level 0 for the smaller `w`). The multiset is then taken over
    the permutation part of `G`, since parity acts trivially on the `d`'s.
  - Outputs that are exact rational linear combinations of kept outputs, as polynomials in `w`, are dropped
    (Level 1). The `E_{c,j}` are expanded exactly, as integer polynomials, for `c <= linear_reduction_max_c`. This is
    chosen automatically from a work budget, because the number of monomials grows fast. Candidate relations are
    found modulo a prime and then verified exactly. Anything not verified is kept, so a failure can only cost
    compression, never injectivity. Every relation is recorded in `meta` and can be replayed with
    `meta.reconstruct_base`.

`homogenise`, `root_normalise`, the complex field and Level 2 reductions are not implemented.

Exact linear dependencies are rare, so the safe reductions remove very little. For `n = 5`, `k = 4`, `G = S(5)`:

| | Bursarial outputs | after safe reductions | dropped |
|---|---|---|---|
| without `P` | 138060 | 137979 | 81 (checked for `c <= 4`) |
| with `P` (determinants omitted) | 101760 | 101701 | 59 (checked for `c <= 5`, none found at `c = 5`) |

Almost all the redundancy is therefore nonlinear, and needs Level 2.

The unreduced outputs are a reference to validate reductions against, not something to train on. Their degree in `w`
goes up to the multiset size, and they overflow easily. At `n = 5`, `k = 4`, unit-scale Gaussian inputs give
outputs up to about `1e191` in float64, and `inf`/`nan` in float32.

## Use

```python
import torch
from symembed import GroupSpec, EmbeddingOptions, build_invariant_embedding

spec = GroupSpec(species=(("a", "b"), ("p",)), parity=False)       # S(2) x S(1)
embed, meta = build_invariant_embedding(k=2, group_spec=spec, options=EmbeddingOptions(dtype=torch.float64))

X = torch.randn(1000, 3, 2, dtype=torch.float64)   # (..., n, k); X[..., i, :] has the i-th label of spec.labels
Y = embed(X)                                       # (..., meta.D), differentiable

meta.to_json()            # machine-readable description (slots, group, signed permutations, outputs, certification)
meta.latex()              # LaTeX tables (needs amsmath + longtable)
meta.expanded_output(21)  # small cases: output 21 expanded in the d's and e's (sympy)
```

- `embed` is a `torch.nn.Module` following IWPC's `Encoding` conventions: `input_shape` and `output_shape` are integer
  buffers, and it acts on tensors of shape `(..., *input_shape)`. It does not depend on IWPC.
- `embed.flattened()` takes input of shape `(..., n*k)` instead, with component `c` of vector `i` at position `i*k + c`.

## Exploration

`python -m symembed.explore.simple_cases` (run from the repo root, with the venv's Python) writes very simple cases
out in full into `docs/explorations/`, for inspection by eye:

- `simple_cases.pdf` shows the slots of `w`, the action of `G`, and every output expanded in the `d`'s and `e`'s,
  marked kept or dropped, with the identity used to drop it. Expansions longer than 100 terms show only their term
  count.
- `simple_cases_full.txt` has every expansion, with no limit.

The cases are S(2) on (a, b), S(3) on (a, b, c) and S(2) x S(1) on (a, b, p), each with k = 2 and with k = 1.

## Layout

| file | contents |
|---|---|
| `spec.py` | `GroupSpec`, `EmbeddingOptions` |
| `slots.py` | slots of `w`, orbit grouping, group elements, signed-permutation table (pure Python, independent of symatom) |
| `bursar.py` | Bursarial coefficients from `{g·w}` in torch |
| `reductions.py` | exact expansion of low-`c` outputs, exact linear-dependency search (Level 1) |
| `module.py` | `InvariantEmbedding` (the evaluator) and its flattened variant |
| `meta.py` | `EmbeddingMeta`: JSON, LaTeX, symbolic expansion, possible degrees in `L` |
| `build.py` | `build_invariant_embedding` |
| `explore/` | scripts that write cases out in full for inspection |

## Conventions

- `q_g(x) = sum_j (g·w)_j x^j`: slot 0 multiplies `x^0`. The legacy
  `Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset.py` puts component 0 on the highest power instead. The
  tests allow for this when cross-checking against it.
- Outputs are `E_{c,j}` for `c = 1..|G|` and `j = 0..c(m-1)`, ordered by `c` and then `j`.

## Tests

```
/Users/lester/github/Echinocoder/venv/bin/pytest symembed
```

The tests cover:
- the worked example of the notes, symbol for symbol;
- the exact expansion against sympy, and the linear algebra on synthetic cases, including a relation too large to
  reconstruct, which must be kept;
- the safe reductions, by replaying the dropped outputs and comparing them numerically and symbolically;
- comparison against full symbolic expansion;
- a cross-check against the legacy Bursarial embedder;
- invariance under permutations, rotations, and reflections (when `P` is in `G`);
- shapes and dtypes;
- `gradcheck`;
- the metadata.
