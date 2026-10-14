# symembed

Continuous, algebraic embeddings of `n` vectors in `k` dimensions that are invariant under a labelling group `G` and under
rotations, and that separate orbits on **every** input.

- `G` is a direct product of symmetric groups, one per species, optionally times parity `P`. Parity reflects every
  vector, so it negates every determinant.
- The specification and the mathematics are in [`docs/invariant_multisets.tex`](docs/invariant_multisets.tex). Section
  numbers below refer to it.

## Status

Only the **unpruned Bursarial baseline** is implemented (certification Level 0, Section 4.3).

- Layer 1 (Section 1): `w` holds all the dots `d` and determinants `e`, and `G` acts on it by signed permutations.
- Layer 2 (Section 3): the Bursarial embedding of the multiset `{g·w}`.
- Layer 3 (reductions): not implemented, so `reductions="none"` is the only option.
- `homogenise` and `root_normalise` are not implemented either.

The baseline is a reference against which reductions can be validated, not something to train on. Its outputs have
degree up to `|G|` in `w`, and they overflow easily. For `n = 5`, `k = 4`, `G = S(5)` (`D = 138060`) and unit-scale
Gaussian inputs:

- the largest outputs reach about `1e191` in float64;
- float32 gives `inf`/`nan`.

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

## Layout

| file | contents |
|---|---|
| `spec.py` | `GroupSpec`, `EmbeddingOptions` |
| `slots.py` | slots of `w`, orbit grouping, group elements, signed-permutation table (pure Python, independent of symatom) |
| `bursar.py` | Bursarial coefficients from `{g·w}` in torch |
| `module.py` | `InvariantEmbedding` (the evaluator) and its flattened variant |
| `meta.py` | `EmbeddingMeta`: JSON, LaTeX, symbolic expansion, possible degrees in `L` |
| `build.py` | `build_invariant_embedding` |

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
- comparison against full symbolic expansion;
- a cross-check against the legacy Bursarial embedder;
- invariance under permutations, rotations, and reflections (when `P` is in `G`);
- shapes and dtypes;
- `gradcheck`;
- the metadata.
