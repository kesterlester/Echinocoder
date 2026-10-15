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

- **`embedder="bursar_involution"`**: involution packing (Section 3.3), which is Level 0.
  - An element `g0` of `G` acting on `w` as an involution, chosen with the smallest `|trace|`, packs `w` into
    `(m + |tr|)/2` complex coordinates on which `g0` acts as complex conjugation.
  - The multiset of packed vectors is then closed under conjugation, so the complex Bursarial polynomial has real
    coefficients, and only those are output. That roughly halves the size: for example 11 → 5, 174 → 90, and
    138060 → 79980 for `S(5)`, `k = 4`.
  - It combines with `reductions="safe"`.
  - The core, `involution.packing_from_matrix`, is generic in the involution `T`, so it applies to any multiset
    known to be invariant under a fixed linear involution.
  - On the simple cases, full polynomial closure leaves the same number of outputs for the packed and plain
    embeddings: packing halves the raw count for free, but does not change what survives that reduction (see
    `docs/explorations/polynomial_closure.txt`).

- **`embedder="bursar_graded"`**: graded packing (Section 3.4), which is Level 0 and subsumes involution packing.
  - It uses a dihedral subgroup `<h, s>` of `G`, where `s h s = h^-1`. Each packed coordinate placed at power
    `x^j` is multiplied by `exp(2 pi i j / r)` under `h`. So `B(zeta x, y) = B(x, y)`, and every coefficient
    `E_{c,j}` whose `j` is not a multiple of `r` vanishes identically. `s` acts as complex conjugation, so the
    remaining coefficients are real.
  - `G`-invariant coordinates can be output directly. If `h` then multiplies every remaining coordinate by one
    scalar, the grading is by `c` instead (scalar mode).
  - The `(h, s)` with the fewest outputs is found automatically, trying one `h` per conjugacy class.
  - Coefficients whose power of `x` is not a sum of `c` occupied positions (possible when some positions are
    empty padding) are also omitted, since they are zero for every input.
  - Typical sizes: 11 → 4 for `S(2)`, `k = 2` (your zip); 174 → 32 for `S(3)`, `k = 2`; 138060 → 13142 for
    `S(5)`, `k = 4`.
  - With `reductions="safe"`, the `c = 1` outputs are reduced exactly (Level 1): each is `(p_j S) w`, with
    integer coefficients, where `S` is the sum of the group matrices. Any that is a rational combination of the
    invariant outputs or of earlier ones is dropped. For example `S(5)`, `k = 4` goes from 13142 to 13140.
  - For `c >= 2` the coefficients involve `exp(2 pi i / r)`, so the exact Level 1 search does not yet reach them.
    `meta` shows the complex packed coordinates with numerical coefficients. Their exact definition is a projection
    of slot functionals onto weight spaces; the weight-0 coordinates and the invariants are exact integer rows.
  - `bursar_involution` is kept as the simpler special case (`h` the identity).

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

# Involution packing (about half the outputs, still Level 0), optionally with the safe reductions:
embed2, meta2 = build_invariant_embedding(2, spec, EmbeddingOptions(embedder="bursar_involution", reductions="safe"))
meta2.packed_coordinates()   # [(real part, imaginary part), ...] of each packed coordinate, in the d's and e's
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

`python -m symembed.explore.polynomial_closure` writes `docs/explorations/polynomial_closure.txt`. For each case and
for both embedders it keeps only the outputs that are not exact polynomials in those already kept, both formally in
`w` and as functions of `L` (where the Gram and Plücker relations hold). Every identity is verified exactly.

## Layout

| file | contents |
|---|---|
| `spec.py` | `GroupSpec`, `EmbeddingOptions` |
| `slots.py` | slots of `w`, orbit grouping, group elements, signed-permutation table (pure Python, independent of symatom) |
| `bursar.py` | Bursarial coefficients from `{g·w}` in torch, real and complex |
| `involution.py` | involution packing: generic packing for any linear involution, slot-form packing, choice of `g0` |
| `graded.py` | graded packing: weight-space packing for a dihedral `<h, s>`, search, evaluator, metadata |
| `reductions.py` | exact (Gaussian-integer) expansion of low-`c` outputs, exact linear-dependency search (Level 1) |
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
- involution packing: that `g0` acts as conjugation on the packed coordinates, for every involution and for a generic
  non-permutation involution; that the imaginary parts vanish, symbolically and numerically; the output counts;
  invariance; and the safe reductions applied on top;
- graded packing: that every candidate packing grades correctly, is real and is injective; that the omitted
  coefficients vanish; the output counts; invariance; `gradcheck`; that `h` = identity reproduces involution packing;
- comparison against full symbolic expansion;
- a cross-check against the legacy Bursarial embedder;
- invariance under permutations, rotations, and reflections (when `P` is in `G`);
- shapes and dtypes;
- `gradcheck`;
- the metadata.
