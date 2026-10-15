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
  - It uses a dihedral subgroup `<h, s>` of `G`, where `s h s = h^-1`. A packed coordinate of weight `k` (multiplied
    by `zeta^k`, `zeta = exp(2 pi i / r)`, under `h`) is placed at a power `x^j` with `j = k + shift (mod r)`. Then
    `B(x, y)` and `B(zeta x, zeta^shift y)` agree up to a constant, and every coefficient `E_{c,j}` with
    `j != shift * c (mod r)` vanishes identically. `s` acts as complex conjugation, so the remaining coefficients
    are real.
  - The shift, and the choice of weight `k` or `r - k` for each conjugate pair of coordinates, minimise the length
    `W` of `q` (the number of coefficients is about `sum_c c (W - 1) / r`). For `S(3) x S(2)`, `k = 4` this takes
    `W` from 10 to 8 and the output count from 124 to 100.
  - `G`-invariant coordinates can be output directly. If `h` then multiplies every remaining coordinate by one
    scalar, the grading is by `c` instead (scalar mode).
  - The `(h, s)` with the fewest outputs is found automatically, trying one `h` per conjugacy class.
  - Coefficients whose power of `x` is not a sum of `c` occupied positions (possible when some positions are
    empty padding) are also omitted, since they are zero for every input.
  - So are coefficients killed by a sign (character) symmetry of individual coordinates (Section 3.4, "Exact zeros
    from symmetries of single coordinates"). For example, when a group element multiplies the coordinate at
    position 0 by -1, then `E_{c,0}` vanishes for every odd `c` (`S(3)`, `k = 3`: 20 → 17 outputs). Among the
    shortest layouts, the coordinates placed at the two ends are chosen to maximise these zeros (`S(5)`, `k = 4`:
    60 of them).
  - Typical sizes: 11 → 4 for `S(2)`, `k = 2` (your zip); 174 → 24 for `S(3)`, `k = 2`; 1494 → 100 for
    `S(3) x S(2)`, `k = 4`; 138060 → 13082 for `S(5)`, `k = 4`.
  - With `reductions="safe"`, the `c = 1` outputs are reduced exactly (Level 1): each is `(p_j S) w`, with
    integer coefficients, where `S` is the sum of the group matrices. Any that is a rational combination of the
    invariant outputs or of earlier ones is dropped. For example `S(5)`, `k = 4` goes from 13082 to 13081.
  - The packing is built exactly in the cyclotomic integers `Z[zeta_L]`, with `L = lcm(r, 4)`. These are handled
    by `cyclotomic.py`, with each number stored as a short integer vector. The floating-point coordinates are derived
    from the exact ones, and grading and reality are checked exactly at build time.
  - For `2 <= c <= c_max` the outputs are also expanded exactly. The build verifies exactly that the omitted
    coefficients vanish and the kept ones are real. It then drops outputs that are linear combinations of kept ones
    with coefficients in `Q(zeta_L)`, found by "restriction of scalars": each output `P` contributes the integer
    rows `P, zeta P, ...`, so the integer relation search applies unchanged. Gains are small, for example
    `S(4)`, `k = 2` goes from 613 to 612.
  - Beyond `c_max` (unless `linear_reduction_max_c` is set), outputs that vanish identically are dropped too. They
    are proposed by exact evaluation modulo two primes at random points, and each is then verified by an exact
    expansion of that one coefficient, restricted to the positions that can contribute to it (cheap near the ends
    of the range of `j`). For example `S(3) x S(2) x P`, `k = 4` goes from 100 to 99 (`Re E_{5,32} = 0`).
  - `meta` still shows the complex packed coordinates with numerical coefficients, although exact ones are now
    available in `meta.packing.exact`.
  - `bursar_involution` is kept as the simpler special case (`h` the identity).

- **`embedder="bursar_nested"`**: nested coset packing (Section 3.5), which is Level 0.
  - `MS` is the disjoint union, over right cosets `Hg`, of the `H`-orbits of `g.w`. Each `H`-orbit is embedded
    by graded packing for `H`, and the multiset of those vectors by a plain Bursarial embedding.
  - It is roughly 2x smaller than flat graded packing when there are many cosets: for `S(5)`, `k = 4`, 13082 →
    5455, with `H = S(2) x S(3)` and 10 cosets.
  - It is used only when it wins; otherwise the graded embedding is returned. The exact Level 1 searches are not
    applied to it.

- **`embedder="bursar_cascade"`**: cascaded graded embedding (Section 3.6), which is Level 0.
  - For a dihedral (or cyclic) subgroup `D = <h, s>`, the packed coordinates `z(w)` determine the whole `D`-orbit
    (roots of unity and conjugation), so no group element is ever applied to `w`. The orbit is embedded in stages
    along a chain of cyclic subgroups of `<h>`: each stage is the graded Bursarial polynomial of `p` points (a
    small prime or factor `p`), with positions `= weight + shift (mod the remaining order)` so that its outputs are
    again eigen-coordinates; a final stage embeds `{Y, conj Y}` by `Re Y` and `Q conj(Q)`.
  - When `D` is smaller than `G`, the `D`-orbits over right cosets are nested under a plain Bursarial embedding.
  - Sizes: `S(3) x S(2)`, `k = 4`: 63 (graded 100; hand-crafted 54); with `P`: 52; `S(4)`, `k = 4`: 153;
    `S(5)`, `k = 4`: 3420 (nested graded 5455).

- `EmbeddingOptions(root_normalise=True)` replaces each output that is provably homogeneous of degree `r` in `L` by
  `sign(f)|f|^(1/r)`, which scales linearly with `L` and tames the dynamic range for float32. Outputs that are
  numerically zero at random inputs (typically identically zero through the Gram relations) are left as they are,
  because a root would amplify their rounding noise. It is not available for `bursar_nested`, and it combines with
  `parity_split`.

`homogenise`, the complex field and Level 2 reductions are not implemented.

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
- `EmbeddingOptions(parity_split=True)` gives every output definite parity under reflections. Outputs that an
  exact counting argument proves pure are kept; every other output `o` becomes `(o(L) + o(pL))/2` and
  `(o(L) - o(pL))/2`. For graded embeddings the argument also uses coordinates on which parity acts as complex
  conjugation (`(d's) + i (e's)`), so for `S(3) x S(2)`, `k = 4` exactly the mixed outputs are split (100 → 184).
  It does nothing when parity is in `G`, and is not available for `bursar_nested`.
- `EmbeddingOptions(metric=(1, -1, -1, -1))` uses Minkowski dot products. Any diagonal metric of `+1`/`-1` entries
  works; the default is Euclidean. Only the dot products change. Every construction and count downstream of `w` is
  identical, and injectivity on the multiset is unaffected.
  - For indefinite metrics some orbits are not closed, for example that of a null vector, whose closure contains 0.
    So the step back from the multiset to `L` holds only up to orbit closures. No continuous map can do better, and
    continuous invariant functions need no more. See Remark "Indefinite metrics" in Section 1 of the notes.

### Collider events: the beam frame (`build_beam_embedding`, Section 6 of the notes)

```python
from symembed import BeamSpec, EmbeddingOptions, build_beam_embedding
from symembed.beam import beam_inputs_from_lab, beam_inputs_from_beams, beam_inputs_from_pt_eta_phi_m

spec = BeamSpec(species=(("j1", "j2", "j3"),),     # three exchangeable visible particles
                masses=(None,),                    # per species: None (masses vary) or a fixed common mass
                parity=False, met=False, mass_variable="m")
embed, meta = build_beam_embedding(spec, EmbeddingOptions(embedder="bursar_cascade"))   # 47 outputs
Y = embed(beam_inputs_from_lab(four_momenta))      # (..., 3, 4) lab (E, px, py, pz) -> (..., meta.D)
```

- The beams `p, q` fix a frame: `w = (m_i, pt_i, transverse dots | pz_i, transverse wedges)`. What is left of the
  Lorentz group is rotations about the beam axis (removed by `w`) and the swap of the beams, which negates exactly
  the bracket (`pz`'s and wedges); identical particles act by signed permutations. `sqrt(s)` and fixed masses are
  parameters, not outputs.
- `parity=True` makes the embedding blind to handedness; with `reductions="safe"` the wedges are then omitted (the
  transverse Gram matrix fixes the transverse momenta up to `O(2)`). `met=True` adds a transverse-only object.
- `swap="group"` (default) treats the swap as a group element, so every embedder can use it; `swap="projective"`
  stores the bracket up to sign with Echinocoder's RP embedder, either before a plain Bursarial embedding
  (`embedder="bursar"`, a baseline) or as the first stage of a cascade (`embedder="bursar_cascade"`).
- Input: beam-frame `(m, px, py, pz)` per object, shape `(..., n_objects, 4)`. Front ends convert lab four-momenta,
  `(pt, eta, phi, m)`, or four-momenta in any frame together with the beams (crossing angles; also returns
  `sqrt(s)`).

## Exploration

`python -m symembed.explore.simple_cases` (run from the repo root, with the venv's Python) writes very simple cases
out in full into `docs/explorations/`, for inspection by eye:

- `simple_cases.pdf` shows the slots of `w`, the action of `G`, and every output expanded in the `d`'s and `e`'s,
  marked kept or dropped, with the identity used to drop it. Expansions longer than 100 terms show only their term
  count.
- `simple_cases_full.txt` has every expansion, with no limit.

The cases are S(2) on (a, b), S(3) on (a, b, c) and S(2) x S(1) on (a, b, p), each with k = 2 and with k = 1.

`python -m symembed.explore.simple_cases_graded` writes the graded embeddings of the same cases out in full
(`docs/explorations/simple_cases_graded.pdf` and `simple_cases_graded_full.txt`). Each shows the chosen `(h, s)`,
the invariants, the packed coordinates and every output, all with exact coefficients such as `sqrt(3)/2`.

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
| `graded.py` | graded packing: weight-space packing for a dihedral `<h, s>`, search, evaluator, metadata, exact reductions |
| `nested.py` | nested coset packing: subgroup and coset search, evaluator, metadata |
| `cascade.py` | cascaded graded embedding: stage planning along chains of cyclic subgroups, search, evaluator, metadata |
| `beam.py` | beam-frame front end for collider events: `BeamSpec`, slots and signed permutations, front ends, builder |
| `parity.py` | option `parity_split`: exact purity check, even/odd parts |
| `rootnorm.py` | option `root_normalise`: exact homogeneity check, signed roots with finite gradient at 0 |
| `cyclotomic.py` | exact arithmetic in `Z[zeta_L]` (integer vectors in the power basis) |
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
