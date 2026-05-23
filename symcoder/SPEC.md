# symcoder — Developer Specification

symcoder is the numerical encoding layer above symatom.  It takes a `Plan`
(context + operations) and a physics event (dict of label→vector) and produces
a permutation-invariant float64 vector via a two-phase pluggable encoder
hierarchy.

---

## Dependencies

- **symatom** — provides `Plan`, `Context`, `VectorType`, `Operation`,
  `ArgumentSymmetry`, `Atom`, `repS`, `canonical_pair_flavours`,
  `FlavouredOperator`, `PairFlavour`.
- **numpy** — all numerical output is `np.ndarray` of dtype `float64`.
- **Echinocoder polynomial embedder** — loaded at import from
  `Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py` in the
  repo root.  Not a package; located relative to `symcoder/`.

---

## Public API

### `EvaluableOperation`

Subclass of `symatom.Operation` (frozen dataclass).  Adds an `eval_fn` field.

```python
EvaluableOperation(
    name:              str,
    rank:              int,
    parity:            int,          # +1 or -1
    argument_symmetry: ArgumentSymmetry,
    eval_fn:           Callable,     # ([array, ...]) -> float
)
```

**Contracts:**
- `eval_fn(vectors)` receives a list of arrays in `atom.labels` order and
  returns a plain Python `float`.  It must return the raw value *before* the
  atom's sign is applied; `evaluate()` applies the sign.
- `eval_fn` is excluded from equality and hashing.  Two `EvaluableOperation`
  objects with identical `(name, rank, parity, argument_symmetry)` compare
  equal regardless of `eval_fn`.
- For `ANTISYMMETRIC` operations, `eval_fn` must itself be antisymmetric under
  argument transposition.

---

### `evaluate(atom, event) -> float`

```python
evaluate(atom: Atom, event: dict[str, array]) -> float
```

Returns `atom.sign * atom.operation.eval_fn([event[lbl] for lbl in atom.labels])`.

Raises `TypeError` if `atom.operation` is not an `EvaluableOperation`.

---

### `encode(plan, event, orbit_factory, phase2_factory) -> np.ndarray`

```python
encode(
    plan:           Plan,
    event:          dict[str, array],
    orbit_factory:  OrbitEncoderFactory | None = None,
    phase2_factory: Phase2EncoderFactory | None = None,
) -> np.ndarray   # dtype float64
```

Main encoder.  Returns a 1-D float64 array.  The length is identical to
`sum(s.length for s in describe_encoding(plan, orbit_factory, phase2_factory))`.

Passing `orbit_factory=None` skips Phase 1; passing `phase2_factory=None` skips
Phase 2.

---

### `encode_and_describe(plan, event, orbit_factory, phase2_factory) -> (np.ndarray, EncodingTree)`

```python
encode_and_describe(
    plan:           Plan,
    event:          dict[str, array] | None,
    orbit_factory:  OrbitEncoderFactory | None = None,
    phase2_factory: Phase2EncoderFactory | None = None,
) -> tuple[np.ndarray, EncodingTree]
```

Single code path for both encoding and describing.  Passing `event=None` runs
describe-only mode: values are zeros, tree is identical to `describe_encoding()`.

---

### `describe_encoding(plan, orbit_factory, phase2_factory) -> EncodingTree`

```python
describe_encoding(
    plan:           Plan,
    orbit_factory:  OrbitEncoderFactory | None = None,
    phase2_factory: Phase2EncoderFactory | None = None,
) -> EncodingTree
```

Pure function (no event data).  Returns an `EncodingTree` whose segments mirror
the flat output of `encode()`.  Delegates to `encode_and_describe(event=None)`.

---

### `decode_alignment(plan, phase1_results, all_pair_decoded, atol) -> AnnotatedMultisetOfRepSEvalVectors`

```python
decode_alignment(
    plan:             Plan,
    phase1_results:   dict[(str, tuple), AnnotatedMultisetOfReals],
    all_pair_decoded: list[AnnotatedMultisetOfRealPairs],
    atol:             float = 0.0,
) -> AnnotatedMultisetOfRepSEvalVectors
```

Alignment decoder.  Recovers the G-orbit of `repS` evaluations from Phase 1
and Phase 2 decoded outputs.  See `alignment_decoder.py` for the full algorithm.

**Note:** Requires Phase 2 encoding with `use_complementarity_drop=False` for
plans containing antisymmetric cross-FO blocks (see `test_alignment_decoder.py`
for explanation).

---

## Encoder Factories

### `OrbitEncoderFactory` — Phase 1

```python
OrbitEncoderFactory(factories: list[AtomOrbitEncoderFactory])
```

Top-level Phase 1 factory.  For each `FlavouredOperator` in
`repS(plan.context, plan.operations)`, calls each sub-factory's `assess()`,
collects all returned encoders, and selects the one with smallest `output_dim`
(fewest stored reals).  Ties broken by `priority` (higher wins).

**Built-in sub-factories:**

| Factory | Output dim | Notes |
|---|---|---|
| `HalfSortEncoderFactory` | `orbit_size / 2` | Physical `(atom, −atom)` partition; priority 1.0 |
| `SortEncoderFactory` | `orbit_size` | Universal fallback; priority 0.5 |

`PolyEncoderFactory` exists as a stub but raises `NotImplementedError`; omit it
from the factory list.

**Typical instantiation:**

```python
orbit_factory = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
```

---

### `Phase2EncoderFactory` — Phase 2

```python
Phase2EncoderFactory(factories: list[OverlapBlockEncoderFactory])
```

Top-level Phase 2 factory.  Delegates to its `OverlapBlockEncoderFactory`
list; the first factory that builds a non-empty encoder wins.

---

### `OverlapBlockEncoderFactory`

```python
OverlapBlockEncoderFactory(
    row_pair_factories: list[PairOrbitEncoderFactory],
    use_complementarity_drop: bool = True,
)
```

Groups `PairFlavour`s into overlap blocks (fixed `op_u, flavour_u, op_v,
flavour_v`).  Within each block, applies complementarity drop if enabled
(removes the largest non-self-pair PairFlavour; deducible from the others by
multiset complementation).  For each remaining PairFlavour, calls each
row-pair factory's `assess()` and selects the encoder with smallest
`output_dim`.

**`standard_row_pair_factories()`** returns the default ordered list:

| Factory | Encoder | Output dim | Orbit type | Notes |
|---|---|---|---|---|
| `SelfPairEncoderFactory` | `NullPairEncoder` | 0 | NULL_SELF | Self-pairing, full overlap |
| `SelfPairNegEncoderFactory` | `SelfPairNegEncoder` | n | TYPE_NEG | σ-compressed self-pairing |
| `SelfPairType11EncoderFactory` | `SelfPairType11Encoder` | n | TYPE_11 | σ-compressed self-pairing |
| `NegPairEncoderFactory` | `NegPairEncoder` | 2n | TYPE_NEG | Non-self-pairing |
| `Type22PairEncoderFactory` | `Type22PairEncoder` | 2n | TYPE_22 | Both signs flip |
| `Type21PairEncoderFactory` | `Type21PairEncoder` | 2n | TYPE_21 | u-sign flips |
| `Type12PairEncoderFactory` | `Type12PairEncoder` | 2n | TYPE_12 | v-sign flips |
| `Type11PairEncoderFactory` | `Type11PairEncoder` | 2n | TYPE_11 | Universal fallback |

Here n = count of positive-sign representatives in the pair orbit.
`min(output_dim)` selection ensures σ-compressed encoders are preferred
automatically when applicable.

**Typical instantiation:**

```python
phase2_factory = Phase2EncoderFactory([
    OverlapBlockEncoderFactory(standard_row_pair_factories())
])
```

---

## Orbit Types (TYPE_11/12/21/22/NEG)

Orbit types classify the achievable sign combinations in a pair orbit.  They
are STRUCTURAL properties of the orbit — determined by the group acting on the
atom labels — NOT by the `argument_symmetry` of the operations.

| Type | Achievable sign set | Output reals | σ-compression possible? |
|---|---|---|---|
| TYPE_11 | `{(+,+)}` | 2n (standard) / n (σ) | Yes for self-pairing blocks |
| TYPE_12 | `{(+,+),(+,−)}` | 2n | No |
| TYPE_21 | `{(+,+),(−,+)}` | 2n | No |
| TYPE_22 | `{(+,+),(+,−),(−,+),(−,−)}` | 2n | Not yet analysed |
| TYPE_NEG | `{(+,+),(−,−)}` | 2n (standard) / n (σ) | Yes for self-pairing blocks |

**σ-compression (self-pairing blocks):** For blocks where `op_u == op_v` and
`F_u == F_v`, the swap symmetry `(u_k, v_j) ↔ (v_j, u_k)` makes the z-value
multiset σ-closed (`z → i·conj(z)` for TYPE_11; `{z²}` τ-closed for TYPE_NEG),
halving storage from 2n to n reals.  See `DOCS/sigma_symmetry_analysis.md`.

**Never use SS/SA/AS/AA** terminology for orbit types — those describe
ArgumentSymmetry of the two operations and do not determine the orbit type.

---

## Output Types for Decoding

```python
AnnotatedMultisetOfReals         # Phase 1 orbit decoder output
AnnotatedMultisetOfRealPairs     # Phase 2 row-pair decoder output
AnnotatedMultisetOfRepSEvalVectors  # alignment decoder output
```

Each stores a list treated as a multiset (order not guaranteed) together with
the algebraic atoms that produced the values.

---

## EncodingTree

Hierarchical descriptor of the encoded output.

```python
EncodingTree
    .phase1: Phase1Tree
        .orbits: list[SegmentInfo]     # one per FO orbit
    .phase2: Phase2Tree
        .blocks: list[OverlapBlockNode]
            .segments: list[SegmentInfo]  # one per PairFlavour (incl. NULLs)
```

Supports flat iteration: `for s in tree`, `tree[i]`, `len(tree)` all delegate to
`tree.flat()` which yields all `SegmentInfo` leaves in the same order as the
encoded output.

---

### `SegmentInfo`

Frozen dataclass.  Key fields:

| Field | Type | ORBIT | ASSOC | NULL |
|---|---|---|---|---|
| `kind` | `str` | `"ORBIT"` | `"ASSOC"` | `"NULL_SELF"` or `"NULL_COMP"` |
| `start` | `int` | ✓ | ✓ | ✓ |
| `length` | `int` | encoder-dependent | encoder-dependent | `0` |
| `op_u` | `str` | ✓ | ✓ | ✓ |
| `flavour_u` | `tuple[int,...]` | ✓ | ✓ | ✓ |
| `op_v` | `str\|None` | `None` | ✓ | ✓ |
| `flavour_v` | `tuple[int,...]\|None` | `None` | ✓ | ✓ |
| `overlap` | `tuple[int,...]\|None` | `None` | ✓ | ✓ |
| `symmetry_class` | `str\|None` | `None` | `"11"/"12"/"21"/"22"/"neg"/"11_sigma"/"neg_sigma"` | same |
| `sign_compressed` | `bool\|None` | `True`/`False` | `None` | `None` |
| `method_name` | `str\|None` | `"sort"/"half_sort"` | `"11"/"12"/…` | `"null_self"` |
| `notional_length` | `int\|None` | full orbit size | uncompressed pair count | dropped size |

**Derived property:** `stop = start + length`.

**Methods:**
- `__str__()` — fixed 9-token format (pipe through `column -t` for alignment).
- `to_dict()` — JSON-serialisable dict.

---

## Internal Modules

| Module | Responsibility |
|---|---|
| `eval.py` | `EvaluableOperation`, `evaluate` |
| `encode.py` | `encode`, `encode_and_describe`, `describe_encoding` |
| `describe.py` | `SegmentInfo`, `EncodingTree`, `Phase1Tree`, `Phase2Tree`, `OverlapBlockNode` |
| `alignment_decoder.py` | `decode_alignment` |
| `decoded_types.py` | `AnnotatedMultisetOfReals`, `AnnotatedMultisetOfRealPairs`, `AnnotatedMultisetOfRepSEvalVectors` |
| `pairs.py` | `eval_pair_orbit` (brute-force reference), `_is_self_pair` (internal) |
| `sign_correlation.py` | `_complex_to_reals` (internal) |
| `encoders/` | All encoder/factory classes (see sub-modules below) |

**`encoders/` sub-modules:**

| Module | Contents |
|---|---|
| `_base.py` | `AtomOrbitEncoder`, `AtomOrbitEncoderFactory`, `PairOrbitEncoder`, `PairOrbitEncoderFactory`, `OrbitSpec`, `PairOrbitSpec`, `EncodingResult` |
| `orbit_encoder.py` | `OrbitEncoder`, `OrbitEncoderFactory` |
| `sort_encoder.py` | `SortEncoder`, `SortEncoderFactory`, `HalfSortEncoder`, `HalfSortEncoderFactory` |
| `poly_encoder.py` | `PolyEncoder`, `PolyEncoderFactory` (stub — `NotImplementedError`) |
| `row_pair_encoders.py` | All row-pair encoder/factory classes; `standard_row_pair_factories()` |
| `overlap_block.py` | `OverlapBlockEncoder`, `OverlapBlockEncoderFactory` |
| `phase2_encoder.py` | `Phase2Encoder`, `Phase2EncoderFactory` |

---

## Key Invariants for Tests

1. `len(encode(plan, event, of, pf)) == sum(s.length for s in describe_encoding(plan, of, pf))`
   for all valid plans, events, and factory pairs.
2. `encode(plan, event, ...)` is identical for any permutation of labels within a group.
3. `encode(plan, event, ...)` returns an empty array when all groups are empty or no
   operations are given.
4. `describe_encoding(plan, of, pf)` is a pure function: identical trees for any
   two calls with the same plan and factories.
5. `NULL_SELF` and `NULL_COMP` segments have `length == 0` and `start == stop`.
6. Non-null segments tile the output array contiguously without gaps:
   `non_null_segs[i].stop == non_null_segs[i+1].start` for consecutive non-null segs.
7. `ORBIT` segments (`tree.phase1`) all precede `ASSOC`/`NULL` segments (`tree.phase2`).
8. Within each overlap block, at most one segment has `kind="NULL_COMP"`;
   zero or more have `kind="NULL_SELF"`.
9. For TYPE_11 / TYPE_NEG self-pairing blocks, `method_name` is `"11_sigma"` or
   `"neg_sigma"` (σ-compressed), and `output_dim == len(representatives)`.
10. `decode_alignment` requires Phase 2 encoding with `use_complementarity_drop=False`
    for plans with antisymmetric cross-FO blocks.
