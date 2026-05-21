# symcoder — Developer Specification

symcoder is the numerical encoding layer above symatom.  It takes a `Plan`
(context + operations) and a physics event (dict of label→vector) and produces
a permutation-invariant float64 vector.  The methodology is documented in full
in `docs/encoder.tex`.

---

## Dependencies

- **symatom** — provides `Plan`, `Context`, `VectorType`, `Operation`,
  `ArgumentSymmetry`, `Atom`, `repL`, `canonical_pair_flavours`,
  `FlavouredOperator`, `PairFlavour`.
- **numpy** — all numerical output is `np.ndarray` of dtype `float64`.
- **Echinocoder polynomial embedder** — loaded at module import from
  `Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py` in the
  repo root.  Not a package; located relative to `symcoder/encode.py`.

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
  returns a plain Python `float` (or float-castable).  It must return the raw
  value *before* the atom's sign is applied; `evaluate()` applies the sign.
- `eval_fn` is excluded from equality and hashing.  Two `EvaluableOperation`
  objects with identical `(name, rank, parity, argument_symmetry)` compare
  equal regardless of `eval_fn`.
- For `ANTISYMMETRIC` operations, `eval_fn` must itself be antisymmetric under
  argument transposition (the sign machinery assumes this).

---

### `evaluate(atom, event) -> float`

```python
evaluate(atom: Atom, event: dict[str, array]) -> float
```

Returns `atom.sign * atom.operation.eval_fn([event[lbl] for lbl in atom.labels])`.

Raises `TypeError` if `atom.operation` is not an `EvaluableOperation`.

---

### `encode(plan, event) -> np.ndarray`

```python
encode(plan: Plan, event: dict[str, array]) -> np.ndarray  # dtype float64
```

Main encoder.  Returns a 1-D float64 array whose length equals
`sum(s.length for s in describe_encoding(plan))`.

**Two-phase structure:**

**Phase 1 — ORBIT segments** (one per distinct `FlavouredOperator` in
`repS(plan.context, plan.operations)`):
- **No sign-paired atoms** (all currently implemented SYMMETRIC and
  UNSTRUCTURED operations): store
  `sorted({eval(u, event) for u in fo.atoms_one_per_sign()})`.
  Length: `fo.count_of_atoms_one_per_sign()` reals.
- **Sign-paired atoms** (all currently implemented ANTISYMMETRIC operations,
  sign compression): store
  `sorted({abs(eval(u, event)) for u in fo.atoms_one_per_sign() if u.sign == +1})`.
  Length: `fo.count_of_atoms_one_per_sign() // 2` reals.
  The negative-sign twins carry no independent information.
- FlavouredOperators are deduplicated by `(name, rank, parity,
  argument_symmetry, flavour.counts)`; the first occurrence wins.

**Phase 2 — ASSOC segments** (one per non-dropped `PairFlavour`):
- PairFlavours from `canonical_pair_flavours` are grouped into **overlap
  blocks** (fixed `op_u, flavour_u, op_v, flavour_v`).
- Within each block the association with the largest `pf.count()` is dropped
  (deducible by complementation; see `docs/encoder.tex` §4.3).
- For each remaining association, evaluate `z_k = eval(u_k) + i·eval(v_k)` for
  all atom-pairs with `u.sign == v.sign == +1` (`pf.count()` values), then
  embed via `_embed_compressed` and unpack to `2 * pf.count()` reals.

**Invariants:**
- Output is identical for any two events related by a permutation of labels
  within any group.
- Output length depends only on `plan`, not on `event`.
- Empty plan (no operations, or all-empty groups) returns `np.array([], dtype=float64)`.

---

### `encode_brute(plan, event) -> np.ndarray`

Reference implementation.  Applies Phase-2 sign compression (5a) but does
**not** apply the Phase-2 overlap-block dropping (5b).  Used for cross-checking
only.  Not suitable for production use on large plans.

---

### `describe_encoding(plan) -> list[SegmentInfo]`

```python
describe_encoding(plan: Plan) -> list[SegmentInfo]
```

Pure function (no event data).  Returns one `SegmentInfo` per output segment,
in the same order as `encode()`.  The total `sum(s.length for s in result)`
equals `len(encode(plan, any_event))`.

**Segment ordering:** all `ORBIT` segments first, then `ASSOC`/`NULL` segments
grouped by overlap block in canonical sort order.

---

### `SegmentInfo`

Frozen dataclass.  Fields:

| Field | Type | ORBIT | ASSOC | NULL |
|---|---|---|---|---|
| `kind` | `str` | `"ORBIT"` | `"ASSOC"` | `"NULL_SELF"` or `"NULL_COMP"` |
| `start` | `int` | ✓ | ✓ | ✓ |
| `length` | `int` | `fo.count_of_atoms_one_per_sign()` or `fo.count_of_atoms_one_per_sign()//2` | `2*pf.count()` | `0` |
| `op_u` | `str` | ✓ | ✓ | ✓ |
| `flavour_u` | `tuple[int,...]` | ✓ | ✓ | ✓ |
| `op_v` | `str\|None` | `None` | ✓ | ✓ |
| `flavour_v` | `tuple[int,...]\|None` | `None` | ✓ | ✓ |
| `overlap` | `tuple[int,...]\|None` | `None` | ✓ | ✓ |
| `symmetry_class` | `str\|None` | `None` | `"SS"/"SA"/"AS"/"AA"` | same |
| `sign_compressed` | `bool\|None` | `True`/`False` | `None` | `None` |

**Derived property:** `stop = start + length`.

**Methods:**
- `__str__()` — fixed 9-token format (pipe through `column -t` for alignment):
  `[start:stop]  kind  op_u  op_v  variant  u=(...)  v=(...)  shared=(...)  len=N`.
  The `variant` column carries the symmetry class (`SS`/`SA`/`AS`/`AA`) for
  pair rows, or `SC` (sign-compressed) / `.` for `ORBIT` rows.  Unused fields
  show `.`.  An optional `  |  example` tail follows when `example` is set.
- `to_dict()` — JSON-serialisable dict.  `ORBIT` segments include
  `sign_compressed`; `ASSOC`/`NULL` include `op_v`, `flavour_v`, `overlap`,
  `symmetry_class`.  All segments include `kind`, `start`, `stop`, `length`,
  `op_u`, `flavour_u`.

---

## Internal Modules

| Module | Responsibility |
|---|---|
| `eval.py` | `EvaluableOperation`, `evaluate` |
| `encode.py` | `encode`, `encode_brute`, embedding helpers |
| `describe.py` | `SegmentInfo`, `describe_encoding` |
| `pairs.py` | orbit/pair evaluation helpers |

**`pairs.py` functions** (not public, used by `encode.py`):

- `eval_single_orbit(fo, plan, event)` — all atom evals (both signs for ANTISYMMETRIC).
- `eval_single_orbit_compressed(fo, plan, event)` — Phase-1 values: `|eval|`
  of sign=+1 atoms for ANTISYMMETRIC; plain `eval` for SYMMETRIC.
- `eval_pair_orbit(pf, plan, event)` — full orbit z-values (reference).
- `eval_pair_orbit_positive(pf, plan, event)` — z-values for sign=(+1,+1) subset only.

**`encode.py` helpers** (not public):

- `_embed_compressed(z_pos, pf)` — embed `pf.count()` complex z-values into
  `pf.count()` complex polynomial coefficients, exploiting sign symmetry.
  Returns complex array; caller unpacks via `_complex_to_reals`.
- `_complex_to_reals(c)` — unpack n complex to 2n reals: `[re0,im0,re1,im1,…]`.
- `_sort_encode(values)` — sort a real array; Phase-1 primitive.
- `_overlap_block_key(pf)` — key for `itertools.groupby` to group PairFlavours
  into overlap blocks.  Must agree with `canonical_pair_flavours` sort order.
- `_encoding_canonical_key(pf)` — deduplication key for `encode_brute`.

---

## Key Invariants for Tests

1. `len(encode(plan, event)) == sum(s.length for s in describe_encoding(plan))`
   for all valid plans and events.
2. `encode(plan, event)` is identical for any permutation of labels within a group.
3. `encode(plan, {})` returns an empty array when all groups are empty.
4. `encode_brute` and `encode` produce the same total information (verified by
   permutation-invariance tests; exact equality not guaranteed due to different
   Phase-2 structure).
5. `describe_encoding` is a pure function: repeated calls return equal lists.
6. `ORBIT` segments precede all `ASSOC`/`NULL` segments.
7. Within each overlap block in Phase 2, at most one segment has `kind="NULL_COMP"`;
   zero or more segments have `kind="NULL_SELF"` (one per self-pairing association).
   Any block with at least one non-self-pairing association has exactly one `NULL_COMP`.
8. `NULL_SELF` and `NULL_COMP` segments have `length == 0` and `start == stop`.
9. Non-null segments (`ORBIT` and `ASSOC`) tile the output array contiguously without gaps.
10. `sign_compressed` is `True` iff the segment is `ORBIT` and its operation is
    `ANTISYMMETRIC`.
