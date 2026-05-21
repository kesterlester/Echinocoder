# NULL_COMP Decoder — Working Notes

*Written 2026-05-21 as pre-implementation research notes, before coding began.*

---

## Background: what NULL_COMP compression is

The original user called this "Step 4 compression" before we named it NULL_COMP.
The full derivation is in `DOCS/encoder.tex` §4.3 ("Dropping the Largest Remaining
Association").

An **overlap block** groups all PairFlavours that share `(op_u, flavour_u, op_v, flavour_v)`
and differ only in their overlap counts.  Within each block, after NULL_SELF entries are
removed, the PairFlavour with the **largest count** is not encoded — it is labelled
`NULL_COMP` in the metadata and contributes 0 reals to the output.

The reason this is lossless: the full Cartesian product of Phase 1 values partitions
exactly into the sub-multisets belonging to each association in the block, so knowing all
but one of them — plus Phase 1 — determines the remaining one by multiset complement.

---

## The mathematical guarantee

Let:
- `U` = Phase 1 decoded values for `(op_u, flavour_u)` — a multiset of |U| reals
- `V` = Phase 1 decoded values for `(op_v, flavour_v)` — a multiset of |V| reals
- `Z_full` = the **multiset Cartesian product** `{(u_i, v_j) : u_i ∈ U, v_j ∈ V}`,
  which has `|U| × |V|` elements (with repetitions if any u_i or v_j coincide)

Key fact from `encoder.tex` §4.3 (proved there):
> The associations of block B **partition** Z_full into disjoint sub-multisets Z_pf,
> one per PairFlavour in the block.

Therefore:
```
Z_NULL_COMP = Z_full  ∖  Z_NULL_SELF  ∖  ⊔_{ASSOC pf} Z_pf
```
where `∖` is multiset subtraction and `⊔` is multiset union.

This identity holds for ALL events, including degenerate ones where some eval values
coincide — because everything is done in multiset arithmetic, so repeated values appear
with the correct multiplicity.

---

## The decoding algorithm

### Inputs available at decode time

Inside `OverlapBlockEncoder.decode(values, phase1_results)`:

1. `self._selections` — the full list of PairSelections, including the comp-dropped one.
   - `sel.is_comp_drop == True` → the NULL_COMP entry
   - `sel.encoder.output_dim == 0` → NULL_SELF entries  
   - otherwise → ASSOC entries

2. `phase1_results` — dict keyed by `(op_name, tuple(flavour.counts))` → `AnnotatedMultisetOfReals`

3. Block identity: `pf0 = self._selections[0].pf` gives `op_u`, `flavour_u`, `op_v`, `flavour_v`
   (same for all selections in the block by definition of an overlap block)

### Step-by-step

```
U = phase1_results[(op_u.name, tuple(flavour_u.counts))].values   # list of floats
V = phase1_results[(op_v.name, tuple(flavour_v.counts))].values   # list of floats

Z_full = Counter()
for u in U:
    for v in V:
        Z_full[(u, v)] += 1

# Remove NULL_SELF contributions (already decoded as [(a,a) for a in shared Phase 1 values])
for null_self_decoded in [d for sel, d in zip(selections, decoded_list) if sel.encoder.output_dim == 0]:
    for pair in null_self_decoded.pairs:
        Z_full[pair] -= 1

# Remove all ASSOC contributions (already decoded and polished)
for assoc_decoded in [d for sel, d in zip(selections, decoded_list) if sel.encoder.output_dim > 0]:
    for pair in assoc_decoded.pairs:
        Z_full[pair] -= 1

# Remainder is NULL_COMP
null_comp_pairs = [pair for pair, count in Z_full.items() for _ in range(count)]
```

The `atom_pairs` for the NULL_COMP result can be obtained by calling
`plan.orbit_enumerator.orbit_elements(comp_sel.pf, plan.context)` on the comp-dropped
selection's pf — or left empty (since the output is a multiset with no claimed bijection).

---

## Two subtraction strategies: exact and fuzzy

### The key invariant: polishing makes Strategy 1 exact for *any* Phase 1 encoder

After polishing, the Phase 2 decoded u-values are not merely *close to* Phase 1 values —
they *are* the Phase 1 values, copied bit-for-bit from `phase1_results[key].values`.
Z_full is also built from those same Python float objects.  Counter lookup therefore
finds exact matches regardless of how noisy the Phase 1 encoder is, because both sides of
the comparison come from the same source.

The noise level of the Phase 1 encoder affects only how close Phase 1 values are to
ground truth — it does not affect the exactness of the match between polished Phase 2
pairs and Z_full.

This means **Strategy 1 is correct whenever polishing is enabled, even with a polynomial
(noisy) Phase 1 encoder.** The earlier note that "polynomial Phase 1 breaks exact
matching" was wrong.

### Strategy 1 — exact Counter subtraction (works whenever `polish_outputs=True`)

```python
# Z_full built from phase1_results values (the same objects polishing uses)
# After polishing, every decoded pair is an exact (Phase 1 value, Phase 1 value) tuple
# → Counter lookup is exact
Z_full[(u_i, v_j)] ...   # exact float equality always works
```

Prerequisite: polishing enabled (which is the default).  **Implement this first.**

### Strategy 2 — fuzzy multiset subtraction (needed when `polish_outputs=False`)

When a caller disables polishing (`polish_outputs=False`), Phase 2 decoded values retain
their raw root-finding noise (~1e-7 for large orbits).  These values do not appear in
Z_full, so Counter subtraction fails.  Strategy 2 handles this case.

**Why would anyone disable polishing?**  Primarily as a diagnostic: running without
polishing and then applying Strategy 2 reveals the raw Phase 2 decode error.  This lets
the library report something like "your encoding is good to one part in 10⁸, because we
can round-trip decode to one part in 10⁹."  It is an intentional quality-measurement
mode, not a degraded mode.

**The fuzzy algorithm:**

1. Build `Z_full_list` = list of all |U|×|V| pairs (with repetitions).
2. Build `known_list` = all ASSOC + NULL_SELF decoded pairs (unpolished).
3. Use `scipy.optimize.linear_sum_assignment` on a Chebyshev-distance cost matrix between
   `known_list` and `Z_full_list` to find the best bijective match.
4. The unmatched entries of `Z_full_list` are the NULL_COMP pairs.

Cost function: Chebyshev distance `max(|u_decoded - u_full|, |v_decoded - v_full|)`.
Expected distances: ~1e-7 for large Phase 2 orbits, ~1e-9 for small ones.

**Note on Phase 1 pluggability.**  A future polynomial Phase 1 encoder would:
- encode a single orbit as polynomial coefficients, with decode via root-finding;
- have ~1e-9 noise (real roots, smaller degree than Phase 2 polynomials);
- be **continuous and differentiable** w.r.t. the event (vs. sort-based which is only
  continuous — sort order changes create non-differentiable kinks at swap events).

Users who never decode (the majority) may prefer differentiability over sort precision.
Under such an encoder, *with* polishing Strategy 1 still works exactly (polishing
replaces Phase 2 values with Phase 1 values, whatever those are).  *Without* polishing,
Strategy 2 would be needed, and the residual error would reflect Phase 2 noise (~1e-7),
not Phase 1 noise (~1e-9).

**Implementation order**: implement Strategy 1 first (simpler, covers the polished case).
Add Strategy 2 later as the diagnostic/unpolished path.

### Testing Strategy 2

**Option A — add a polynomial Phase 1 encoder.**  A `PolyPhase1Encoder` that embeds a
sorted orbit as polynomial coefficients and decodes by root-finding.  This is a genuine
future feature (for differentiability), so building it serves both purposes.  However, it
requires Strategy 2 to be tested in the *unpolished* path, because with polishing, even a
noisy Phase 1 gives exact Counter subtraction.

**Option B — inject artificial noise into Phase 2 decoded outputs in the test framework**
(i.e. bypass polishing and inject noise into the raw decoded values).  This gives direct
control over the noise level and does not rely on any encoder's root-finding noise being
conveniently sized for the specific test event.  Conceptually odd but well-controlled.

Both options valid; Option B is the quicker practical choice since it requires no new
encoder code.  Option A becomes the natural regression test once a PolyPhase1Encoder
exists.

---

## Concrete example: Block (dot,(2,0)) × (dot,(2,0))

From `standard_row_pair_factories()` on the test fixture (3 electrons, 2 muons):

```
Block (dot,(2,0)) x (dot,(2,0)):
    [NULL_COMP] overlap=(1,0) dim=6 notional=24 n_reps=6   ← SelfPairNeg(σ) encoder
    [NULL_SELF] overlap=(2,0) dim=0 notional=6  n_reps=3
```

- `op_u = op_v = dot`, `flavour_u = flavour_v = (2,0)` (choose 2 electrons from 3)
- `U = V = {dot(a,b), dot(a,c), dot(b,c)}` — 3 Phase 1 values, call them d1, d2, d3
- `Z_full` = 3×3 = 9 pairs: all (d_i, d_j)
- `NULL_SELF` (overlap=(2,0), full overlap): diagonal pairs `{(d1,d1),(d2,d2),(d3,d3)}` — 3 pairs
- `NULL_COMP` (overlap=(1,0), 1 shared electron): off-diagonal pairs — should be 6 pairs
  - `{(d1,d2),(d2,d1),(d1,d3),(d3,d1),(d2,d3),(d3,d2)}`
- Verification: 3 (NULL_SELF) + 6 (NULL_COMP) = 9 = |Z_full| ✓

This block has **only one non-self association** (overlap=(1,0)), so that one is always
the NULL_COMP (it's the only candidate, as noted in `encoder.tex` §4.3 remark).  Nothing
is actually ASSOC-encoded in this block; the entire non-self content is reconstructed from
Phase 1.

---

## Concrete example: Block (dot,(1,1)) × (dot,(1,1)) — the multi-ASSOC case

```
Block (dot,(1,1)) x (dot,(1,1)):
    [NULL_COMP] overlap=(0,0) dim=12 notional=48 n_reps=12
    [ASSOC]     overlap=(0,1) dim=12 notional=48 n_reps=12
    [ASSOC]     overlap=(1,0) dim=6  notional=24 n_reps=6
    [NULL_SELF] overlap=(1,1) dim=0  notional=12 n_reps=6
```

- `U = V` = 6 values: {dot(a,p), dot(a,q), dot(b,p), dot(b,q), dot(c,p), dot(c,q)}
- `Z_full` = 36 pairs
- `NULL_SELF` (full overlap): 6 diagonal pairs
- `ASSOC (0,1)` (shared muon, different electrons): encoded in 12 reals → 12 decoded pairs
- `ASSOC (1,0)` (shared electron, different muons): encoded in 6 reals → decoded pairs
- `NULL_COMP (0,0)` (no shared labels): recovered by complement from the above

Decoding order matters: both ASSOC results must be decoded (and polished) first, then
both subtracted from Z_full, along with NULL_SELF, to recover NULL_COMP.

---

## API design note for `OverlapBlockEncoder.decode()`

**Current API**: returns one `AnnotatedMultisetOfRealPairs` per *non-comp-dropped*
selection, skipping the NULL_COMP entry entirely.

**Required for NULL_COMP decoding**: the returned list should include all selections
(including the NULL_COMP one), or at minimum the caller needs access to the NULL_COMP
decoded result.

**Proposed API change** (to discuss with user before implementing):

Option A — return all selections, including NULL_COMP:
```python
results: list[AnnotatedMultisetOfRealPairs]  # one per selection, in order
# comp-dropped selection's result is reconstructed by complement
```
This would change the length of the returned list (currently `len(active_sels)`, 
after change `len(self._selections)`).  The existing test would need updating.

Option B — return same length but include NULL_COMP reconstruction as a separate return value:
```python
return decoded_list, null_comp_decoded
```

Option C — add a `decode_null_comp` method called separately after `decode()`.

**My inclination**: Option A is cleanest — always return one result per selection,
with NULL_COMP reconstructed in place.  The test change is straightforward.

---

## Prerequisite checklist before implementing

**Strategy 1 (exact, sort-based Phase 1):**
- [x] Phase 1 decoders working (SortEncoder, HalfSortEncoder)
- [x] Phase 2 leaf decoders working (all 6 types including σ-variants)
- [x] Polishing implemented (prerequisite for exact Counter subtraction)
- [x] OverlapBlockEncoder.decode() working for ASSOC + NULL_SELF (use_complementarity_drop=False)
- [ ] NULL_COMP reconstruction via Counter subtraction in OverlapBlockEncoder.decode() (use_complementarity_drop=True)
- [ ] Test for NULL_COMP using use_complementarity_drop=True

**Strategy 2 (fuzzy, polynomial Phase 1 or noisy Phase 1):**
- [ ] Strategy 1 complete and tested (prerequisite)
- [ ] Either: PolyPhase1Encoder implemented OR noise-injection test harness available
- [ ] Fuzzy bipartite matching subtraction implemented
- [ ] Test verifying Strategy 2 recovers correctly under injected Phase 1 noise

---

## IEEE floating-point edge cases: are they a problem for Strategy 1?

**Signed zeros (+0.0 vs -0.0):** Not a problem.  Python (and numpy) define `+0.0 == -0.0`
and `hash(+0.0) == hash(-0.0)`, so they are treated as the same Counter key.  Even the
full chain `0.0 (Python) → numpy.float64 → negate → float()` produces a value that
compares equal and hashes equal to the original.

**numpy.float64 vs Python float as dict keys:** Not a problem.  Python's data model
requires equal objects to hash equally, and numpy honours this: `np.float64(x) == float(x)`
and `hash(np.float64(x)) == hash(float(x))` for all finite values.  A Counter built with
numpy.float64 keys is correctly looked up with Python float keys (or vice versa).

**The Type21/22 negation chain:** Irrelevant after polishing.  Negation is a pre-polish
operation; polishing discards those intermediate values and replaces them with values
taken directly from `phase1_results[key].values` — the same objects used to build Z_full.
No conversion chain operates on those values.

**NaN:** The one genuine concern, but not from an IEEE-signed-zero angle.  In CPython,
dict lookup checks `key is stored_key` before `key == stored_key`, so a NaN can be found
as a Counter key — but *only* if the identical Python object is used for both store and
lookup.  A NaN created by a second operation (same bit pattern, different object) would
not find the first entry.  However, NaN can only arise from user error in `eval_fn`
(division by zero, inf−inf, etc.), not from the encoding/decoding pipeline under normal
operation.  The appropriate response is a defensive assertion at the start of NULL_COMP
decoding:

```python
assert not any(np.isnan(u) or np.isnan(v) for u in U for v in V), \
    "NaN in Phase 1 decoded values — check eval_fn for divide-by-zero or similar"
```

**Overall verdict:** Strategy 1 (exact Counter) is safe for all normal floating-point
values produced by well-behaved `eval_fn` implementations.  No special handling of IEEE
edge cases is required beyond the NaN guard.

---

## Gotchas to watch for

1. **Counter subtraction must not go negative**: if the partition identity holds and 
   polishing is exact, this cannot happen.  Add an assertion `assert all(v >= 0 for v
   in Z_full.values())` after subtraction as a debugging aid.

2. **Blocks with no NULL_COMP**: some blocks have `all(not s.is_comp_drop for s in
   self._selections)` — e.g. if every association is NULL_SELF.  No reconstruction
   needed; just return decoded ASSOC + NULL_SELF results.

3. **Self-pairing blocks (op_u == op_v, flavour_u == flavour_v)**: U == V, Z_full is a
   square matrix of pairs.  The diagonal is NULL_SELF.  Off-diagonal is split between
   ASSOC(s) and NULL_COMP.  The σ-compressed encoders (SelfPairType11, SelfPairNegσ)
   appear here; their decoded pairs may include off-diagonal values from both the upper
   and lower triangle.

4. **Non-self-pairing blocks (op_u ≠ op_v or flavour_u ≠ flavour_v)**: U and V are
   different multisets.  No NULL_SELF entries.  Z_full = rectangular (|U| × |V|) pairs.
   One entry will be ASSOC or NULL_COMP (and with a single non-self association block, it
   is always NULL_COMP — so Z_full itself IS the NULL_COMP, directly from U × V).

5. **Single-association blocks**: the one non-self-pair association is ALWAYS NULL_COMP.
   No ASSOC data was encoded.  So `Z_full - NULL_SELF` directly gives NULL_COMP.  This
   is a common case (see the `(dot,(0,2)) × *` entries in the test fixture above).

6. **atom_pairs for NULL_COMP**: the atom-pairs are knowable (call the orbit enumerator
   on `comp_sel.pf`) but not required for the test (which only checks values).  Can be
   populated lazily.
