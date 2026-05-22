# Decoder Specification

*Status: draft — being written alongside implementation planning.*

---

## Purpose and scope

Every object in the encoder hierarchy that has an `encode()` method **may** have a
corresponding `decode()` method.  This document specifies the *principle* those decoders
must follow, before any implementation is written.

### Why decoders are optional, and why they exist at all

The final library's intended use case does not need decoders: only the encoded output is
ever consumed by external callers, and those callers never attempt to invert it.

Decoders exist solely for **validation**.  It is straightforward to show that an encoder's
output is invariant under the intended group action, but invariance alone does not prove
the encoder is correct: an encoder that always outputs the zero vector is trivially
invariant, yet catastrophically wrong.  A working decoder is the only way to demonstrate
that encoding is *lossless beyond the intended ambiguity* — i.e. that the only
information discarded is precisely the G × H equivalence.

Because decoders are optional, an encoder that lacks a decoder is still a valid
contribution to the library, provided its correctness can be established by other means.
It simply cannot participate in the round-trip decode-validation test suite.

This document covers **Step A** only: decoding back to *invariant space* (atom evaluation
values).  Step B — reconstructing actual vectors from invariant values, up to the action
of the continuous group H (O(3), SO(1,3), …) — is deferred and will be specified
separately once Step A is complete and tested.

---

## Core principle: federated, symmetric decoding

The encoder hierarchy is:

```
OrbitEncoder            (Phase 1 top)
  └─ AtomOrbitEncoder   (one per FlavouredOperator orbit: SortEncoder, HalfSortEncoder, …)

Phase2Encoder           (Phase 2 top)
  └─ OverlapBlockEncoder (one per overlap block)
       └─ PairOrbitEncoder (one per row-pair: Type11, Type12, …, SelfPairType11, …)
```

The decoding hierarchy mirrors it exactly.  Each encoder knows only about its own piece:
it does not look up or down the tree except at the single well-defined junction described
below (overlap block ↔ Phase 1).  This federated structure means each decoder can be
tested independently before composing into the full decoder.

---

## The decoded-output type: annotated multiset

The fundamental output of every `decode()` call at Step A is an **annotated multiset**:

> A list of values (reals, or real-pairs, or similar) that is *semantically* a multiset
> (i.e. order carries no information and must not be relied upon), together with a
> reference to the *known, fixed algebraic set of atoms* (or atom-pairs) that produced
> those values, where the bijection between values and atoms is **unknown** and is
> precisely the G-ambiguity that the encoder was designed to erase.

### Output type corresponds to the mathematical object, not the encoder variant

Because encoding is pluggable — the same orbit can be encoded by SortEncoder or
HalfSortEncoder, and the same orbit of pairs can be encoded by Type11PairEncoder or
SelfPairType11Encoder — the decoder's output type must be determined by *what was
encoded*, not *how it was encoded*.  Every decoder for a Phase 1 orbit returns the same
type; every decoder for a Phase 2 orbit-of-pairs returns the same type.  This means
decoders must always fully decompress: a HalfSortDecoder reconstructs the full orbit
(both positive and negative values), and a σ-compressed pair decoder reconstructs the
full orbit of pairs.  The encoder variant is an implementation detail that the decoder's
output does not expose.

### Concrete type names

Two concrete types are expected at the outset (additional types may be added later):

- `AnnotatedMultisetOfReals` — for all Phase 1 orbit decoders; holds a list of floats
  (the full orbit, fully decompressed) and a reference to the algebraic atom set.
- `AnnotatedMultisetOfRealPairs` — for all Phase 2 row-pair decoders; holds a list of
  (float, float) tuples (the full orbit of pairs, fully decompressed) and a reference
  to the algebraic atom-pair set.

### API

Keep the API minimal initially.  Direct access to the internal value list and atom
reference is sufficient; do not add convenience methods speculatively.  But do add
them if it becomes apparent that some mildy non-trivial / fancy list-operations are 
being repeated in many places. If such repeats start to appear, THEN they should be
consolidated into convenience methods to make agorithm intent more readable and code
clearer and with less duplucation.

### Why lists, not Python sets

Floating-point rounding can make genuinely distinct values compare as equal (or vice
versa), so a Python `set` or `frozenset` would silently merge or miscount elements.
Store values in an ordinary list.

### Enforcing multiset semantics

The decoder is free to populate the internal list in any order, including a random
shuffle immediately after construction.  This makes any caller that accidentally depends
on list order fail immediately or behave erratically, surfacing the bug without requiring
explicit checks.

### Role of SegmentInfo

The SegmentInfo tree (from the encoding step) serves as the decoder's guide: it records
how each segment was produced (encoder method, dimensions, symmetry class, etc.) and is
the natural place to also carry the algebraic atom or atom-pair information needed for
decoding.  Extending SegmentInfo to hold atoms (or enough metadata to reconstruct them)
is preferred over requiring the caller to re-run the encoder hierarchy to obtain that
information.  The exact extension is to be determined during implementation.

---

## Phase 1 decoding

### SortEncoder → SortDecoder

The encoder evaluates each atom in the orbit and sorts the results.  The output is
already a canonical representative of the orbit under permutation; the sort order itself
is the artefact of encoding and carries no physical meaning.

The decoder is trivial:
- Input: the slice of the encoded array corresponding to this orbit.
- Output: `AnnotatedMultisetOfReals` — the same values, now explicitly tagged as
  unordered — together with a reference to the atom set.

No numerical work is needed.  The values in the encoded array are the decoded values.

### HalfSortEncoder → HalfSortDecoder

The encoder stores only the positive-signed half of an ANTISYMMETRIC orbit (n values
instead of 2n).  The decoder fully decompresses: it reads the n stored values, appends
their negatives, and returns an `AnnotatedMultisetOfReals` of all 2n values — the same
type and same semantics as the SortDecoder output.  No numerical work is needed.

---

## Phase 2 decoding: row-pair level

### General PairOrbitEncoder → PairOrbitDecoder

The encoder forms complex numbers z_k from pairs of atom evaluations and encodes the
multiset {z_k} as polynomial coefficients, erasing all ordering information about the
pairs.

The decoder:
- Input: the slice of the encoded array corresponding to this row-pair segment.
- Numerical step: recover the multiset of complex values {z_k} from the polynomial
  coefficients (root-finding).
- Algebraic step: from {z_k}, recover the corresponding multiset of real pairs
  (u_value, v_value) by inverting the encoder's (u, v) → z map.  This inversion is
  specific to each encoder variant.
- Output: `AnnotatedMultisetOfRealPairs`, together with the atom-pair set.

Each concrete encoder variant (Type11, Type12, Type21, Type22, Neg, SelfPairType11,
SelfPairNeg, …) must document its (u, v) ↔ z map and provide the corresponding
inversion.  For the σ-variants, all information needed for inversion is available from
the SegmentInfo for that segment.

The σ-compressed variants (SelfPairType11, SelfPairNeg) store fewer real values than
their non-compressed counterparts, but the mathematical object encoded is the same: an
orbit of atom-pairs.  Their decoders must fully decompress — recovering the complete
multiset of (u_value, v_value) pairs — and return the same `AnnotatedMultisetOfRealPairs`
type as any other row-pair decoder.  The σ-compression is an encoder implementation
detail invisible to the decoded output.

The atom-pairs are known algebraically; the bijection between decoded pairs and
atom-pairs is unknown (G-ambiguity, as always).

### NULL_SELF entries

A NULL_SELF pair was dropped because its z-values are entirely determined by Phase 1:
z_k = (1+i)·a_k, where a_k is the Phase 1 evaluation of the shared atom.  The
decoder for a NULL_SELF entry is almost trivial: read the Phase 1 decoded values for the
relevant orbit, and if those values are [u1, u2, u3, …] recover the real-pairs as
[(u1,u1), (u2,u2), (u3,u3), …].  No root-finding or complex numbers are needed.

---

## Phase 2 decoding: overlap block level

### ASSOC and NULL_SELF pairs

The overlap block calls each row-pair decoder in turn.  ASSOC decoders are fully
self-contained (see above).  NULL_SELF decoders receive the Phase 1 decoded output for
the shared orbit from the overlap block.

### NULL_COMP pairs (complementarity reconstruction)

**Ownership.** The complementarity drop is a decision made by the OverlapBlockEncoder,
not by any individual row-pair encoder.  By choosing to drop one of its constituent
pairs, the OverlapBlock takes on a corresponding obligation: it must be able to prove
that the drop was harmless by demonstrating, at decode time, that the full input it
received can be recovered from only the data it *did* choose to encode.

The dropped row-pair encoder bears no special decode responsibility.  The NULL_COMP
segment has no data in the encoded array; its reconstruction is the OverlapBlock's own
business.

**What the NULL_COMP decoder test checks.**  The test is *not* "does the reconstructed
multiset match what the dropped encoder's `decode()` would have returned?"  That would
require running the dropped encoder's decode path, which defeats the point.

The test is: "given only the encoded output of the non-dropped pairs (plus Phase 1
decoded results), does the OverlapBlock decoder reconstruct the complete set of atom-pair
evaluations that were originally presented to the overlap block?"  The reconstructed
output is compared directly against ground-truth atom-pair evaluations (evaluate the
atoms on the original event), not against any intermediate decoder output.

This also means NULL_COMP decoder tests must use `use_complementarity_drop=True` (the
default): with the drop switched off, there is no NULL_COMP segment and nothing to test.
The leaf-level decoder tests (`Type22`, `NegPairEncoder`) that switch the drop off do so
purely to exercise those row-pair decoders in isolation; they cannot double as NULL_COMP
tests.

**Mechanism.**  The OverlapBlock, having decoded all ASSOC and NULL_SELF pairs and
holding the Phase 1 decoded output, recovers the dropped pair's multiset of
(u_value, v_value) pairs by reversing the complementarity relation used during encoding.
The exact formula is deferred pending confirmation from the encoding source code.

### Output polishing: replacing noisy Phase 2 values with Phase 1 values

**Why polishing is needed.**  Phase 1 decoders (SortEncoder, HalfSortEncoder) recover
their values by reading them directly from the sorted array — no root-finding, near-zero
numerical error.  Phase 2 decoders recover their values by polynomial root-finding, which
accumulates error that grows with polynomial degree.  For high-degree orbits (n > 6) the
root-finding error can reach ~1e-9.

Both the Phase 1 and Phase 2 decoded multisets contain evaluations of the *same*
underlying algebraic atoms (one for u, one for v).  The Phase 1 decoded u-values are
simply a more accurate copy of the u-component of the Phase 2 decoded pairs.  Polishing
replaces the noisier Phase 2 u-values with the cleaner Phase 1 u-values (and similarly for
v), using an optimal matching to handle the fact that the two multisets may be in different
orders.

**The polishing algorithm.**  Given decoded Phase 2 pairs `{(u_k', v_k')}` and Phase 1
decoded values `{u_j''}` and `{v_j''}`:

1. Solve for the permutation π* that minimises Σ(u_k' − u_{π(k)}'')² using
   `scipy.optimize.linear_sum_assignment` (the Hungarian algorithm).
2. Replace each u_k' with u_{π*(k)}''.
3. Solve independently for the permutation ρ* that minimises Σ(v_k' − v_{ρ(k)}'')².
4. Replace each v_k' with v_{ρ*(k)}''.
5. The polished output is `{(u_{π*(k)}'', v_{ρ*(k)}'')}`.

The u and v assignments are solved *independently*.  The resulting pairs may therefore
have π* ≠ ρ*, meaning that u- and v-components from different decoded roots end up in the
same output pair.  This is not a concern (see "No algebraic wrongness" below).

**Correctness property.**  Polishing is always better than not polishing.

In the zero-error limit the Phase 1 and Phase 2 values agree exactly, and the Hungarian
matching finds the identity permutation; polishing is a no-op.  When numerical errors are
nonzero, polishing replaces imprecise Phase 2 values with more precise Phase 1 values.
Even in the extreme case where errors are large enough that the optimal assignment is the
"wrong" one (does not correspond to the true mathematical bijection between roots and
atoms), the assigned Phase 1 values are still closer to the ground truth than the original
noisy Phase 2 values, because the Phase 1 values *are* the ground truth up to their own
(much smaller) error.  The cause of any residual inaccuracy is always numerical rounding,
not the polishing step.

**No algebraic wrongness.**  At no point in the decoder is any individual floating-point
value linked to a specific algebraic atom — the output is always a multiset.  The π* and
ρ* permutations are numerical bookkeeping, not algebraic statements.  There is no
algebraic notion of a "wrong" pairing.

The ASSOC encoding does preserve a numerical correlation between u_k and v_k (they come
from the same complex root z_k = u + iv).  Independent polishing of u and v may break
this specific numerical correlation, but the output is still a valid annotated multiset of
real pairs that faithfully represents the original orbit up to the G-ambiguity.  The
correlation loss is acceptable because the output type makes no claim about per-element
correspondences.

**The `polish_outputs` flag.**  Each leaf row-pair decoder's `decode()` method accepts a
`polish_outputs=True` keyword argument.  Set to False to receive raw polynomial-decoded
values (useful for debugging the root-finding step in isolation, and for tests that
specifically check unpolished output).  Default is True.

**Where polishing lives.**  Polishing is performed inside the *leaf row-pair decoders*
(not at the OverlapBlock or Phase2 level), because:
- The leaf decoder has direct semantic access to which Phase 1 orbit its u- and v-atoms
  belong to.
- Polishing is logically an enhancement of the leaf decoder's output.
- Placing it at the leaf level keeps OverlapBlock and Phase2 decoders simple.

The leaf decoder receives the relevant Phase 1 decoded multisets as optional parameters
(`u_phase1: AnnotatedMultisetOfReals | None`, `v_phase1: AnnotatedMultisetOfReals | None`)
when polishing is enabled.

**NULL_SELF exception.**  NULL_SELF pairs are reconstructed directly from Phase 1 decoded
values (no polynomial inversion at all), so they are already at Phase 1 accuracy.  No
further polishing is needed or meaningful for NULL_SELF entries.

---

### Interface between Phase 1 and Phase 2

The top-level decoder runs Phase 1 first, producing an `AnnotatedMultisetOfReals` for
each orbit.  These are then passed (dependency injection) to each OverlapBlockDecoder
that needs them.  The OverlapBlockDecoder passes the relevant Phase 1 results into each
leaf row-pair decoder (for NULL_SELF reconstruction and, when `polish_outputs=True`, for
output polishing).

This coupling is therefore deeper than a single top-level handshake: Phase 1 decoded
values flow all the way down to the leaf row-pair decoders.  The direction is strictly
top-down (Phase 1 → Phase 2 → OverlapBlock → leaf), and the Phase 1 data is always
already-decoded when passed down, so no circularity arises.

---

## Phase 2 decoding: top level (Phase2Encoder → Phase2Decoder)

The Phase 2 decoder iterates over all overlap block decoders, supplying each with the
relevant Phase 1 decoded results, and collects their outputs.  No additional logic is
required at this level.

---

## Full encoding tree decoding

The EncodingTree decoder:
1. Decodes Phase 1: calls each AtomOrbitDecoder, collects `AnnotatedMultisetOfReals`
   objects.
2. Decodes Phase 2: calls Phase2Decoder (which calls OverlapBlockDecoders), passing in
   the Phase 1 decoded results where needed.
3. Returns a decoded tree mirroring the EncodingTree structure, with each leaf replaced
   by its annotated multiset.

---

## Alignment decoding: recovering the orbit of repS evaluations

### What the alignment decoder produces

The alignment decoder is the **root of the decoding tree**.  Just as an OverlapBlockDecoder
assembles the outputs of its constituent PairOrbitDecoders into a coherent block result,
the alignment decoder assembles the outputs of all Phase 1 orbit decoders and all
OverlapBlock decoders — one annotated multiset per Phase 1 orbit and one per Phase 2
association, each giving G-ambiguous evaluations of a specific set of atoms or
atom-pairs — into a single richer object: the **G-orbit of the repS evaluation vector**.

Concretely: `repS` is the canonical and **information complete** ordered list of all atoms across all
FlavouredOperators.  For a given event E, each element g ∈ G produces a permuted (and
possibly sign-modified) version of that list, yielding a vector of evaluation values
`eval(g · repS, E)`.  The collection of all such vectors, for all g ∈ G, is the G-orbit
of `eval(repS, E)`.  This is the alignment decoder's output, represented as a
**multiset** of these vectors (column order has no meaning; the orbit is a multiset, not a
sequence).

This object contains all G-invariant information present in the encoding and no more.
In particular, it does **not** contain a bijection between atoms and values — no such
bijection can be recovered, because the encoder deliberately discards it.  The orbit is
the quotient of the full assignment by G, and that is the correct and complete answer.

### The table picture

Represent the G-orbit of `eval(repS, E)` as a table:

- **Rows** are indexed by positions in `repS`.  Each row is identified with one Phase 1
  orbit (or, for larger rank, a specific position within one FlavouredOperator's atoms).
  Every element in a given row is a value that the Phase 1 decoder for that orbit could
  have returned.
- **Columns** are the elements of the G-orbit — one per group element g (up to the
  stabiliser of `repS`).  Columns form a **multiset**; their order carries no meaning.
- The entry at (row r, column g) is `eval(g · repS_r, E)`, where `repS_r` is the r-th
  atom of `repS`.

The number of columns is `|G| / |stab(repS, E)|`.  For a generic event,
`stab(repS, E) = {e}` and there are exactly `|G|` columns.

### What the OverlapBlock outputs provide: letterbox views

Each OverlapBlock for association `(fo_u, fo_v, s)` provides a **2-row letterbox view**
of the table: a multiset of `(u_val, v_val)` pairs, corresponding to two specific rows
of the table seen simultaneously.

Because a single Phase 1 orbit may have fewer atoms than the number of table columns
(its stabiliser within G is non-trivial), the same `(u_val, v_val)` pair may appear in
the OverlapBlock output with multiplicity greater than 1.  The multiplicity factor for
a given atom-pair `(u, v)` is `|stab(u, v)| / |stab(repS, E)|`; for generic events this
is simply `|stab(u, v)|`, a fixed algebraic quantity that can be precomputed from the
atom-pair structure alone.  After dividing out this multiplicity, the OverlapBlock output
is a faithful 2-row projection of the table.

### Row information content and join order

The alignment difficulty for each row is entirely determined by the number of
**distinct** values it contributes to the table: a row with many repeated values tells
the algorithm little about how to distinguish columns from one another.  The general
principle is:

> **Add rows in decreasing order of distinct-value count.**  A row with `d` distinct
> values out of `|G|` table columns carries `log₂(d!)` bits of column-discriminating
> information.  Rows with the most distinct values are the best seeds and the best
> early extensions; rows with the fewest distinct values are the worst and should be
> deferred until the columns have already been largely distinguished by earlier rows.

This single principle handles the entire spectrum uniformly:

- **High-duplication rows** (few distinct values, many repeats) contribute little per
  step and are naturally added late, at which point most ambiguity has already been
  resolved and their contribution is correspondingly cheap.
- **Extreme case — constant rows** (Phase 1 orbit size 1, all `|G|` column entries
  identical): these contribute zero discriminating information and are trivially handled
  at any point in the join — their value is just copied into every partial column at
  whatever step they are reached.  They may optionally be handled outside the main join
  loop as a micro-optimisation, but this is not required: the general mechanism
  handles them correctly regardless.
- **Fully distinct rows** (all `|G|` entries different): these carry the maximum
  possible information and are the ideal seeds; they fully determine column identities
  in a single step with no branching.

The catastrophic failure mode — spending `|G|!` effort on a row that contributes nothing
— is avoided automatically by this ordering, without needing a special-case detector for
any particular level of duplication.

*Illustrative example.*  In the standard `{electrons=(a,b,c), muons=(p,q)}` setup,
`G = S_3 × S_2`, `|G| = 12`.  The atom `dot(p, q)` has Phase 1 orbit size 1: `S_3`
does not touch muon labels and `S_2` acts on `{p,q}` symmetrically, leaving
`dot(p, q)` invariant.  Its row contains the single value `dot(p, q, E)` repeated 12
times — zero distinct values, zero discriminating information.  Under the ordering rule
it is deferred to last (or placed outside the join entirely); either way the cost for
this row is O(|G|) copy operations.  The concern is not this row in isolation but
ensuring the same reasoning applies uniformly to rows with orbit size 2, 3, etc., which
the distinct-value ordering rule achieves without additional case logic.

### The alignment algorithm: multi-way relational join

Formally, the alignment problem is a **multi-way relational join**.  Each OverlapBlock
output (after multiplicity correction) is a binary relation `R_{ij}` on
`(row_i values × row_j values)`.  The target table T is the unique k-ary relation whose
binary projections are exactly `{R_{ij}}`.

**Uniqueness (lossless-join property).**  Any candidate table T′ consistent with all
`R_{ij}` must equal T as a multiset.  This is because every column of T′ must lie in the
G-orbit of `eval(repS, E)`, and all `|G| / |stab(repS, E)|` elements of that orbit are
already accounted for in T — there is no room for a spurious column.  The only genuine
non-uniqueness is when two columns of T are identical (non-trivial `stab(repS, E)`), in
which case those columns are indistinguishable and no alignment is needed or possible.

**Algorithm.**

1. **Pre-place constant rows.**  Identify all rows r where the Phase 1 orbit of `fo_r`
   has size 1.  Insert the single value into every partial column.  Remove these rows
   from further consideration.

2. **Seed from the most-informative row.**  Among the remaining rows, choose row r*
   with the most distinct Phase 1 values (equivalently, the smallest maximum repetition
   multiplicity).  Its values seed the partial columns: each distinct value in row r*
   starts one equivalence class of partial columns; the multiplicity of that value
   determines how many columns are in that class.  If all Phase 1 values of row r* are
   distinct, each seeds exactly one column, and the rest of the algorithm proceeds
   without any branching.

3. **Extend one row at a time.**  At each step, choose the next row r to add — prefer
   the row whose pairwise OverlapBlock with the already-placed rows has the most distinct
   values in r, minimising residual ambiguity.  For each partial column (or equivalence
   class of partial columns), look up the known row r* value in `R_{r*, r}` to read off
   the corresponding row r value.

4. **Resolve ambiguities.**  Where a repeated value in an already-placed row creates
   multiple candidate extensions, keep all candidates as branches.  Prune branches using
   additional blocks `R_{i,j}` connecting two already-placed rows i and j: any branch
   that implies a `(row_i value, row_j value)` pair absent from `R_{i,j}` is eliminated.
   Continue until all branches either collapse to a single candidate or are confirmed as
   genuinely identical columns.

5. **Multiplicity bookkeeping.**  The number of table columns represented by each
   unresolved branch is the number of group elements that map to indistinguishable
   columns under all available pairwise constraints — i.e. the size of the residual
   stabiliser.  Record this as the column multiplicity.

**Complexity.**  In the generic (no-collision) case, branching never occurs: all
ambiguities are resolved in step 4 with at most one candidate per partial column.  Total
cost is O(`|repS|` × `|G|`) extensions plus O(`|repS|`²) pairwise lookups — trivially
fast for realistic physics group sizes (`|G|` ~ 10–1000).

When collisions do occur, the branching factor is bounded by the size of the collision
equivalence classes, not by `|G|!`.  The algorithm never enumerates all orderings of the
full column set.

### Numerical robustness: pairings may be wrong even when values are exact

After polishing, **every value in every row is an exact Phase 1 value** — the same
floating-point number will appear in all OverlapBlocks that include that row.  This much
is guaranteed: row-value identity holds exactly.

What polishing does *not* guarantee is that the **pairings within an OverlapBlock are
correct**.  Polishing uses independent Hungarian assignments for the u- and v-rows.  If
two Phase 1 values in a row are very close (near-collision), the Hungarian step may swap
the assignment — mapping the noisy root that should have matched value A onto value B and
vice versa.  The result is a polished OverlapBlock output whose individual values are all
exact Phase 1 values, but whose internal pairings may be locally wrong.

Concretely: suppose the true column alignment has

```
row1row2: { (0, 2.0000001), (5, 1.9999999), (1, 9) }
```

and the Phase 1 values for row 2 are `{2.0, 2.0, 9}` (a genuine Phase 1 near-collision
at 2.0).  Polishing maps both noisy row-2 values to `2.0`, giving
`{ (0, 2.0), (5, 2.0), (1, 9) }` — correct values, but now the two `2.0` entries are
ambiguous.  Meanwhile `row1row3` independently gives `{ (0, 4), (5, -3), (1, 7) }`.  If
the `row2row3` decoder happened to pair row2=2.0 with row3=−3 rather than row3=4, the
join on row2 will produce a conflict with `row1row3`.

**Exact match is a guide, not gospel.**  The alignment algorithm must treat each
OverlapBlock's pairings as strong but fallible evidence, and resolve conflicts by
consulting the full set of letterbox views rather than trusting any single one.

**Conflict resolution principle.**  When a candidate table T produced by the join is
contradicted by one or more R_{ij} views — i.e. a column of T implies a
`(row_i value, row_j value)` pair that does not appear in `R_{ij}` — the algorithm
should not attempt to diagnose *why* the conflict arose.  Instead:

1. Define a **disagreement loss**: for each `(row_i value, row_j value)` pair implied by
   a candidate table T, find the closest pair present in R_{ij} and accumulate the
   squared Euclidean distance between them, summed across all i, j.  A loss of this form
   (analogous to χ²) is recommended as a natural starting point: it weights each
   discrepancy by its magnitude, so many tiny disagreements (distance 0.00001 each) are
   correctly rated as far less serious than a single large one (distance 10000).  Other
   differentiable or rank-based loss functions are possible and may be substituted if
   motivated; the key requirement is that the loss is zero when T is consistent with all
   R_{ij} and increases smoothly with the degree of inconsistency.
2. Enumerate alternative column assignments within the ambiguous equivalence classes
   (those sharing identical row values in at least one row, which bound the search space).
3. Choose the assignment that minimises the disagreement loss.

No causal diagnosis of which OverlapBlock was wrong is needed or attempted.  The
algorithm's job is only to find the alignment that best satisfies all the evidence
simultaneously.

In all cases the correct table is one of the `|G| / |stab(repS, E)|` elements of the
G-orbit of `eval(repS, E)`.  A polishing swap shifts to a neighbouring candidate within
that orbit; the majority-vote correction restores the right one.

### Output type

```
AnnotatedMultisetOfRepSEvalVectors
```

A multiset of tuples, each of length `|repS|`, where position r holds the evaluation of
`repS[r]` (the r-th atom of `repS`) for that group element.  The tuple positions
correspond to known algebraic atoms; the tuple values are floats.  The multiset of
tuples is the G-orbit, and column order within the multiset carries no meaning.

### Correctness check

In round-trip tests, the event E is known.  Ground-truth columns are enumerated directly
as `{ eval(g · repS, E) : g ∈ G }` (as a multiset) and compared against the alignment
decoder's output.  The test passes if and only if the two multisets are equal (using
fuzzy float comparison; exact equality if integer events are used).

### Libraries

For realistic group sizes (`|G|` ~ 10–1000), the join is straightforward to implement in
pure Python using dict-based hash joins — no external library is required.  If group
sizes grow to the point where a more efficient join implementation is warranted, standard
libraries already in the environment (`numpy`, `pandas`) may be used.  Any additional
external library dependency must be low-fragility: a single well-maintained package with
minimal transitive dependencies is acceptable; a dependency tree of dozens of packages is
not.

---

## What the decoded output represents (summary)

After Step A, the caller holds:

1. **The collection of annotated multisets** (one per encoder in the tree), each saying:
   > "These values were produced by evaluating this known set of atoms [or atom-pairs],
   > in some unknown order."

2. **The G-orbit of repS evaluations** (alignment decoder output): the maximal
   G-invariant information that the encoding preserves — the full correlational structure
   of evaluations across all atoms, up to the action of G.

Neither output contains a bijection between atoms and values.  The G-ambiguity (the
unknown permutation) has been deliberately erased by the encoder.  What remains is the
quotient by G — not a specific representative, but the entire orbit.

The fact that values are H-invariant scalars means H has already been quotiented out.
Together: the decoded output represents the original event up to the action of G × H,
which is exactly the ambiguity the encoding was designed to preserve.

Step B (vector-space reconstruction) begins from the alignment decoder's output and
constructs actual vectors, choosing a canonical representative under H.  It is deferred.

---

## Numerical precision and testing

Root-finding is approximate.  Equality checks in decode-validation tests must therefore
use fuzzy comparison (tolerance-based).

For development and debugging, a useful alternative is to populate events with integers
or rationals only.  Evaluable operations on integer inputs produce exact integer or
rational outputs, which survive encoding exactly.  Decoding can then be tested with
exact equality (no fuzzy logic), as a single-switch simplification.  End-to-end
round-trip tests written this way are still fully persuasive as correctness checks.

The one-line switch between exact-integer mode and real-valued mode should be designed
in from the start of the decode test suite.

---

## Deferred items

- **Step B (vector-space reconstruction)**: recovering actual vectors from invariant
  values, up to H.  Specified separately after Step A is complete.

- **SegmentInfo extension**: the precise fields to add to SegmentInfo to carry atom /
  atom-pair information for decoder use.

- **Alignment decoder: join-order heuristic**: the greedy "most distinct values first"
  row-ordering is described in principle; the exact heuristic (tie-breaking, handling of
  multi-group FOs with partial repetition) is to be settled during implementation.

- **Alignment decoder: stabiliser multiplicity precomputation**: the algebraic formula
  for `|stab(u, v)|` as a function of atom-pair flavour is known in principle (it follows
  from the orbit-stabiliser theorem applied within each group factor); the concrete
  implementation is deferred to the alignment decoder coding phase.
