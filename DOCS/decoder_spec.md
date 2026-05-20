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

Within each overlap block, one row-pair was dropped by the complementarity rule because
its contribution is algebraically determined by:
- The Phase 1 decoded multiset for orbit u (or v, depending on block structure), and
- The decoded multisets of all other (ASSOC and NULL_SELF) pairs in the same block.

The overlap block, having decoded all other pairs and holding the Phase 1 decoded output,
reconstructs the NULL_COMP pair's `AnnotatedMultisetOfRealPairs` by reversing the
complementarity relation used during encoding.

The exact complementarity formula is deferred pending confirmation from the encoding
source.

### Interface between Phase 1 and Phase 2

The top-level decoder runs Phase 1 first, producing an `AnnotatedMultisetOfReals` for
each orbit.  These are then passed (dependency injection) to each OverlapBlockDecoder
that needs them.  The OverlapBlockDecoder does not hold a reference to the Phase 1
decoder directly; it receives already-decoded Phase 1 results.

This is the *only* point in the decoder hierarchy where two levels must communicate.

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

## What the decoded output represents (summary)

After Step A, the caller holds a collection of annotated multisets, one per encoder in
the tree.  Each says:

> "These values (stored as a list, treated as a multiset) were produced by evaluating
> this known set of atoms [or atom-pairs], in some unknown order."

The unknown order is the G-ambiguity.  The fact that values are H-invariant scalars
means H has already been quotiented out.  Together: the decoded output represents the
original event up to the action of G × H, which is exactly the ambiguity the encoding
was designed to preserve.

Step B (vector-space reconstruction) begins from these annotated multisets and
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

- **NULL_COMP complementarity formula**: the exact algebraic relation to reverse during
  decoding.  To be added once confirmed from the encoding source code.

- **Step B (vector-space reconstruction)**: recovering actual vectors from invariant
  values, up to H.  Specified separately after Step A is complete.

- **SegmentInfo extension**: the precise fields to add to SegmentInfo to carry atom /
  atom-pair information for decoder use.
