# symatom — Specification

**Status**: draft, pre-implementation  
**Last revised**: May 2026 (rev 6 — added OrbitEnumerator strategy to Section 6; orbit_size/orbit_elements/OrbitEnumerator to Section 8.7)  

---

## 1. Purpose and scope

`symatom` is a symbolic algebra layer for representing and manipulating
*atoms* — multilinear operations applied to labelled vectors drawn from
named types — together with the symmetry machinery needed to canonicalise
them and enumerate their orbits under group actions.

`symatom` is a **building block**. It does not perform any encoding, produce
any numerical output, or know anything about permutation-invariant embeddings.
Those concerns belong to packages that sit above it.

### Explicitly out of scope for this spec

- Numerical evaluation of atoms given concrete vector coordinates (a future
  evaluation hook is anticipated; see Section 9).
- Any specific canonicalisation algorithm (the spec defines the *contract*
  only).
- The polynomial zipping machinery that sits above repS / repL.
- Mixed-symmetry operations (deferred; see Section 4.3).

---

## 2. Concepts and definitions

### 2.1 Vector label

A **vector label** is an opaque, comparable, hashable symbol — in practice a
short string such as `"a"`, `"b"`, `"p"` — that names one vector in an event.
Labels have no intrinsic geometry; they are purely symbolic.

### 2.2 Vector type

A **vector type** is a named, ordered collection of vector labels all of the
same physical type (e.g. all electron momenta, all muon momenta).  

Properties:
- `name: str` — a human-readable identifier, e.g. `"electrons"`.
- `labels: tuple[label, ...]` — the ordered tuple of labels belonging to this
  type, e.g. `("a", "b", "c", "d")`.
- `size: int` — `len(labels)`.

The **symmetric group** $S_n$ (where $n$ = `size`) acts on a vector type by
permuting its labels. This is the only group action a vector type declares.

Multiple vector types may coexist in a calculation. The combined symmetry
group G acting on a collection of vector types is the direct product of their
individual symmetric groups; this object is represented by **TheGroup**
(Section 2.6).

A label belongs to exactly one vector type.

### 2.3 Operation

An **operation** is a named, reusable description of a type of multilinear
contraction. It does **not** know about any specific set of vectors, spatial
dimension, or numerical implementation.

Required properties:
- `name: str` — e.g. `"dot"`, `"eps3"`, `"eps4"`.
- `rank: int` — the number of vector arguments the operation takes (≥ 1).
- `parity: int` — either `+1` or `−1`. Declares how the operation's value
  transforms under spatial parity (inversion of all coordinates). A dot
  product has `parity = +1`; a scalar triple product has `parity = −1`.
- `argument_symmetry: ArgumentSymmetry` — declares how the operation's value
  transforms under permutations of its **own** arguments (see Section 4).

Operations are registered in an **operation registry** that is local to a
`Plan` (Section 6); they are not global singletons.

### 2.4 Atom

An **atom** is an operation applied to an ordered tuple of vector labels,
together with an overall sign.

Properties:
- `operation: Operation`
- `labels: tuple[label, ...]` — length must equal `operation.rank`. Labels may
  come from different vector types. **Note:** the constructor reorders the
  stored labels into a canonical form depending on `argument_symmetry` (see
  Section 3.2); the caller should not assume the stored order matches the
  order passed in.
- `sign: int` — either `+1` or `−1`.

The `sign` field is an intrinsic part of the atom's identity. It is the
**only** place a sign lives; canonicalisation absorbs all sign changes it
introduces directly into the atom's `sign` field and never returns an
additional external sign to the caller. This means two atoms that differ only
in sign — e.g. `(+1, eps, (a,b,p))` and `(−1, eps, (a,b,p))` — are distinct
atoms, not the same atom with an external multiplier.

An atom is **well-formed** if:
1. `len(labels) == operation.rank`,
2. every label belongs to a known vector type in the current context,
3. `sign == +1` whenever `operation.argument_symmetry != ANTISYMMETRIC`, and
4. all labels in the tuple are **distinct**.

Rule 3 is enforced at construction time (raises an error on violation). It
expresses the fact that only antisymmetric operations can ever produce
sign-flipped atoms under argument reordering or orbit action; for symmetric
and unstructured operations `sign = −1` has no legitimate meaning within this
framework and its presence indicates a bug. In particular, atoms representing
dot products must always carry `sign = +1`.

Rule 4 is enforced at construction time (raises an error on violation). Each
operation is understood to take *r* **distinct** vector arguments. A squared
length `a·a` is not a rank-2 dot applied to `(a, a)`; it is a separate rank-1
operation (e.g. `lenSq`) that happens to evaluate the same expression. The
distinct-label rule also ensures that the orbit structure is clean (see
Section 7.2).

An atom is **purely symbolic**. It has no numerical value until an evaluation
hook (Section 8) is invoked.

**Important distinction — internal symmetry vs orbit action:**  
The `argument_symmetry` of an operation describes what happens when the atom's
*own* labels are permuted among themselves — e.g. `dot("a","b") = dot("b","a")`
because dot is symmetric. This is entirely separate from what happens when the
*vector type pool* acts on the atom — e.g. the permutation `a↔c` acting on
`dot("a","b")` replaces `"a"` with `"c"` to give `dot("c","b")`, which is a
*different atom* with a potentially unrelated value. The orbit machinery
(Section 7) tracks the second kind; `argument_symmetry` only informs
canonicalisation within a single atom.

### 2.5 Atom tuple

An **atom tuple** is a finite, ordered sequence of atoms. Atom tuples are the
primary objects that get canonicalised and orbit-computed.

A 1-tuple is just a single atom. The pair-of-rows objects motivating this
package are 2-tuples.

### 2.6 TheGroup

A **TheGroup** represents the full symmetry group G = S_{n₁} × … × S_{nₘ},
the direct product of one symmetric group per vector type.  It is the
group-theory counterpart to the collection of VectorTypes: the VectorTypes
name and enumerate the labels; TheGroup encodes the symmetry that acts on them.
TheGroup acts on both individual atoms and on atom-pairs.

TheGroup stores `types: tuple[VectorType]` and exposes:

*Single-atom methods* (short names — the simpler, more fundamental case):

- `orbit(u) -> list[Atom]` — all atoms in the G-orbit of `u`.
- `orbit_brute(u) -> list[Atom]` — same, brute-force reference implementation.
- `stabiliser_size(u) -> int` — `|Stab_G(u)|`, computed algebraically.
- `orbit_size(u) -> int` — `|G| / |Stab_G(u)|`, computed algebraically.
- `in_orbit(candidate, rep) -> bool` — True iff `candidate ∈ G·rep`.
- `in_orbit_brute(candidate, rep) -> bool` — same, brute-force reference.

*Atom-pair methods* (`_pair` suffix):

- `orbit_pair(u, v) -> list[(Atom, Atom)]` — all atom-pairs in the G-orbit of `(u, v)`.
- `orbit_brute_pair(u, v) -> list[(Atom, Atom)]` — same, brute-force reference.
- `stabiliser_size_pair(u, v) -> int` — `|Stab_G(u, v)|`, computed algebraically.
- `orbit_size_pair(u, v) -> int` — `|G| / |Stab_G(u, v)|`, computed algebraically.
- `in_orbit_pair(candidate, rep) -> bool` — True iff `candidate ∈ G·rep` (both args are `(Atom, Atom)` tuples).
- `in_orbit_brute_pair(candidate, rep) -> bool` — same, brute-force reference.

*Sign correlation (pair only):*

- `sign_correlation_type(u, v) -> SignCorrelationType` — classifies how the G-orbit of `(u, v)` couples the signs of the two atoms.  Algebraic O(n).
- `sign_correlation_type_brute(u, v) -> SignCorrelationType` — same, brute-force reference.

*Group structure:*

- `order() -> int` — `|G| = n₁! × … × nₘ!`.

All `*_brute` methods are permanent O(∏ n_g!) reference implementations kept
for cross-validation; they must not be removed.

TheGroup stores VectorType objects (rather than just sizes) because the
brute-force enumeration methods need actual label sequences to generate
concrete permutations.  The algebraic methods only use per-type label *counts*;
a future step could reduce TheGroup to `TheGroup(sizes: tuple[int])`, but for
now the VectorType tuple is carried directly.

---

## 3. Argument symmetry

### 3.1 Supported symmetry types

The following `ArgumentSymmetry` values are supported in this version:

| Value | Meaning |
|---|---|
| `SYMMETRIC` | The operation's value is unchanged under any permutation of its arguments. E.g. a dot product. |
| `ANTISYMMETRIC` | The operation's value is multiplied by the sign of the permutation under any permutation of its arguments. E.g. an epsilon contraction. |
| `UNSTRUCTURED` | No declared symmetry. The framework makes no assumptions and performs no sign adjustments when permuting arguments. |

### 3.2 Use in canonicalisation

`ArgumentSymmetry` is used in two places: at atom construction time, and
inside canonicalisation.  In both cases the purpose is the same — to put an
atom's own argument list into a canonical order (e.g. alphabetical) and to
absorb any resulting sign change into the atom's `sign` field (e.g. an odd
permutation of an `ANTISYMMETRIC` operation's arguments contributes a factor
of `−1`).  The `Atom` constructor enforces this at creation so that two atoms
that are argument-permutations of each other are immediately equal as Python
objects.  The orbit machinery never uses `ArgumentSymmetry` directly.

### 3.3 Future extension

Mixed symmetry (e.g. symmetric in some argument slots, antisymmetric in
others) is **not** supported in this version but is anticipated as a future
need (e.g. when compound atoms — pairs of atoms treated as first-class atoms —
are introduced). The `ArgumentSymmetry` type is therefore designed as an
extensible type, not a closed enum, so that a `MixedSymmetry` variant can be
added without breaking existing code.

---

## 5. Context

A **context** bundles the label-side and group-theory-side information for a
computation into a single immutable object.  It is a frozen dataclass with two
public attributes, one from each side:

- `types: tuple[VectorType]` (§2.2) — the ordered tuple of vector types in
  scope; the label world.
- `the_group: TheGroup` (§2.6) — the symmetry group G = S_{n₁} × … × S_{nₘ}
  constructed automatically from `types`; the group-theory world.

`the_group` is the authoritative group object for orbit enumeration, stabiliser
computation, and sign-correlation queries.  Callers use `ctx.the_group` directly
rather than constructing TheGroup instances ad hoc.

Context also provides label-lookup utilities:

- `type_of(label) -> VectorType` — returns the type a label belongs to.
- `all_labels` property — the full ordered tuple of labels in scope.

Computations that need a context receive it as an explicit parameter.

---

## 6. Plans

A **plan** bundles together the configuration needed for a particular
computation:

- A **context** (Section 5): the label types in scope and their symmetry group.
- An **operation registry**: the set of named operations available in this
  plan.
- An **orbit enumerator** (Section 8.8): the strategy used to enumerate
  G-orbits of atom-pairs.

A plan is the unit of "local configuration". Two plans may differ in any of
these components. Passing different plans to the same computation function
produces results under different configurations without any global state being
affected.

The orbit enumerator defaults to `DirectOrbitEnumerator` (Section 8.8).
`BruteForceOrbitEnumerator` remains available and is used in tests as the
reference implementation.

---

## 7. Orbits and stabilisers

Given a canonical atom tuple `t` and a plan, the **orbit** of `t` is the set
of all distinct atom tuples reachable by applying elements of the full symmetry
group G = `plan.context.the_group` (see Section 5.1).  Each element of the
orbit is itself in canonical form.

Required operations:

- `orbit(t, plan) -> list[atom_tuple]` — returns the orbit as an ordered list
  of canonical atom tuples (order is deterministic but not otherwise
  specified).
- `stabiliser_size(t, plan) -> int` — returns the size of the stabiliser
  subgroup of `t` (the number of group elements that map `t` to itself). By
  the orbit-stabiliser theorem this equals `|G| / |orbit(t)|`.
- `orbit_and_stabiliser_size(t, plan) -> (list[atom_tuple], int)` — combined
  for efficiency.

**Orbit symmetry vs internal symmetry** (repeated from Section 2.4 for
emphasis): the orbit machinery acts by relabelling vector labels across all
atoms in the tuple simultaneously. It does *not* permute arguments within a
single atom; that is the job of the canonicalisation step.

### 7.1 Counting canonical classes

The distinct-label rule (Section 2.4 Rule 4) makes it straightforward to
count how many distinct canonical forms exist for a single-atom tuple given
an operation of rank *r* and a context with types of sizes *n₁, n₂, …, nₘ*.

Because each group $S_{n_i}$ acts transitively on the $k$-element subsets of
its label set, two single-atom tuples that draw the same number of labels from
each type are always in the same orbit. Two tuples that draw *different*
numbers from some type are always in different orbits (the type index of each
label is fixed by the context). Therefore the canonical classes are indexed
exactly by the **label-composition tuples** $(k_1, k_2, \dots, k_m)$ with:

$$k_i \ge 0, \quad k_i \le n_i, \quad \sum_i k_i = r.$$

The number of canonical classes is the number of such tuples. For example,
with `eps3` (rank 3), electrons (size 4) and muons (size 2) there are exactly
three classes: $(3,0)$, $(2,1)$, $(1,2)$ — corresponding to atoms drawn
entirely from electrons, two electrons and one muon, and one electron and both
muons respectively. The class $(0,3)$ is excluded because there are only two
muon labels.

This count holds for both `SYMMETRIC` and `ANTISYMMETRIC` operations. For
`ANTISYMMETRIC` operations the distinct-label rule is also consistent with
antisymmetry (an antisymmetric form vanishes whenever two arguments are
equal).

### 7.2 Atom utilities

The following utilities on individual atoms are required:

- `are_negatives(a, b) -> bool` — returns `True` if and only if `a` and `b`
  have the same operation and the same labels but opposite signs. Only
  meaningful (and only ever `True`) when `operation.argument_symmetry ==
  ANTISYMMETRIC`; always `False` for symmetric or unstructured operations,
  consistent with the well-formedness rule that such atoms always carry
  `sign = +1`. Callers can use this to ask "would a parity flip (or an odd
  argument permutation) map atom `a` to atom `b`?" without needing to perform
  the transformation themselves.

---

## 8. Flavour, Ingredients, and the repL / repS representations

### 8.1 Motivation

Given a context (e.g. Electrons `{a,b,c,d}` + Muons `{p,q}`) and a set of
operations (e.g. `mass`, `dot`, `eps3`), we want to enumerate all atoms that
can be formed by applying each operation to every valid combination of labels.
These atoms are the basis elements from which representations of an event are
built.

The enumeration naturally falls into clusters determined solely by the operation
and by how many labels come from each species. These clusters — not the
individual atoms — are the mathematically meaningful objects. Knowing the
cluster tells you immediately how many atoms it contains and how to generate
them; the individual atoms within the cluster are almost secondary.

### 8.2 Flavour

A **Flavour** is an ordered tuple of non-negative integers, one per vector
type in the context, encoding how many arguments an operation draws from each
type. For a context with types `G_1, …, G_m` and an operation of rank *r*:

- `Flavour((k_1, …, k_m))` means *k_i* arguments come from type *G_i*.
- Validity requires `k_i ≥ 0`, `k_i ≤ |G_i|` (type size), and
  `Σ k_i = r`.

Example: with types `[Electrons, Muons]` and `eps3` (rank 3), the valid
Flavours are `(3,0)`, `(2,1)`, and `(1,2)`. The Flavour `(0,3)` is excluded
because there are only two muon labels.

A Flavour is **context-agnostic**: it is just a tuple of counts. The context
is carried by the `FlavouredOperator` (Section 8.4).

*Why "Flavour"?* The metaphor is culinary: the species composition of a cluster
of atoms is its "taste" — all electrons and no muons is a very different
flavour from a mix of electrons and muons, just as different particle types
play different physical roles. Within one Flavour, all atoms taste the same;
individual atoms differ only in which specific labels were chosen (the
Ingredients, Section 8.3).

### 8.3 Ingredients

An **Ingredients** is an ordered list of labels consistent with a given
Flavour. For example, `[a, b, p]` is one Ingredients of Flavour `(2,1)`;
`[a, c, p]` is another.

Key properties:

- **Ordering matters.** Ingredients are passed directly to the operation as
  positional arguments. For an `ANTISYMMETRIC` operation, the permutation
  that sorts the ingredients into canonical order contributes a sign, which
  is absorbed into the resulting atom's `sign` field. This property becomes
  important when mixed-symmetry operations are introduced.
- **Transient.** Ingredients are used during atom construction and are not
  persisted. They are not a named type in the public API; they appear
  implicitly in the generation algorithm.

### 8.4 FlavouredOperator

A **FlavouredOperator** bundles one operation with one Flavour in a given
context. It is a lazy recipe for a cluster of atoms: it knows how many atoms
it contains and can generate them on demand, but it does not materialise them
unless asked.

Properties:

- `operation` — the `Operation` instance.
- `flavour` — the `Flavour` instance.
- `context` — the `Context` (so it knows the label pools).
- `signed: bool` — `True` for repL semantics, `False` for repS (Section 8.5).

Required methods:

- `.count_of_atoms_one_per_sign() -> int` — returns the number of atoms this cluster contains, by
  pure combinatorics. No atoms are generated.
- `.atoms_one_per_sign() -> iterator` — lazily generates all atoms for this cluster. Labels
  are generated as **combinations** (one per choice of labels from each
  type's pool), not as permutations-then-deduplicate. Permutation-based
  generation is acceptable in unit tests to cross-check, but the direct
  combination strategy is preferred in production because species pools may
  be large.
- `.matches_ignoring_sign(atom) -> bool` — returns `True` if the atom belongs to the
  vocabulary of this `FlavouredOperator` (ignoring sign). For `SYMMETRIC` and `ANTISYMMETRIC`
  operations this can usually be determined without full enumeration (see
  below).
- `.canonical_representative() -> Atom` — returns the first atom from
  `.atoms_one_per_sign()` as the orbit representative for this `FlavouredOperator`.  All
  atoms in a `FlavouredOperator` share the same operation and Flavour and
  therefore lie in the same G-orbit, so any choice of representative is equally
  valid.

Mixed-symmetry operations are **not supported** in this version and are not
yet handled by `FlavouredOperator`. Support is deferred to a future version.

**Membership without full enumeration.** For `SYMMETRIC` operations,
membership reduces to: is the operation correct, do the labels match the
flavour, and is the sign `+1`? For `ANTISYMMETRIC` operations, membership
also allows `sign = −1` in repL mode. Full enumeration is not required for
either case.

### 8.5 Sign semantics: repS vs repL

Within a `FlavouredOperator`:

- **repS** (`signed=False`): one atom per label combination (one per
  Ingredients orbit under the operation's own argument symmetry). Exactly
  one sign is included per combination — always `sign = +1` on entry, but
  the `Atom` constructor immediately applies internal argument
  canonicalization (Section 2.4): for `SYMMETRIC` and `ANTISYMMETRIC`
  operations the stored labels are globally sorted (with any sign flip
  absorbed into the sign field); for `UNSTRUCTURED` operations the labels
  are kept in group-concatenation order, sorted within each group but not
  across groups.
- **repL** (`signed=True`): same as repS, plus the additive inverse for each
  `ANTISYMMETRIC` atom. So both `sign = +1` and `sign = −1` appear for each
  combination. For `SYMMETRIC` and `UNSTRUCTURED` operations, repL and repS
  are identical.

A convenient mnemonic: in repL, every antisymmetric atom travels with its
negative twin; in repS, the twin is suppressed.

The `count()` values follow directly:

| operation symmetry | repS count | repL count |
|---|---|---|
| SYMMETRIC or UNSTRUCTURED | `∏ C(nᵢ, kᵢ)` | same |
| ANTISYMMETRIC | `∏ C(nᵢ, kᵢ)` | `2 × ∏ C(nᵢ, kᵢ)` |

### 8.6 repL and repS functions

`repL(context, operations)` and `repS(context, operations)` are functions
(not types) that enumerate all valid `(operation, Flavour)` combinations for
the given context and return a list of `FlavouredOperator` instances.

"Valid" means: each species contributes no more labels than its pool size
allows, i.e. `k_i ≤ |G_i|` for all *i*, and `Σ k_i = rank(operation)`.

The returned list covers every combination of every operation with every
valid Flavour. The individual atoms in any cluster are available lazily via
`.atoms_one_per_sign()` and are never materialised unless explicitly requested.

### 8.7 PairFlavour

A **PairFlavour** is to a pair of atoms what a `Flavour` is to a single atom:
it is the orbit-type of the pair under the full symmetry group G.

**Theorem**: two atom-pairs `(u, v)` and `(u', v')` (with all labels distinct,
as guaranteed by the well-formedness rules) lie in the same G-orbit if and
only if they have the same PairFlavour.  This makes the set of distinct
PairFlavours a direct index of the distinct G-orbits of atom-pairs — which is
exactly the set of non-redundant row-pairs needed for the encoding layer above.

A PairFlavour bundles:

- `op_u: Operation`, `flavour_u: Flavour` — operation and flavour of atom *u*
- `op_v: Operation`, `flavour_v: Flavour` — operation and flavour of atom *v*
- `overlap: tuple[int, ...]` — one non-negative integer per type: the number
  of labels shared between *u* and *v* from that type.

**Validity** (enforced at construction, context-free):

- `len(flavour_u.counts) == len(flavour_v.counts) == len(overlap)`
- `0 ≤ overlap[i] ≤ min(flavour_u.counts[i], flavour_v.counts[i])` for each *i*

The type-size upper bound `flavour_u.counts[i] + flavour_v.counts[i] −
overlap[i] ≤ |G_i|` requires knowing the type sizes and is not validated at
construction; `canonical_pair_flavours` enforces it implicitly by only
generating overlaps that satisfy it.

**Canonical ordering**: PairFlavour imposes a canonical ordering on the
`(op_u, flavour_u)` vs `(op_v, flavour_v)` sides at construction time, so
that swapping the two sides always produces the same object.  This mirrors
the fact that the `SimpleCanonicaliser` sorts atoms within a tuple, meaning
the G-orbit of `(u, v)` and the G-orbit of `(v, u)` are identified by the
canonicaliser — and therefore need not be encoded separately.

Note: sign information (`sign = +1` vs `sign = −1` for ANTISYMMETRIC atoms)
is deliberately absent from PairFlavour.  Sign-related encoding redundancies
— such as the fact that zipping `(rowA, rowB)` and zipping `(rowA, −rowB)`
carry the same information — are encoding-layer concerns, not group-theory
concerns, and are handled separately above `symatom`.

**`count(type_sizes: tuple[int, ...]) -> int`**: returns the number of
ordered atom-pairs `(u, v)` with this PairFlavour, given the sizes of the
vector types.  Pure combinatorics; no atoms are materialised.  For type *i*
with size *nᵢ*, flavour counts *kᵤᵢ* and *kᵥᵢ*, and overlap *sᵢ*:

    ways_i = C(nᵢ, sᵢ) × C(nᵢ − sᵢ, kᵤᵢ − sᵢ) × C(nᵢ − kᵤᵢ, kᵥᵢ − sᵢ)

The total count is the product over all types.  This counts *ordered* pairs;
when `op_u == op_v` and `flavour_u == flavour_v`, each unordered pair
`{u, v}` with `u ≠ v` appears twice.

**`pair_flavour_of(atom_u, atom_v, context) -> PairFlavour`**: computes the
PairFlavour of a concrete atom-pair.

**`canonical_pair_flavours(fo_list, context) -> list[PairFlavour]`**: returns
all distinct PairFlavours for the given list of `FlavouredOperator`s and
context, generated directly from `(FlavouredOperator, FlavouredOperator,
overlap)` triples without materialising atom-pairs.  For each ordered pair of
FlavouredOperators `(fo_u, fo_v)`, the valid overlap range for type *i* is:

    s_lo = max(0, k_u_i + k_v_i − n_i),   s_hi = min(k_u_i, k_v_i)

PairFlavour's own canonical ordering ensures that `(fo_u, fo_v)` and
`(fo_v, fo_u)` contribute identical PairFlavours; a `seen` set handles
deduplication.  Returns a deterministically sorted list.

**`brute_force_canonical_pair_flavours(fo_list, context) -> set[PairFlavour]`**:
reference / cross-check implementation.  Materialises all atom-pairs and
extracts their PairFlavours.  *O*(*N*²) where *N* = total atoms.  Not suitable
for large contexts; intended only for unit-test cross-validation against
`canonical_pair_flavours`.

**`orbit_size(type_sizes: tuple[int, ...]) -> int`**: returns the number of
distinct atom-pairs in the G-orbit of the canonical representative of this
PairFlavour.  This equals `count(type_sizes)` multiplied by 2 for each
ANTISYMMETRIC operation in the pair.  The extra factor arises because a
permutation that reorders the arguments of an antisymmetric atom flips its
sign — so both the `sign=+1` and `sign=−1` variants of every label combination
appear as distinct elements of the true G-orbit.

**`orbit_elements(context) -> list[tuple[Atom, Atom]]`**: materialises and
returns all atom-pairs in the G-orbit of the canonical representative.  Uses
`FlavouredOperator.atoms_one_per_sign()` (with `signed=True`) to enumerate all concrete
atoms for each side, then filters by `pair_flavour_of(u, v, context) == self`.
The length of the returned list always equals `orbit_size(type_sizes)`.  This
is the primary entry point for the encoding layer: iterate over the returned
pairs, evaluate both atoms numerically, and form the complex numbers
`z_k = eval(u') + i·eval(v')` to pass to the polynomial embedder.

### 8.8 OrbitEnumerator

Enumerating the G-orbit of an atom-pair is a performance-critical operation:
for large particle multiplicities (*n* ~ 30 per species) the cost per orbit
dominates the encoding.  `symatom` therefore exposes orbit enumeration as a
pluggable **strategy**, following the same pattern as `Canonicaliser`.

**`OrbitEnumerator`** (abstract base class, `symatom/orbit_enum.py`): defines
the interface:

```python
def orbit_elements(self, pf: PairFlavour, context: Context) -> list[tuple[Atom, Atom]]
```

**`BruteForceOrbitEnumerator`**: delegates to `pf.orbit_elements(context)` —
the cross-product + `pair_flavour_of` filter described above.  Correct by
inspection.  O(`atoms_u` × `atoms_v`).  **Kept permanently** as the reference
implementation for testing.

**`DirectOrbitEnumerator`**: O(orbit\_size) direct combinatorial construction
— no filtering step anywhere.  For each type *i*, enumerates C(*n*ᵢ, *s*ᵢ)
choices of shared labels, then C(*n*ᵢ − *s*ᵢ, *k*ᵤᵢ − *s*ᵢ) u-only labels,
then C(*n*ᵢ − *k*ᵤᵢ, *k*ᵥᵢ − *s*ᵢ) v-only labels; takes the Cartesian
product across types; assembles label tuples and constructs atoms directly.
The `Atom` constructor sorts labels and adjusts sign via `perm_sign`, so
unsorted concatenation order is handled correctly.  Both `sign=+1` and
`sign=-1` variants are emitted for ANTISYMMETRIC operations, giving all four
sign combinations for ANTISYMMETRIC×ANTISYMMETRIC pairs.  This is the default
`orbit_enumerator` in `Plan` and the production implementation for large *n*.

The parametrised tests in `symatom/tests/test_orbit_enum.py` run both
enumerators side-by-side; additional cross-comparison tests assert that
`DirectOrbitEnumerator` produces the same element set as `BruteForceOrbitEnumerator`
for every `PairFlavour`.

---

## 9. Evaluation hook (anticipated, not specified here)  <!-- was §8 -->

It is anticipated that each `Operation` instance will carry an optional
**evaluation method** of the form

```
evaluate(labels: tuple[label, ...], vectors: dict[label, array]) -> float
```

that, given a mapping from labels to concrete numerical vectors, returns the
numerical value of the atom. The sign of the atom is applied by the caller.

The interface of this method is **not** specified in this version of the spec.
It is noted here so that the data structures chosen for `Operation` leave room
for it without requiring a refactor.

---

## 11. What this spec does not decide

The following are left to the implementation:

- The concrete data structures for labels, atoms, and atom tuples (named
  tuples, dataclasses, frozen dataclasses, etc.).
- Whether orbits are computed eagerly or lazily.
- The file and module structure within the `symatom` package.
- How operations are registered (decorator, explicit registry call, etc.).
- The concrete form of the evaluation hook.

---

## 12. Relation to the broader project

`symatom` sits below a not-yet-named encoding package that will use it to:

1. Build repS and repL representations of events.
2. Identify non-redundant canonical atom-tuple pairs.
3. Compute the orbit of each such pair.
4. Feed those orbits into a polynomial zipping step based on the Echinocoder
   polynomial encoder.

`symatom` has no knowledge of any of this. It provides atoms, groups,
canonicalisation, and orbits. The encoding logic is entirely the responsibility
of the layer above.
