# symatom — Specification

**Status**: draft, pre-implementation  
**Last revised**: May 2026 (rev 2 — sign lives in atom, not external to canon)  

---

## 1. Purpose and scope

`symatom` is a symbolic algebra layer for representing and manipulating
*atoms* — multilinear operations applied to labelled vectors drawn from
named groups — together with the symmetry machinery needed to canonicalise
them and enumerate their orbits under group actions.

`symatom` is a **building block**. It does not perform any encoding, produce
any numerical output, or know anything about permutation-invariant embeddings.
Those concerns belong to packages that sit above it.

### Explicitly out of scope for this spec

- Numerical evaluation of atoms given concrete vector coordinates (a future
  evaluation hook is anticipated; see Section 8).
- Any specific canonicalisation algorithm (the spec defines the *contract*
  only).
- The repS / repL representations and the polynomial zipping machinery that
  uses them.
- Mixed-symmetry operations (deferred; see Section 4.3).

---

## 2. Concepts and definitions

### 2.1 Vector label

A **vector label** is an opaque, comparable, hashable symbol — in practice a
short string such as `"a"`, `"b"`, `"p"` — that names one vector in an event.
Labels have no intrinsic geometry; they are purely symbolic.

### 2.2 Vector group

A **vector group** is a named, ordered collection of vector labels all of the
same physical type (e.g. all electron momenta, all muon momenta).  

Properties:
- `name: str` — a human-readable identifier, e.g. `"electrons"`.
- `labels: tuple[label, ...]` — the ordered tuple of labels belonging to this
  group, e.g. `("a", "b", "c", "d")`.
- `size: int` — `len(labels)`.

The **symmetric group** $S_n$ (where $n$ = `size`) acts on a vector group by
permuting its labels. This is the only group action a vector group declares.

Multiple vector groups may coexist in a calculation. The full symmetry group
acting on a collection of vector groups is the direct product of their
individual symmetric groups.

A label belongs to exactly one vector group.

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
  come from different vector groups.
- `sign: int` — either `+1` or `−1`.

The `sign` field is an intrinsic part of the atom's identity. It is the
**only** place a sign lives; canonicalisation absorbs all sign changes it
introduces directly into the atom's `sign` field and never returns an
additional external sign to the caller. This means two atoms that differ only
in sign — e.g. `(+1, eps, (a,b,p))` and `(−1, eps, (a,b,p))` — are distinct
atoms, not the same atom with an external multiplier.

An atom is **well-formed** if:
1. `len(labels) == operation.rank`, and
2. every label belongs to a known vector group in the current context.

An atom is **purely symbolic**. It has no numerical value until an evaluation
hook (Section 8) is invoked.

**Important distinction — internal symmetry vs orbit action:**  
The `argument_symmetry` of an operation describes what happens when the atom's
*own* labels are permuted among themselves — e.g. `dot("a","b") = dot("b","a")`
because dot is symmetric. This is entirely separate from what happens when the
*vector group pool* acts on the atom — e.g. the permutation `a↔c` acting on
`dot("a","b")` replaces `"a"` with `"c"` to give `dot("c","b")`, which is a
*different atom* with a potentially unrelated value. The orbit machinery
(Section 7) tracks the second kind; `argument_symmetry` only informs
canonicalisation within a single atom.

### 2.5 Atom tuple

An **atom tuple** is a finite, ordered sequence of atoms. Atom tuples are the
primary objects that get canonicalised and orbit-computed.

A 1-tuple is just a single atom. The pair-of-rows objects motivating this
package are 2-tuples.

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

`ArgumentSymmetry` is used *only* inside canonicalisation, to put an atom's
own argument list into a canonical order (e.g. alphabetical). Any sign change
that results from reordering the arguments (e.g. an odd permutation of an
`ANTISYMMETRIC` operation's arguments contributes a factor of `−1`) is
multiplied into the atom's `sign` field. The orbit machinery never uses
`ArgumentSymmetry` directly.

### 3.3 Future extension

Mixed symmetry (e.g. symmetric in some argument slots, antisymmetric in
others) is **not** supported in this version but is anticipated as a future
need (e.g. when compound atoms — pairs of atoms treated as first-class atoms —
are introduced). The `ArgumentSymmetry` type is therefore designed as an
extensible type, not a closed enum, so that a `MixedSymmetry` variant can be
added without breaking existing code.

---

## 4. Canonicalisation

### 4.1 What canonicalisation is

Canonicalisation maps an atom tuple to a **canonical form**: a chosen
representative of its equivalence class under the action of the full symmetry
group. Two atom tuples are considered equivalent if one can be obtained from
the other by:

1. Applying a permutation drawn from the product of the vector groups' symmetric
   groups (relabelling vectors), **and**
2. Adjusting signs consistently with the `argument_symmetry` and `parity`
   declarations of each operation involved.

The canonical form is itself a valid atom tuple. All sign changes introduced
during canonicalisation are absorbed into the `sign` fields of the atoms in
the canonical form; no sign is returned separately to the caller.

### 4.2 The canonicalisation contract

Any conforming canonicalisation implementation must satisfy all of the
following:

**C1 — Idempotent**: `canon(canon(x)) == canon(x)` for all atom tuples `x`.

**C2 — Representative**: `x` and `canon(x)` lie in the same orbit under the
full symmetry group.

**C3 — Consistent**: `canon(x) == canon(y)` if and only if `x` and `y` lie in
the same orbit.

**C4 — Sign consistency**: `canon(x)` returns a canonical atom tuple (not a
`(tuple, sign)` pair). All sign information is absorbed into the atoms'
`sign` fields. The consistency requirement is: if atom tuple `y` was obtained
from atom tuple `x` by applying a group element `g` (relabelling vector
labels), then the net sign of `canon(y)` relative to `canon(x)` equals the
sign by which `g` transforms `x` into `y`. Concretely, if evaluating `x` and
`y` on concrete vectors would give values related by a factor `s(g,x) ∈
{+1,−1}`, then the product of the `sign` fields of the atoms in `canon(y)`
equals `s(g,x)` times the product of the `sign` fields of the atoms in
`canon(x)`.

**C5 — Deterministic**: `canon(x)` returns the same result on every call with
the same input.

### 4.3 What the spec does NOT mandate

- Any particular ordering of labels (alphabetical, appearance-order, etc.).
- Any particular ordering of atoms within a tuple.
- Any particular strategy for breaking ties.

These are implementation choices. Different implementations satisfying C1–C5
are all conforming, and the test suite (Section 9) tests the contract, not the
specific canonical form produced.

### 4.4 Canonicalisation plans

A **canonicalisation plan** is an explicit object encapsulating a chosen
canonicalisation implementation together with any configuration it requires.
Plans are **local**: they are passed explicitly to every function that needs
to canonicalise. There is no global or default plan.

This allows multiple plans to be run side-by-side in the same process — for
example, running a repS-based plan and a repL-based plan simultaneously to
cross-validate results.

**Vocabulary.** A plan implicitly declares a **vocabulary**: the set of signed
atoms that are considered named, first-class members of the representation.
This is what distinguishes a repL-style plan from a repS-style plan:

- A *repL-style* vocabulary includes both `(+1, op, labels)` and
  `(−1, op, labels)` for antisymmetric operations. A vector-group permutation
  that would flip the sign of such an atom merely maps it to another named
  vocabulary member — no sign "escapes" the vocabulary.
- A *repS-style* vocabulary includes only `(+1, op, labels)` for antisymmetric
  operations. A sign-flipping permutation produces `(−1, op, labels)`, which
  is a valid signed atom but not a positively-signed vocabulary member.

Both styles use the same signed-atom representation and the same
canonicalisation and orbit machinery. The vocabulary is a property of the
layer above `symatom`; it is recorded here only to explain why two plans that
differ only in vocabulary are expected to produce encodings of the same final
size once all canonicalisations and redundancy-eliminations are applied.

---

## 5. Vector group context

A **context** is a collection of vector groups that are in scope for a given
computation. It provides:

- Lookup: given a label, which group does it belong to?
- Membership: the full set of labels in scope.
- Group product: the full symmetry group is the direct product of the
  symmetric groups of all groups in the context.

A context is immutable once created. Computations that need a context receive
it as an explicit parameter.

---

## 6. Plans

A **plan** bundles together the configuration needed for a particular
computation:

- A **context** (Section 5): which vector groups are in scope.
- A **canonicalisation implementation**: any conforming implementation of the
  canonicalisation contract (Section 4.2).
- An **operation registry**: the set of named operations available in this
  plan.

A plan is the unit of "local configuration". Two plans may differ in any of
these components. Passing different plans to the same computation function
produces results under different configurations without any global state being
affected.

---

## 7. Orbits and stabilisers

Given a canonical atom tuple `t` and a plan, the **orbit** of `t` is the set
of all distinct atom tuples reachable by applying elements of the full symmetry
group (the product of the symmetric groups of all vector groups in the
context). Each element of the orbit is itself in canonical form.

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

---

## 8. Evaluation hook (anticipated, not specified here)

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

## 9. Test contracts

The test suite must include tests for the canonicalisation contract properties
C1–C5 (Section 4.2) against any conforming implementation. These tests should:

- **C1**: Apply `canon` twice and assert the result equals applying it once.
- **C2**: Assert that `canon(x)` and `x` produce the same orbit (i.e. `x`
  appears in `orbit(canon(x), plan)`).
- **C3**: For pairs `(x, y)` known to be in the same orbit, assert
  `canon(x) == canon(y)`. For pairs known to be in different orbits, assert
  `canon(x) != canon(y)`.
- **C4**: For `y` obtained from `x` by a known group element `g` with known
  sign `s(g,x)`, assert that the product of the `sign` fields of atoms in
  `canon(y)` equals `s(g,x)` times the product of the `sign` fields of atoms
  in `canon(x)`.
- **C5**: Call `canon(x)` twice in the same session and assert identical
  results.

Tests should be written against the **contract**, using the plan mechanism to
inject whichever canonicalisation implementation is under test. This ensures
that swapping in a new implementation is immediately testable.

---

## 10. What this spec does not decide

The following are left to the implementation:

- The concrete data structures for labels, atoms, and atom tuples (named
  tuples, dataclasses, frozen dataclasses, etc.).
- The specific canonicalisation algorithm.
- Whether orbits are computed eagerly or lazily.
- The file and module structure within the `symatom` package.
- How operations are registered (decorator, explicit registry call, etc.).
- The concrete form of the evaluation hook.

---

## 11. Relation to the broader project

`symatom` sits below a not-yet-named encoding package that will use it to:

1. Build repS and repL representations of events.
2. Identify non-redundant canonical atom-tuple pairs.
3. Compute the orbit of each such pair.
4. Feed those orbits into a polynomial zipping step based on the Echinocoder
   polynomial encoder.

`symatom` has no knowledge of any of this. It provides atoms, groups,
canonicalisation, and orbits. The encoding logic is entirely the responsibility
of the layer above.
