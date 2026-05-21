# The σ-Symmetry of Self-Pairing Associations

## Summary

Every surviving (non-null) ASSOC segment in a **self-pairing block** (same
operation, same flavour: `op_u == op_v`, `F_u == F_v`) carries a hidden
factor-of-2 redundancy.  The redundancy arises from a swap symmetry on
atom-pairs that forces the z-value multiset to be closed under
`z → i·conj(z)` (denoted σ).  Exploiting this halves the storage of every
such ASSOC.

Note on scope: NULL_SELF and NULL_COMP segments are already zero-length and
are unaffected (halving zero gives zero, which is already optimal).  Half-sort
compression acts on Phase 1 ORBIT segments and is entirely independent —
σ-compression and half-sort operate on disjoint parts of the encoding.

Self-pairing blocks cannot contain TYPE_12 or TYPE_21 orbits; only TYPE_11,
TYPE_NEG, or TYPE_22 can arise (see Section 2a).  This follows from the
swap-closure of the achievable sign set and is independent of argument
symmetry.

---

## 1. The Swap Symmetry

Fix a self-pairing block: `op_u = op_v = op`, `F_u = F_v = F`.
The u-orbit and v-orbit are the same orbit of atoms.

**Lemma.** If `(u_k, v_j)` is an atom-pair with both atoms having `sign = +1`
that appears in any association of this block, then `(v_j, u_k)` is also an
atom-pair with both atoms having `sign = +1` in the **same** association.

*Proof.*
- Overlap count: the number of shared labels between `u_k` and `v_j` equals
  the number shared between `v_j` and `u_k` (it is a symmetric count), so the
  swapped pair falls in the same association.
- Signs: `u_k` and `v_j` each have `sign = +1` by assumption; in the swapped
  pair the roles are exchanged but the signs are unchanged. ∎

**Corollary (σ-closure).** For any atom-pair representative `(u_k, v_j)` with
both signs `+1`:

```
z_{kj} = eval(u_k) + i·eval(v_j)
z_{jk} = eval(v_j) + i·eval(u_k) = i · conj(z_{kj})
```

Both `z_{kj}` and `i·conj(z_{kj})` appear in the set of z-values for
positive-sign representatives.  Therefore the multiset `{z_k}` is closed under
`σ: z ↦ i·conj(z)`.

Note: σ² = identity (apply twice: `i·conj(i·conj(z)) = i·(−i)·z = z`), so σ
is an involution.  Generic elements have σ-orbit of size 2; the fixed point
condition is `z = i·conj(z)` i.e. `Re(z) = Im(z)`, i.e. `z ∝ (1+i)`, which
requires `eval(u_k) = eval(v_j)` — a measure-zero condition for generic events.

> **Note on TYPE_11:** The lemma assumes both atoms have `sign = +1`.
> For TYPE_11 self-pairing this holds by definition: the achievable sign set
> is `{(+,+)}`, so every orbit element has `sign = +1`.  For all currently
> implemented SYMMETRIC operations this is also the empirical case, but
> TYPE_11 is defined by the orbit's achievable sign set, not by argument
> symmetry; future operation types may produce TYPE_11 orbits with other
> argument symmetries.
>
> **Note on TYPE_NEG (form-1 and form-2 both arise):** `_try_neg_partition`
> selects reps with `u.sign = +1`; `v.sign` may be `+1` (form-1) or `−1`
> (form-2).  Both forms arise empirically in self-pairing blocks.
>
> For form-2 rep `(+u, −v)`: `z = eval(u) − i·eval(v) = a − ib`.
> The swap partner is `(+v, −u)`: `z' = eval(v) − i·eval(u) = b − ia`.
> Then `z'^2 = b^2 − a^2 − 2abi = −conj(a^2 − b^2 − 2abi) = −conj(z^2) = τ(z^2)` ✓.
>
> So for form-2, the **swap** partner contributes `τ(z^2)`, not `i·conj(z)`.
> The σ-closure on `{z_k}` does not hold for form-2, but `{z_k^2}` is still
> τ-closed (anti-conj-closed) because of the swap symmetry.  The TYPE_NEG
> σ-encoder works correctly by embedding `{z_k^2}` directly without needing
> σ-closure on `{z_k}`.

---

## 2. Which Associations Does This Affect?

### 2a. Which orbit types can appear in self-pairing blocks?

**TYPE_12 and TYPE_21 are rigorously excluded.** The argument is a swap-closure
proof on the achievable sign set.

In a self-pairing block the roles of u and v are interchangeable: any group
permutation σ that maps `(u_k, v_j)` to another orbit element simultaneously
maps `(v_j, u_k)` to another orbit element in the same association (by the
symmetric overlap-count argument of Section 1).  Therefore the achievable sign
set must be **closed under the swap `(s_u, s_v) → (s_v, s_u)`**.

Checking each type:

| Type | Achievable set | Closed under swap? |
|------|----------------|--------------------|
| TYPE_11 | `{(+,+)}` | ✓ |
| TYPE_NEG | `{(+,+),(−,−)}` | ✓ |
| TYPE_22 | `{(+,+),(+,−),(−,+),(−,−)}` | ✓ |
| TYPE_12 | `{(+,+),(+,−)}` | ✗  (swap gives TYPE_21) |
| TYPE_21 | `{(+,+),(−,+)}` | ✗  (swap gives TYPE_12) |

**Conclusion: self-pairing blocks can only be TYPE_11, TYPE_NEG, or TYPE_22.**
TYPE_12 and TYPE_21 are impossible regardless of operation, group size, or
overlap.

### 2b. Which types actually arise?

- **TYPE_11**: achievable sign set `{(+,+)}` — neither sign can ever flip.
  For all currently implemented operations this arises when both operations
  are SYMMETRIC (e.g. dot × dot), but TYPE_11 is defined by the orbit
  structure, not by argument symmetry.
- **TYPE_NEG**: ANTISYMMETRIC operation, confirmed in tested contexts (single
  group, overlap below maximum). Physically: the shared labels force correlated
  sign flips.
- **TYPE_22**: ANTISYMMETRIC operation, expected for large enough groups where
  u and v have sufficiently non-overlapping labels that their sign flips can be
  driven independently. **Not yet confirmed empirically** — see Open Question 2.

> **⚠ SUSPECT:** The claim "ANTISYMMETRIC self-pairing always gives TYPE_NEG"
> is only confirmed for small groups (4 electrons, overlap=(2,)).  TYPE_22 is
> expected to appear for, e.g., eps3((3,)) × eps3((3,)) with overlap=(1,) in a
> group of ≥5 labels, where the exclusive labels of u and v are disjoint and
> can be permuted independently.

### 2c. Scope of σ-compression

σ-compression has been argued for TYPE_11 and TYPE_NEG self-pairing (Sections
4 and 5).  **TYPE_22 self-pairing is not yet covered** — its achievable set is
symmetric under swap but the σ-compression argument has not been developed for
it.  See Open Question 2.

- **Non-self-pairing blocks**: no swap symmetry. No σ-closure.

---

## 3. Algebraic Structure

Let `P(x) = ∏(x + z_k)` be the characteristic polynomial of the orbit
representatives (n elements, roots at `{-z_k}`).  Since `{z_k}` is σ-closed,
the roots come in σ-pairs `(-z, -iz*)`.

**Claim:** `conj(P(i·conj(x))) = (-1)^n · P(x)`.

> **⚠ SUSPECT — proof not complete in this document.**  A pair-by-pair
> calculation for the SS/TYPE_11 rotation approach below gives a clean
> numerical verification (Section 4), which is taken as sufficient evidence.
> A fully rigorous proof of this functional equation has not been written out.

---

## 4. The Rotation Trick (TYPE_11 Self-Pairing Case)

**Goal:** compress n complex z-values (σ-closed) to n reals.

**Step 1.** Rotate: `w_k = z_k · e^{−iπ/4} = z_k · (1−i)/√2`.

**Step 2.** The set `{w_k}` is conjugate-closed.

*Proof (numerical check).* Take z = 2+i, so z* = 2−i, σ(z) = i(2−i) = 1+2i.
Rotate: w = (2+i)·(1−i)/√2 = (3−i)/√2.  conj(w) = (3+i)/√2.
σ(z)·e^{−iπ/4} = (1+2i)·(1−i)/√2 = (3+i)/√2 = conj(w). ✓

So indeed `σ(z)·e^{−iπ/4} = conj(w)`, confirming `{w_k}` is conjugate-closed.

**Step 3.** The polynomial `Q(t) = P(t·e^{iπ/4})` has real coefficients (since
its roots `{-w_k}` are conjugate-closed).  A real monic polynomial of degree n
has n independent real coefficients.

**Step 4.** Embedding algorithm:
```python
rotation = np.exp(-1j * np.pi / 4)          # (1-i)/√2
w = z_reps * rotation                        # n conjugate-closed values
coeffs, _, _ = zip_embed(w)                  # degree-n poly; roots are conj-paired
r = coeffs.real                              # n reals (imag parts ≈ 0)
return r[0::2] + 1j * r[1::2]               # pack as n/2 complex → n reals total
```

Output: **n reals** (vs. current 2n from Type11PairEncoder).  Factor-2
compression applying to TYPE_11 self-pairing pairs that currently have no
compression at all.

---

## 5. The TYPE_NEG Self-Pairing Case

The NegPairEncoder (for both self-pairing and non-self-pairing TYPE_NEG orbits)
already forms `w_k = z_k²` and embeds the n values directly as n complex
coefficients → 2n reals.

For a self-pairing TYPE_NEG block, both form-1 and form-2 reps arise (see
Section 1 note).  In both cases the swap symmetry guarantees `{z_k^2}` is
anti-conj-closed (τ-closed):

- Form-1 rep `(+u,+v)`: `z = a+ib`.  Swap gives `z' = b+ia`.
  `z'^2 = −conj(z^2) = τ(z^2)` ✓
- Form-2 rep `(+u,−v)`: `z = a−ib`.  Swap gives `z' = b−ia`.
  `z'^2 = −conj(z^2) = τ(z^2)` ✓

Since `{z_k^2}` is anti-conj-closed (n elements, n even), its characteristic
polynomial has the coefficient structure:

> For anti-conj-closed multiset of n elements: the polynomial's **even-indexed**
> coefficients are purely imaginary, **odd-indexed** coefficients are purely real.
> These n/2 + n/2 = **n** independent reals suffice to represent the polynomial.

**Algorithm (implemented):**
```python
w = z_reps ** 2                                        # n values, anti-conj-closed
coeffs, _, _ = zip_embed(w)
return coeffs.imag[0::2] + 1j * coeffs.real[1::2]     # n/2 complex → n reals
```

Output: **n reals** (vs. current 2n from NegPairEncoder).  No τ-half selection
is needed — `{z_k^2}` is embedded directly.  zip_embed is multiset-aware and
handles degenerate events (repeated w values) correctly.

---

## 6. Degenerate Events

σ-fixed points occur when `eval(u_k) = eval(v_j)`, giving `z_k = a(1+i)` and
`w_k = a(1−i)/√2` (real, for the TYPE_11 rotation).  In the TYPE_11 case,
`w_k.real` survives but the conjugate-pairing degenerates.

For TYPE_NEG: note that `τ(w) = −conj(w) = −a+bi` for `w = a+bi`.  The τ-pair
`{w, τ(w)}` has the **same imaginary part** but **opposite real parts**.  The
degenerate case is `Re(w) = 0` (purely imaginary w), where `τ(w) = w` (τ-fixed
point).  These are measure-zero.

Both degenerate cases are handled correctly by the implemented approach:
**TYPE_11** embeds `{w_k = z_k · e^{−iπ/4}}` directly via `zip_embed` (which is
multiset-aware); **TYPE_NEG** embeds `{z_k^2}` directly via `zip_embed`.
No manual half-selection is needed in either case.

---

## 7. Compression Summary

After σ-compression is implemented:

| Context | Current encoder | Current output | σ-encoder | σ-output |
|---|---|---|---|---|
| TYPE_11 self-pair ASSOC | Type11PairEncoder | 2n | SelfPairType11Encoder | n |
| TYPE_NEG self-pair ASSOC | NegPairEncoder | 2n | SelfPairNegEncoder | n |
| TYPE_12/21 ASSOC | Type12/21PairEncoder | 2n | (no new symmetry known) | — |
| TYPE_22 ASSOC | Type22PairEncoder | 2n | (no new symmetry known) | — |
| Non-self-pair TYPE_NEG ASSOC | NegPairEncoder | 2n | (no swap symmetry) | — |

The min(output_dim) selection in the overlap-block encoder naturally picks
the σ-encoders when they apply, because they produce smaller output_dim.

**Both TYPE_11 and TYPE_NEG self-pairing cases achieve category (b) compression.**
Existing factor-2 compressors (TYPE_12, TYPE_21, NEG for non-self-pairing blocks)
produce 2n reals for n base pairs.  σ-compression produces n reals — strictly
better — for both TYPE_11 and TYPE_NEG self-pairing.  The σ-step is the same
factor-of-2 halving in both cases (current 2n → n); the fact that NegPairEncoder
has NOTIONAL_FACTOR=2 is separate accounting of prior compression and does not
reduce the σ-step's benefit.

---

## 8. Implementation Approach

The σ-encoders should be implemented as new row-pair encoder classes and
factories, following the same pattern as `Type11PairEncoderFactory`,
`NegPairEncoderFactory`, etc. in `symcoder/encoders/row_pair_encoders.py`.

Each factory's `assess()` calls `_pair_orbit_atoms(pf, plan)`, applies the
existing partition helper to confirm the orbit type (TYPE_11 or TYPE_NEG), then
additionally checks `_is_self_pairing_block(pf)` (same operation, same
flavour).  If both conditions hold, it returns the σ-compressed encoder with
smaller `output_dim`; otherwise returns `[]` so the standard encoder is used
as fallback.

No changes to the overlap-block or phase-2 machinery are needed; `min(output_dim)`
selection handles everything automatically.

---

## 9. Open Questions

1. **Form-1 vs form-2 for TYPE_NEG self-pairing.** Does `_try_neg_partition`
   always return form-1 reps `(+u, +v)` for self-pairing blocks?  If so, the
   σ-closure is guaranteed and Section 5 is valid.  If form-2 reps can occur,
   the argument needs revisiting.  This is the most pressing question before
   implementation.

2. **Does TYPE_22 appear in self-pairing blocks?**  In tested contexts (small
   single group) all ANTISYMMETRIC self-pairing orbits are TYPE_NEG.  For
   larger groups with zero-overlap self-pairs, TYPE_22 may arise.  If so, a
   third σ-encoder case would be needed (its structure is not yet analysed).

3. **Are there further symmetries within the n/2 representatives?**  For the
   TYPE_11 case, the n/2 representative values `{w_k : Im(w_k) > 0}` are
   embedded as a generic complex multiset.  No forced further structure has been
   identified.  This is conjectured to be the floor.

4. **Non-self-pairing cross-ASSOC segments.**  No swap symmetry, so no
   σ-closure.  Whether other functional equations (from the group-theoretic
   structure of the overlap partition) give further compression is an open
   question with no current candidate mechanism.
