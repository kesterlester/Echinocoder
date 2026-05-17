# The σ-Symmetry of Self-Pairing Associations

## Summary

Every ASSOC segment in a **self-pairing block** (same operation, same flavour:
`op_u == op_v`, `F_u == F_v`) carries a hidden factor-of-2 redundancy.
The redundancy arises from a swap symmetry on atom-pairs that forces the
z-value multiset to be closed under `z → i·conj(z)` (denoted σ).
Exploiting this halves the storage of every such ASSOC, on top of all
reductions already implemented (NULL_SELF, NULL_COMP, half-sort compression).

Self-pairing blocks contain only TYPE_11 orbits (SYMMETRIC operation) or
TYPE_NEG orbits (ANTISYMMETRIC operation); TYPE_12/21 cannot arise because
both operations are identical.

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

> **⚠ SUSPECT — needs checking:** The lemma assumes the representatives have
> both atoms with `sign = +1`.  For TYPE_11 self-pairing this is automatic
> (SYMMETRIC operations carry no sign degree of freedom, so every orbit element
> has sign +1).  For TYPE_NEG self-pairing the physical partition assigns reps
> with `u.sign = +1`; the v-sign may be `+1` (form-1 orbit: `{(+u,+v),(−u,−v)}`)
> or `−1` (form-2 orbit: `{(+u,−v),(−u,+v)}`).  The σ-closure argument goes
> through cleanly for form-1, but for form-2 the swapped representative
> `(v_j, u_k)` has `v_j.sign = +1` and `u_k.sign = +1`, giving
> `z_{jk} = eval(v_j) + i·eval(u_k)`, which is NOT `i·conj(z_{kj})` when
> `z_{kj} = eval(u_k) + i·eval(−v_j)`.  It must be verified empirically
> whether TYPE_NEG self-pairing orbits always produce form-1 reps.

---

## 2. Which Associations Does This Affect?

- **Self-pairing blocks**: `op_u = op_v`, `F_u = F_v`. The orbit types that
  can arise are TYPE_11 (SYMMETRIC operation) and TYPE_NEG (ANTISYMMETRIC
  operation, physical partition confirms this for the tested contexts).
  TYPE_12/21 cannot arise because both operations are identical.
- Within such a block: NULL_SELF (max overlap) is already dropped. The
  remaining stored ASSOCs (all overlap levels short of maximum) all enjoy the
  σ-symmetry, subject to the form-1 caveat above.
- **Non-self-pairing blocks**: no swap symmetry. No σ-closure.

> **⚠ SUSPECT — needs checking:** For ANTISYMMETRIC self-pairing with zero
> overlap between u and v (distinct label sets), TYPE_22 may appear instead of
> TYPE_NEG when the group is large enough.  In tested contexts (single group of
> 4 electrons) the only ANTISYMMETRIC self-pairing pairs have overlap=(2,) and
> are confirmed TYPE_NEG.  The claim "ANTISYMMETRIC self-pairing always gives
> TYPE_NEG" has NOT been verified for large groups or zero-overlap cases.

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

> **⚠ This section assumes form-1 reps throughout (see Section 1 caveat).**

The NegPairEncoder (for both self-pairing and non-self-pairing TYPE_NEG orbits)
already forms `w_k = z_k²` and embeds the n values directly as n complex
coefficients → 2n reals.

With σ-closure on `{z_k}` (form-1 reps): `σ(z_k) = iz_k^* ∈ {z_k}`.
Squaring: `σ(z_k)^2 = (iz_k^*)^2 = −(z_k^*)^2 = −conj(z_k^2)`.

So `{w_k = z_k^2}` is closed under `τ: w ↦ −conj(w)`.

The n values `{w_k}` form τ-pairs `{w, −w^*}` of size 2 (generically).
There are n/2 such pairs, each carrying **2 real DoF** (one complex number up
to the 2-fold τ-symmetry).  Total: n/2 × 2 = **n reals**.

**Algorithm:**
```python
w = z_reps ** 2                                        # n values, τ-closed
w_reps = w[w.imag > 0]                                 # n/2 τ-representatives
w_full = np.concatenate([w_reps, -np.conj(w_reps)])    # n anti-conj-closed values
coeffs, _, _ = zip_embed(w_full)
return coeffs.imag[0::2] + 1j * coeffs.real[1::2]     # n/2 complex → n reals
```

Output: **n reals** (vs. current 2n from NegPairEncoder).  This is a
genuine additional factor-of-2 compression on top of the existing TYPE_NEG
compression, applying only to self-pairing blocks.

> **⚠ SUSPECT — floating-point robustness of τ-half selection.**  The filter
> `w.imag > 0` to select τ-representatives must handle τ-fixed points
> (`w.imag == 0`).  See Section 6 on degenerate events.

---

## 6. Degenerate Events

σ-fixed points occur when `eval(u_k) = eval(v_j)`, giving `z_k = a(1+i)` and
`w_k = a(1−i)/√2` (real, for the TYPE_11 rotation).  In the TYPE_11 case,
`w_k.real` survives but the conjugate-pairing degenerates.

τ-fixed points for TYPE_NEG: `w = −w^*` iff `Re(w) = 0`, i.e. `w` is purely
imaginary.

In both cases the degenerate element must be counted once (not twice) among the
representatives.  A robust criterion: `Im(w) > ε` for proper pairs, plus
`Im(w) ≈ 0` elements counted once.  The existing `zip_embed` function is
multiset-aware; the implementation issue is correctly selecting the
representative half.

> **⚠ SUSPECT — the handling of degenerate cases has not been implemented or
> tested.  Correctness at measure-zero events must be verified.**

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
