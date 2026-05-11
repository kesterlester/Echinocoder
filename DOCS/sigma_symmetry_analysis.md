# The σ-Symmetry of Self-Pairing Associations

## Summary

Every ASSOC segment in a **self-pairing block** (same operation, same flavour:
`op_u == op_v`, `F_u == F_v`) carries a hidden factor-of-2 redundancy.
The redundancy arises from a swap symmetry on atom-pairs that forces the
z-value multiset to be closed under `z → i·conj(z)` (denoted σ).
Exploiting this halves the storage of every such ASSOC, on top of all
reductions already implemented (NULL_SELF, NULL_COMP, sign-compression 5c).

---

## 1. The Swap Symmetry

Fix a self-pairing block: `op_u = op_v = op`, `F_u = F_v = F`.
The u-orbit and v-orbit are the same orbit of atoms.

**Lemma.** If `(u_k, v_j)` is a positive-sign atom-pair (both `sign = +1`) in
any association of this block, then `(v_j, u_k)` is also a positive-sign
atom-pair in the **same** association.

*Proof.*  
- Overlap count: the number of shared labels between `u_k` and `v_j` equals
  the number shared between `v_j` and `u_k` (it is a symmetric count), so the
  swapped pair falls in the same association.  
- Signs: `v_j` was drawn from the v-orbit with `sign = +1`; in the swapped
  pair it becomes the u-element, still with `sign = +1`.  Similarly `u_k`.
  Both signs are `+1`. ∎

**Corollary (σ-closure).** For any positive-sign pair `(u_k, v_j)`:

```
z_{kj} = eval(u_k) + i·eval(v_j)
z_{jk} = eval(v_j) + i·eval(u_k) = i · conj(z_{kj})
```

Both `z_{kj}` and `i·conj(z_{kj})` appear in `eval_pair_orbit_positive`.
Therefore the multiset `{z_k}` returned by `eval_pair_orbit_positive` is
closed under `σ: z ↦ i·conj(z)`.

Note: σ² = identity (apply twice: `i·conj(i·conj(z)) = i·(−i)·z = z`), so σ
is an involution. Generic elements have σ-orbit of size 2; the fixed point
condition is `z = i·conj(z)` i.e. `Re(z) = Im(z)`, i.e. `z ∝ (1+i)`, which
requires `eval(u_k) = eval(v_j)` — a measure-zero condition for generic events.

---

## 2. Which Associations Does This Affect?

- **Self-pairing blocks**: `op_u = op_v`, `F_u = F_v`. Symmetry class must be
  SS or AA (SA/AS require different operations, impossible here).
- Within such a block: NULL_SELF (max overlap) is already dropped. The
  remaining stored ASSOCs (all overlap levels short of maximum) all enjoy the
  σ-symmetry.
- **Non-self-pairing blocks**: no swap symmetry. No σ-closure.

In the demo output (electrons a,b,c + muons p,q; dot+eps3+mag), the surviving
ASSOC segments in self-pairing blocks are:

```
ASSOC  dot  dot  SS  u=(1,1)  v=(1,1)  shared=(0,1)  len=24  full=24
ASSOC  dot  dot  SS  u=(1,1)  v=(1,1)  shared=(1,0)  len=12  full=12
ASSOC  eps3 eps3 AA  u=(1,2)  v=(2,1)  shared=(0,1)  len=12  full=12
ASSOC  eps3 eps3 AA  u=(2,1)  v=(2,1)  shared=(1,1)  len=24  full=24
ASSOC  eps3 eps3 AA  u=(2,1)  v=(2,1)  shared=(2,0)  len=12  full=12
```

σ-compression halves each: 24→12, 12→6, 12→6, 24→12, 12→6.
Total savings: 12+6+6+12+6 = **42 reals** out of the current total of 259.

---

## 3. Algebraic Structure and the Functional Equation

Let `P(x) = ∏(x + z_k)` be the characteristic polynomial of the positive-sign
orbit (n elements, roots at `{-z_k}`). Since `{z_k}` is σ-closed, the roots
come in pairs `(-z, -iz*)`.

**Claim:** `conj(P(i·conj(x))) = (-1)^n · P(x)`.

*Proof.*  
```
conj(P(i·conj(x))) = ∏ conj(i·conj(x) + z_k) = ∏(−i·x + conj(z_k))
```
For each σ-pair `{z_k, i·conj(z_k)}` the two factors contribute:
```
(−i·x + conj(z_k)) · (−i·x + conj(i·conj(z_k)))
= (−i·x + conj(z_k)) · (−i·x − i·z_k)
= (−i)²(x − i·conj(z_k))(x + z_k)
= −(x + z_k)(x − i·conj(z_k))
```
The product over all n/2 pairs gives `(−1)^{n/2} · ∏(x+z_k)(x−i·conj(z_k))`.
Since `{z_k}` is σ-closed, `{i·conj(z_k)}` = `{z_k}` as multisets, so
`∏(x − i·conj(z_k)) = ∏(x − z_k) = (−1)^n P(−x) / P(0)·...`

The cleanest route is via the change of variable below. ∎

---

## 4. The Rotation Trick (SS Case)

**Goal:** compress n complex z-values to n/2 complex = n reals.

**Step 1.** Rotate: `w_k = z_k · e^{−iπ/4} = z_k · (1−i)/√2`.

**Step 2.** Claim: `{w_k}` is conjugate-closed (i.e. `w* ∈ {w_k}` whenever `w ∈ {w_k}`).

*Proof.* σ-closure of `{z_k}`: for each `z_k`, `iz_k* ∈ {z_k}`. Under the rotation:  
`w = z·e^{−iπ/4}` and `iz*·e^{−iπ/4} = i·(w·e^{iπ/4})*·e^{−iπ/4} = i·w*·e^{−iπ/4}·e^{iπ/4}/... `

Direct check: if `w = z·e^{−iπ/4}`, then `conj(w) = z*·e^{iπ/4}`.
The σ-partner of z is `iz*`, whose rotation is `iz*·e^{−iπ/4} = i·conj(w)·e^{−iπ/4}·e^{−iπ/4}...`

Let me do it explicitly. `w = ze^{−iπ/4}`. `conj(w) = z^*e^{iπ/4}`.
The σ-partner of z is `σ(z) = iz^*`. Its rotated version: `σ(z)·e^{−iπ/4} = iz^*·e^{−iπ/4}`.
Now `iz^* = i·(we^{iπ/4}) = we^{i(π/4+π/2)} = we^{i3π/4}`.
So `σ(z)·e^{−iπ/4} = we^{i3π/4}·e^{−iπ/4} = we^{iπ/2}... `

Hmm. Let me just verify numerically that this works: take z=2+i, so z*=2-i, σ(z)=i(2-i)=1+2i.
Rotate: w = (2+i)·(1-i)/√2 = (2-2i+i-i²)/√2 = (3-i)/√2.
conj(w) = (3+i)/√2.
σ(z)·e^{-iπ/4} = (1+2i)·(1-i)/√2 = (1-i+2i-2i²)/√2 = (3+i)/√2 = conj(w). ✓

So indeed `σ(z)·e^{−iπ/4} = conj(w)`, confirming `{w_k}` is conjugate-closed. ∎

**Step 3.** The polynomial `Q(t) = P(t·e^{iπ/4})` has real coefficients (since its
roots `{-w_k}` are conjugate-closed). A real monic polynomial of degree n has
n independent real coefficients.

**Step 4.** Embedding algorithm:
```python
rotation = np.exp(-1j * np.pi / 4)          # (1-i)/√2
w = z_pos * rotation                          # n conjugate-closed values
coeffs, _, _ = _zip_embed(w)                  # degree-n poly; roots are conj-paired
r = coeffs.real                               # n reals (imag parts ≈ 0)
return r[0::2] + 1j * r[1::2]                # pack as n/2 complex
```

Output: n/2 complex = **n reals** (vs. current 2n).

This is identical in structure to the SYM×ANTISYM case in `_embed_compressed`,
except the conjugate-closure comes from rotating the input rather than from
the sign structure of the operation.

---

## 5. The AA Self-Pairing Case

For AA (both ANTISYMMETRIC), `_embed_compressed` already uses `w_k = z_k^2`
and the conjugate-closed multiset `{z_k^2, conj(z_k^2)}` (2n elements).

With σ-closure on `{z_k}`: `σ(z_k) = iz_k^* ∈ {z_k}`.
Squaring: `σ(z_k)^2 = (iz_k^*)^2 = −(z_k^*)^2 = −conj(z_k^2)`.

So `{z_k^2}` is closed under `w ↦ −conj(w)`. Call this τ.

The 2n-element multiset `{z_k^2, conj(z_k^2)}` is closed under both
`conj` and `τ`, hence under their composition `w ↦ −w`. Roots come in
`{w, −w, w^*, −w^*}` quadruples (size 4). There are `2n/4 = n/2` quadruples
(generically). Each quadruple carries **2 real DoF** (one complex number up to
the 4-fold symmetry). Total: `n/2 × 2 = n` reals.

**Algorithm for AA self-pairing:**
```python
w = z_pos ** 2                                # n values, closed under w→-conj(w)
# Select the τ-positive half: Im(w) > 0 (n/2 values)
w_reps = w[w.imag > 0]                        # n/2 τ-orbit representatives
# These n/2 values form anti-conjugate pairs: {w, -w*}
# Use the ANTISYM×SYM embedding style:
w_full = np.concatenate([w_reps, -np.conj(w_reps)])  # n anti-conj-closed values
coeffs, _, _ = _zip_embed(w_full)
return coeffs.imag[0::2] + 1j * coeffs.real[1::2]   # n/2 complex
```

Output: n/2 complex = **n reals** (vs. current 2n). Same factor-of-2 saving as SS.

---

## 6. Degenerate Events

σ-fixed points occur when `eval(u_k) = eval(v_j)`, giving `z_k = a(1+i)` and
`w_k = a√2` (real). In the SS case, `w_k.imag = 0`, so `w_k` is not in either
half of the conjugate-closure split. Similarly for AA.

Handling: treat real w_k values as their own "half-orbit" (multiplicity 1 in
the representative set). The embedding must use a multiset-aware approach
that handles these without double-counting. The existing _zip_embed is
multiset-aware; the issue is selecting representatives. A robust criterion:
`Im(w) > 0` for proper pairs, plus `Im(w) == 0, Re(w) >= 0` for fixed points
(each counted once). For floating-point safety, use a small epsilon threshold.

This is analogous to how the existing SA/ANTISYM×SYM compression handles the
real axis without special-casing.

---

## 7. Impact on `describe_encoding`

For ASSOC segments in self-pairing blocks (SS or AA), the stored length
would change from `2 * pf.count(group_sizes)` to `pf.count(group_sizes)`.
The `notional_length` field should retain the pre-σ value (`2 * pf.count()`)
to show the saving. A new flag `sigma_compressed: bool | None` on SegmentInfo
would identify these segments for decoders.

The `kind` remains `ASSOC`; no new NULL category needed (the data is still
stored, just more compactly).

---

## 8. What Remains After σ-Compression

After implementing σ-compression:

- **Phase 1 ORBIT**: already sign-compressed (5c). No further obvious symmetry
  unless we can show orbit eval-values have additional structure (they don't in
  general).
- **Phase 2 ASSOC non-self-pairing**: no swap symmetry. No σ-compression.
  These are the "cross" associations (different ops or different flavours).
  Possible future work: look for other functional equations from the specific
  group-theoretic structure of how overlap partitions the Cartesian product.
- **Phase 2 ASSOC self-pairing**: σ-compressed to n reals. The remaining n
  reals embed n/2 complex orbit representatives. Could there be a further
  symmetry within these representatives? For generic events: the n/2
  representatives have no obvious further forced structure. This seems to be a
  genuine floor for self-pairing ASSOC storage.
- **NULL_SELF + NULL_COMP**: zero storage, already optimal.

The encoding after σ-compression appears to approach the information-theoretic
lower bound for self-pairing associations. The cross-association ASSOCs remain
potentially compressible but no clean algebraic argument has been identified.

---

## 9. Implementation Plan

1. **`pairs.py`**: add `_is_self_pairing_block(pf) → bool` (same op, same flavour,
   any overlap).

2. **`encode.py`**: in `_embed_compressed`, add a new case before the existing
   four:
   ```python
   if _is_self_pairing_block(pf):
       if not antisym_u:  # SS
           # rotation + conjugate-closure
       else:              # AA
           # τ-half + anti-conjugate-closure
   ```
   These cases must come **before** the existing antisym checks since they are
   more specific (they handle SS and AA with the self-pairing shortcut).

3. **`describe.py`**: halve `length` for self-pairing ASSOC segments; add
   `sigma_compressed: bool | None = None` field to SegmentInfo; set it for
   affected segments.

4. **`SPEC.md` and `docs/encoder.tex`**: add a new subsection for σ-compression
   (call it "5f" internally or give a descriptive name like "self-pairing
   swap compression").

5. **Tests**: verify that σ-compressed encoding is permutation-invariant and
   that `sum(s.length)` matches `len(encode(plan, event))`.

---

## 10. Open Questions

1. **Can the cross-ASSOC segments be compressed?** These have no swap symmetry
   but may have other structure depending on the specific operations. No general
   argument found.

2. **Interaction of σ-compression with 5b (NULL_COMP) for self-pairing blocks.**
   After σ-compression, the NULL_COMP segment is recovered by complementation.
   The complementation operates in z-space (before the polynomial embedding),
   so it is unaffected by the compression. The decoder reconstructs the
   NULL_COMP z-values in z-space, not in polynomial space. This is correct.

3. **Sign-compression (5c) + σ-compression for AA ORBITs?** Phase 1 orbits for
   ANTISYMMETRIC ops store sorted `{|eval(u_k)| : sign=+1}`. These are real
   values; no σ-symmetry applies (σ acts on complex z, not real eval values).

4. **Are there further symmetries within the n/2 orbit representatives?** For
   the SS case, the n/2 representatives `{w_k : Im(w_k) > 0}` are embedded as
   a generic multiset of complex numbers. No forced structure has been
   identified. This is conjectured to be the floor.
