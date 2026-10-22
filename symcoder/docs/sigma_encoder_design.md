# σ-Encoder Design Notes (scratch)

## Key resolved questions

### Form-1 vs form-2 for TYPE_NEG: RESOLVED — both forms arise; encoding works for both

The earlier proof that TYPE_NEG self-pairing is always form-1 was incorrect.
Form-2 groups {(+u,−v),(−u,+v)} do arise empirically in self-pairing blocks
(confirmed by test_neg_self_pairing_reps_u_sign_is_positive).

The encoding is correct for both forms because {z_k^2} is τ-closed regardless:
- Form-1 rep (+u,+v): z = a+ib.  Swap gives z' = b+ia.  z'^2 = −conj(z^2) = τ(z^2) ✓
- Form-2 rep (+u,−v): z = a−ib.  Swap gives z' = b−ia.  z'^2 = −conj(z^2) = τ(z^2) ✓

`_try_neg_partition` selects the element with u.sign=+1.  For form-1 this gives
v.sign=+1; for form-2 it gives v.sign=−1.  SelfPairNegEncoder handles both
identically by embedding z^2 directly (no form check needed).

### Self-pairing orbit types: RESOLVED — TYPE_11, TYPE_NEG, TYPE_22 only

The achievable sign set must be closed under the swap (s_u,s_v) → (s_v,s_u)
(because u and v roles are interchangeable in a self-pairing block).
TYPE_12 {(+,+),(+,−)} and TYPE_21 {(+,+),(−,+)} are not swap-closed → excluded.
TYPE_11, TYPE_NEG, TYPE_22 are all swap-closed → possible.

---

## Encoder design

### New helper

```python
def _is_self_pairing_block(pf) -> bool:
    """True if op_u==op_v and F_u==F_v (same operation and flavour, any overlap).
    Note: _is_self_pair() implies this but also requires full overlap (NULL_SELF case)."""
    return (pf.op_u == pf.op_v and
            tuple(pf.flavour_u.counts) == tuple(pf.flavour_v.counts))
```

### SelfPairType11Encoder

Orbit: TYPE_11 self-pairing.  All n reps have sign (+,+).
σ-closure: {z_k, i·conj(z_k)} both present among the n z-values.

Encode:
```python
z = np.array([complex(evaluate(u,e), evaluate(v,e)) for u,v in self._reps])
rotation = np.exp(-1j * np.pi / 4)   # (1−i)/√2
w = z * rotation                      # n conjugate-closed values
coeffs, _, _ = zip_embed(w)           # real polynomial (conj-paired roots)
return coeffs.real                    # n reals
```

output_dim = n (vs 2n for Type11PairEncoder).
_NOTIONAL_FACTOR = 2 (notional = 2n, which is what TYPE_11 without σ stores).
_METHOD_NAME = "11_sigma"

Factory assess():
- if _is_self_pair(pf): return []   (NULL_SELF case)
- if not _is_self_pairing_block(pf): return []
- pairs = _pair_orbit_atoms(pf, plan)
- if not all(u.sign==1 and v.sign==1 for u,v in pairs): return []
  (TYPE_11 means all elements have sign (+,+); if any differ it's not TYPE_11)
- return [SelfPairType11Encoder(pairs, pf, plan)]

Note: all orbit elements are the reps for TYPE_11 (orbit_size = n = len(pairs)).

### SelfPairNegEncoder

Orbit: TYPE_NEG self-pairing.  n reps (u.sign==+1; v.sign may be ±1).
Both form-1 {(+u,+v),(−u,−v)} and form-2 {(+u,−v),(−u,+v)} arise.
In either case {z_k²} is anti-conj-closed (τ-closed) by the swap argument:
the swap partner (+v,−u) or (+v,+u) contributes z'^2 = −conj(z²) = τ(z²).

Encode:
```python
z = np.array([complex(evaluate(u,e), evaluate(v,e)) for u,v in self._reps])
w = z ** 2                            # n values, anti-conj-closed (τ-closed)
coeffs, _, _ = zip_embed(w)
# Anti-conj-closed poly of n elements: even-indexed coeffs purely imaginary,
# odd-indexed coeffs purely real → exactly n independent reals
return _complex_to_reals(coeffs.imag[0::2] + 1j * coeffs.real[1::2])   # n reals
```

output_dim = n (vs 2n for NegPairEncoder).
_NOTIONAL_FACTOR = 4 (notional = 4n; TYPE_NEG orbit size = 2n, TYPE_11-equivalent
would need 4n reals).
_METHOD_NAME = "neg_sigma"

Factory assess():
- if _is_self_pair(pf): return []
- if not _is_self_pairing_block(pf): return []
- reps = _try_neg_partition(_pair_orbit_atoms(pf, plan))
- if reps is None: return []
- return [SelfPairNegEncoder(reps, pf, plan)]

### Registration order update

standard_row_pair_factories() should list σ-factories before their non-σ counterparts.
Since σ-encoders have output_dim=n < 2n (non-σ), min(output_dim) picks σ automatically
regardless of order — but listing them first is clearer:

```python
return [
    SelfPairEncoderFactory(),         # output 0, NULL_SELF
    SelfPairNegEncoderFactory(),       # output n, NEG self-pairing (σ-compressed)
    SelfPairType11EncoderFactory(),    # output n, TYPE_11 self-pairing (σ-compressed)
    NegPairEncoderFactory(),           # output 2n, TYPE_NEG (non-self-pairing)
    Type22PairEncoderFactory(),        # output 2n, TYPE_22
    Type21PairEncoderFactory(),        # output 2n, TYPE_21
    Type12PairEncoderFactory(),        # output 2n, TYPE_12
    Type11PairEncoderFactory(),        # output 2n, TYPE_11 (universal fallback)
]
```

---

## Tests needed

1. test_neg_self_pairing_reps_are_form1
   For all TYPE_NEG pairs in a self-pairing block: v.sign == +1 for all reps.

2. test_sigma_type11_output_dim_is_n
   SelfPairType11Encoder.output_dim == len(reps) == n.

3. test_sigma_neg_output_dim_is_n
   SelfPairNegEncoder.output_dim == len(reps) == n.

4. test_sigma_type11_permutation_invariant
   encode(plan, event) == encode(plan, permuted_event) for a dot-only plan
   using a context where self-pairing dot×dot associations exist.

5. test_sigma_neg_permutation_invariant
   encode(plan, event) == encode(plan, permuted_event) for an eps3-only plan.

6. test_sigma_type11_injectivity (the "bold claim" verification)
   Two different random events give different encodings.  Repeated 100 times.
   If σ-compression were lossy, this would occasionally give equal encodings.

7. test_sigma_neg_injectivity (same for NEG self-pairing)

8. test_sigma_output_dim_smaller_than_non_sigma
   For a plan that has self-pairing ASSOCs, total encode() length with
   standard_row_pair_factories() (including σ-factories) is less than with
   a factory list that omits the σ-factories.

---

## Open: TYPE_22 self-pairing

For large-enough groups, self-pairing ANTISYMMETRIC blocks may produce TYPE_22
orbits.  The σ-compression argument has not been developed for TYPE_22 self-pairing.
No SelfPairType22Encoder is implemented here; Type22PairEncoderFactory handles
these as normal (output 2n).

---

## Implementation status

[x] _is_self_pairing_block helper
[x] SelfPairType11Encoder / SelfPairType11EncoderFactory
[x] SelfPairNegEncoder / SelfPairNegEncoderFactory
[x] standard_row_pair_factories() update
[x] Tests (symcoder/tests/test_sigma_encoders.py — 8 tests, 206 parametrised cases)
[x] sigma_symmetry_analysis.md: updated Sections 1, 5, 6 with corrected analysis
[x] Update project memory
