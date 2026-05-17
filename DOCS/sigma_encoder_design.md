# σ-Encoder Design Notes (scratch)

## Key resolved questions

### Form-1 vs form-2 for TYPE_NEG: RESOLVED — always form-1

TYPE_NEG achievable set = {(+,+),(−,−)}.  A form-2 group {(+u,−v),(−u,+v)} has
per-group achievable contribution {(+,−),(−,+)}.  Combined with the global
{(+,+)}, this yields achievable = {(+,−),(−,+)}, which has both u_flips and
v_flips present → TYPE_22, not TYPE_NEG.

Therefore: `_try_neg_partition` only accepts groups where signs are (s,s) — i.e.
form-1 {(+u,+v),(−u,−v)}.  The σ-closure proof for TYPE_NEG holds without caveat.
The suspect callout in sigma_symmetry_analysis.md can be removed.

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

Orbit: TYPE_NEG self-pairing.  n reps, all form-1: v.sign==+1 guaranteed.
σ-closure: {z_k, i·conj(z_k)} both present.
τ-closure on {z_k²}: τ(w) = −conj(w), and Im(z²) + Im(−conj(z²)) = 0 so
the two τ-partners have opposite imaginary parts → Im > 0 selects exactly one.

Encode:
```python
z = np.array([complex(evaluate(u,e), evaluate(v,e)) for u,v in self._reps])
w = z ** 2                            # n values, τ-closed
w_half = w[w.imag > 0]               # n/2 τ-representatives (generic events)
# For degenerate events (Im(w)=0): τ-partner has same Im → take Re>0 ones too
# Full robustness TODO; assert for now:
assert len(w_half) == self._n // 2
w_full = np.concatenate([w_half, -np.conj(w_half)])   # n anti-conj-closed values
coeffs, _, _ = zip_embed(w_full)
# Anti-conj-closed poly: even-degree imag, odd-degree real
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

[ ] _is_self_pairing_block helper
[ ] SelfPairType11Encoder / SelfPairType11EncoderFactory
[ ] SelfPairNegEncoder / SelfPairNegEncoderFactory
[ ] standard_row_pair_factories() update
[ ] Tests (all 8 above)
[ ] sigma_symmetry_analysis.md: remove form-1 caveat, add proof
[ ] Update project memory
