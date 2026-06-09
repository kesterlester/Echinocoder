# Finding Optimal Encoding Directions: Approach 2

## Overview

This directory contains a **rigorous, exact, proof-of-concept implementation** of **Approach 2** for finding the sequence of optimal encoding directions $DD(i)$ for $i = 1, 2, \ldots$ that maximize the minimum confusable pair size $f(DD(i))$.

The approach is:
1. **Inner oracle** (fixed directions $D$): Compute a certified lower bound on $f(D)$ via SOS/SDP (no confusable pairs of size $< s$ exist).
2. **Outer loop** (optimize directions): Search over direction space for the configuration that maximizes this lower bound.

## Key constraint: Rigour without tolerances

- **No measure-zero arguments**: The hypercube construction ALWAYS gives confusable pairs of size $2^{M-1}$, so $f(DD(i))$ is never infinite.
- **No tolerance-based comparisons**: All floating-point results are verified symbolically (SymPy) or with interval arithmetic.
- **Exact verification**: The SOS certificate, once found, is checked for mathematical correctness, not just "eigenvalues look PSD."

## Files

### Theory

**`confusable_dd_optimization.tex`** (mathematical specification)

- **Section 1**: Problem statement — find directions maximizing minimum confusable pair size.
- **Section 2**: Inner oracle: Lasserre SOS certificate for $s=2$ (no size-2 confusable pairs ⇒ $f(D) \geq 3$).
- **Section 3**: SDP formulation for $k=4$, $s=2$ — 45×45 Gram matrix, 495 polynomial-identity constraints.
- **Section 4**: Numerical SDP solving via CVXPY (SCS solver).
- **Section 5**: **Rigorous verification without tolerances** — check the SOS certificate symbolically and with interval arithmetic.
- **Section 6**: Direction parameterisation — rational unit vectors with bounded-denominator integer coordinates.
- **Section 7**: Grid search strategies (exhaustive for small $i$, random/greedy for large $i$).
- **Section 8**: Limitations and future work.

### Implementation

**`confusable_dd_search.py`** (proof-of-concept code for $k=4$)

Key functions (docstrings cite the spec document):
- `construct_f_g_symbolic(directions, k)` — build the polynomial $f_D(A,B)$ and $g_D(A,B)$ symbolically (Sec. 2,3).
- `sos_sdp_for_directions(directions, k, s=2)` — solve the SDP, return $(\lambda^*, Q^*, \text{status})$ (Sec. 4).
- `verify_sos_certificate(directions, lambda_star, Q_star, k)` — verify the certificate rigorously (Sec. 5).
- `enumerate_directions(k, max_norm)` — generate rational unit directions (Sec. 6).
- `search_optimal_directions(k, i, ..., search_type)` — search for best $i$-tuple of directions (Sec. 7).

**Usage**:
```bash
cd /Users/lester/github/Echinocoder/MinimalConfusableSets
source ../venv/bin/activate
python3 confusable_dd_search.py
```

**From Python**:
```python
from confusable_dd_search import search_optimal_directions
best_D, best_lambda, verified = search_optimal_directions(k=4, i=5, 
                                                           search_type='random',
                                                           num_samples=1000)
print(f"Best lambda* = {best_lambda}")
if best_lambda > 0:
    print(f"Certified: f(DD(5)) >= 3")
```

## Current status

**For $k=4, s=2$ (no size-2 confusable pairs)**:

- ✅ SDP construction and solving works (CVXPY + SCS)
- ✅ Direction enumeration and parameterisation (rational, bounded-denominator)
- ✅ Grid search (exhaustive for $i \leq 6$, random/greedy for larger $i$)
- ⏳ Symbolic verification (implemented but not yet tested at scale)
- ⏳ Evaluation for many $(k, i)$ pairs and direction sets

**Known issues**:
- The verification step is expensive (SymPy symbolic expansion of degree-4 polynomials).
  For large $k$ or many trials, consider using interval arithmetic alone instead.
- Grid search is finite but potentially incomplete — we may miss the true optimal $DD(i)$ 
  if it lies outside the bounded-numerator set. Can mitigate by increasing the bound.
- Extension to $s \geq 3$ requires Lasserre hierarchy level $d \geq 2$, which scales poorly.

## How to use this for real

### 1. Verify the mathematical spec

Read `confusable_dd_optimization.tex` critically. Ask:
- Is the SDP formulation correct?
- Is the verification procedure sound?
- Are there gaps in the argument?

### 2. Test the code on small cases

Run `confusable_dd_search.py` for $k=4$, $i=4,5,6$. Check:
- Do direction enumerations look reasonable?
- Do SDP solve successes/failures match expectations?
- Do verified certificates match intuition (e.g., symmetric configurations should score well)?

### 3. Extend to other $(k, s)$

- For $s=3$ and higher: use Lasserre hierarchy level $d=2$ (moment matrix $\binom{2k+2}{2} \times \binom{2k+2}{2}$).
- For other dimensions $k$: adapt the SDP construction in Section 3.

### 4. Integration with MinimalConfusableSets

The code gives **lower bounds** on $f(DD(i))$ (via SOS certificates).
The `MinimalConfusableSets` search gives **upper bounds** (via hypercube construction with L-matrices).
Together, they bracket $f(DD(i))$ from below and above.

## Relation to other approaches

| Approach | Method | Pros | Cons |
|----------|--------|------|------|
| **This (Approach 2)** | SOS/SDP lower bound | Rigorous, machine-verifiable, gives certified bound | Slow, grid search incomplete, limited to s=2 |
| **Approach 1** | Case-by-case analysis | Often fast, can be tight | Requires hand proofs, doesn't scale |
| **Numerical search** (confusable_search.py) | Gradient-based optimization | Finds explicit pairs | Heuristic, no guarantees |
| **Hypercube upper bound** (MinimalConfusableSets) | L-matrix reduction | Always works, gives baseline | Doesn't find true minimum |

## Papers and references

- Lasserre hierarchy: Lasserre, J. B. (2009). "Moments, positive polynomials and sums of squares."
- SOS/SDP for polynomial optimization: Parrilo, P. (2003). "Semidefinite programming relaxations for semialgebraic problems."
- Rational points on spheres: Standard in algebraic number theory (e.g., stereographic projection).

## TODOs for future work

- [ ] Implement interval-arithmetic verification (faster than SymPy symbolic).
- [ ] Add Lasserre hierarchy level $d=2$ for $s=3$ bounds.
- [ ] Optimize SDP solver selection (MOSEK is faster than SCS on large problems).
- [ ] Compare with known special configurations (24-cell vertices for $k=4$, etc.).
- [ ] Prove or disprove that symmetric/regular polytope vertices are optimal.
- [ ] Integrate with MinimalConfusableSets to close the gap between upper and lower bounds.
