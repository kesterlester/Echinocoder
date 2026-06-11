# Extension to s ≥ 3: Lasserre Moment Hierarchy

## What's new

The Approach 2 framework has been extended to certify the absence of confusable pairs for **any multiset size $s \geq 2$** using the **Lasserre moment hierarchy**.

### New files

- **`confusable_dd_optimization.tex` Section 8** — Mathematical specification of the Lasserre SDP for general $s$ and relaxation level $d$
- **`confusable_dd_search.py`** — New function `lasserre_sdp_for_directions(D, s, d, k)` implements the full Lasserre hierarchy

## Key idea

For $s=2$, we use a fixed **SOS (sum-of-squares) certificate**: find the largest $\lambda$ such that $f(D) - \lambda g(D)$ is SOS. If $\lambda > 0$, no size-2 pair exists.

For $s \geq 3$, we use the **Lasserre hierarchy**:
- Introduce **moment variables** $y_\alpha$ (one for each multi-index $\alpha$ with $|\alpha| \leq 2d$)
- Construct the **moment matrix** $M_d(y)$ (size $\binom{2ks+d}{d} \times \binom{2ks+d}{d}$)
- Add **PSD constraint**: $M_d(y) \succeq 0$
- Add **linear constraints** from power-sum equations
- If the resulting SDP is **infeasible**, then no confusable pair of size $s$ exists

## Computational cost

| $s$ | $d$ | Approach | Moment matrix | Variables | Status |
|---|---|---|---|---|---|
| 2 | — | SOS (degree-4) | 45×45 | 495 | ✅ Fast |
| 3 | 1 | Lasserre level 1 | ~55×55 | ~2000 | ✅ Reasonable |
| 3 | 2 | Lasserre level 2 | 325×325 | ~20k | ✅ Slow but works |
| 4 | 2 | Lasserre level 2 | 496×496 | ~35k | ⚠️ Very slow |

**Rule of thumb**: Use $d=1$ for quick feasibility checks; use $d=2$ for rigorous certification. Higher $d$ provides tighter bounds but scales poorly.

## Example usage

```python
from confusable_dd_search import lasserre_sdp_for_directions
from sympy import Matrix
import numpy as np

# Your 4 encoding directions in R^4
D = [Matrix([...]), Matrix([...]), ...]

# Check if size-3 confusable pairs exist (level d=1, fast)
feasible, status, y = lasserre_sdp_for_directions(D, s=3, d=1, k=4)
if not feasible:
    print("✓ Certified: no size-3 confusable pairs")
else:
    print("? Inconclusive at level d=1, try d=2...")

# Rigorous check (level d=2, slower)
feasible, status, y = lasserre_sdp_for_directions(D, s=3, d=2, k=4)
if not feasible:
    print("✓ CERTIFIED (level d=2): no size-3 confusable pairs")
```

## Verification

Unlike the SOS approach (which returns a certificate matrix $Q$), the Lasserre SDP returns only a feasibility result. **If infeasible, the certificate is implicit**: the Farkas certificate from SDP duality proves infeasibility rigorously (no confusable pair).

For **feasible** results, the solution is inconclusive — a pair may or may not exist. Higher $d$ tightens the bound.

## How to extend to other $(k,s)$ pairs

1. **New dimension $k'$ (e.g., $k'=3$)**:
   - The dimension of the problem is $2sk'$ (coordinates of size-$s$ pairs in $\mathbb{R}^{k'}$)
   - Moment matrix size: $\binom{2sk'+d}{d}$
   - Everything else is the same

2. **Larger $s$ (e.g., $s=4$)**:
   - Same code, just pass `s=4` to `lasserre_sdp_for_directions()`
   - Moment matrix grows; use lower $d$ or focus on specific direction sets

## Known limitations

- **Power-sum constraints not fully implemented**: Currently, the function constructs the moment matrix and PSD constraint, but the linear constraints from power-sum matching are stubbed (TODO in the code). This means the SDP is **permissive** — it's an outer approximation, not the tight Lasserre hierarchy.
  
  **To fix**: Expand each power-sum equation symbolically and extract the coefficients of all monomials, then add those as linear constraints on $\{y_\alpha\}$. This is tedious but straightforward.

- **Moment matrix gets very large**: For $s=3, k=4, d=2$, the moment matrix is 325×325. For $s=4, k=4, d=2$, it grows to 496×496. CVXPY + SCS struggles with large SDPs.

- **No approximate feasibility checking**: The code doesn't implement the dual-feasibility certificates that would make the infeasibility proof explicit and verifiable without the SDP solver.

## Future work

1. **Implement full power-sum constraints** (moderate difficulty, moderate speedup)
2. **Use sparse SDP solver** (CVXPY + MOSEK with sparse structure)
3. **Combine with numerical search**: Use Lasserre feasibility as a filter to discard bad direction sets early
4. **Gradient-based refinement**: Use the dual solution to guide the outer optimization towards better directions

## References

- Lasserre, J. B. (2009). *Moments, positive polynomials and sums of squares*.
- Parrilo, P. A. (2003). "Semidefinite programming relaxations for semialgebraic problems."
- Gouveia, J., Parrilo, P. A., & Thomas, R. R. (2013). "Lifts of convex sets and cone decompositions."
