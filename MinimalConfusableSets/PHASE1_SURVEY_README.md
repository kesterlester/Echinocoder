# Phase 1 Survey: Certified Lower Bounds on f(DD(i)) for k=4

## What this survey computes

For each $i \in [1,16]$ (number of encoding directions in $\mathbb{R}^4$), the survey reports:

1. **Theory-only results** ($i \in [1,6]$):
   - No search needed; proven bounds from the partition argument
   - $i=1,2,3$: $f(DD(i)) = 1$ (directions don't span $\mathbb{R}^4$)
   - $i=4,5,6$: $f(DD(i)) = 2$ (size-2 confusable pairs always exist via partition)

2. **Grid-search results** ($i \in [7,16]$):
   - Sample 500–800 random $i$-tuples of directions (bounded-numerator set)
   - For each, compute the $s=2$ SOS certificate: $\lambda^*(D) = \max \lambda$ s.t. $f(D) - \lambda g(D)$ is SOS
   - Report the **best** $\lambda^*$ found
   - **Certified lower bound**: 
     - If $\lambda^* > 0.001$: proven $f(DD(i)) \geq 3$
     - If $\lambda^* \approx 0$: inconclusive (theoretical bound $f(DD(i)) \geq 3$ still applies)

## Output

**CSV file**: `dd_phase1_results.csv` with columns:
- `i`: number of directions
- `M`: ceil(i/3) = number of bad-bat groups for k=4
- `Hypercube_UB`: upper bound from hypercube construction = $2^{M-1}$
- `Theory_LB`: proven lower bound from partition argument
- `Best_Lambda*`: best SOS certificate found (grid search)
- `Certified_LB`: lower bound we can prove with current code
- `Method`: 'theory_only' or 'grid_search'
- `Verified`: whether the SOS certificate was verified (expensive, skipped for $\lambda^* \approx 0$)
- `Notes`: details on search

**ASCII table**: Summary printed to stdout

## Interpretation

The survey answers: **What is the best we can prove about f(DD(i)) right now?**

Example for $i=7$:
```
i=7: 3 ≤ f(DD(7)) ≤ 4
  - Theory proves f ≥ 3 (no s=2 pairs via partition)
  - Hypercube gives f ≤ 4 (M=3 bad bats)
  - If best_λ* > 0: we also have f ≥ 3 (certified independently)
  - Gap: 1 unit — we don't know if f=3 or f=4
```

For $i=16$:
```
i=16: 3 ≤ f(DD(16)) ≤ 32
  - Gap: 29 units — very wide!
  - Would need s≥3 analysis (Lasserre) to tighten
```

## Computational cost

- **Theory only** ($i=1$–$6$): instant (no search)
- **Grid search** ($i=7$–$16$): 
  - ~500–800 SDP solves per $i$ × 10 values = ~6500 SDPs total
  - Each SDP is small (45×45 moment matrix for s=2)
  - Estimated time: 1–3 hours on modern hardware

## What Phase 1 does NOT do

- Does not search for optimal directions (that's Phase 2)
- Does not implement Lasserre SDP for s≥3 (that's Phase 2)
- Does not use specialized configurations (24-cell vertices, etc.)
- Does not integrate with MinimalConfusableSets hypercube code yet

## Next steps (Phase 2)

Once Phase 1 is done:
1. **Implement full Lasserre SDP** with power-sum constraints
2. **For $i \geq 7$**: try to certify $f(DD(i)) \geq 4$ (no size-3 pairs)
3. **Compare** with MinimalConfusableSets upper bounds from hypercube search

---

## How to run

```bash
cd /Users/lester/github/Echinocoder/MinimalConfusableSets
source ../venv/bin/activate
python3 dd_phase1_survey.py
```

Output goes to `dd_phase1_results.csv` and `dd_phase1_results.txt`.
