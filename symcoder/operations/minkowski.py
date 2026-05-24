"""
symcoder.operations.minkowski
=============================
Standard operations for 3+1-D Minkowski space with metric signature (+−−−).

All vectors are assumed to be 4-component numpy arrays ordered (t, x, y, z).

Operations
----------
dot     rank 2, scalar        v·w = v[0]w[0] − v[1]w[1] − v[2]w[2] − v[3]w[3]
mag_sq  rank 1, scalar        v² = v[0]² − v[1]² − v[2]² − v[3]²
                               (can be negative for spacelike vectors — no sqrt taken)
eps     rank 4, pseudoscalar  det([v1|v2|v3|v4]) — Levi-Civita contraction

Note on mag_sq
--------------
Unlike Euclidean magnitude, the Minkowski norm-squared can be positive
(timelike), negative (spacelike) or zero (lightlike).  Taking a square root
is left to the caller if the sign is known in context.
"""
import numpy as np
from symatom import Operation, ArgumentSymmetry

dot = Operation(
    "dot", rank=2, odd_parity=False,
    argument_symmetry=ArgumentSymmetry.SYMMETRIC,
    eval_fn=lambda v: float(
        v[0][0]*v[1][0] - v[0][1]*v[1][1] - v[0][2]*v[1][2] - v[0][3]*v[1][3]
    ),
    tex=r"\mathbf{#1} \cdot \mathbf{#2}",
)

mag_sq = Operation(
    "mag_sq", rank=1, odd_parity=False,
    argument_symmetry=ArgumentSymmetry.SYMMETRIC,
    eval_fn=lambda v: float(
        v[0][0]**2 - v[0][1]**2 - v[0][2]**2 - v[0][3]**2
    ),
    tex=r"\mathbf{#1}^2",
)

eps = Operation(
    "eps", rank=4, odd_parity=True,
    argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
    eval_fn=lambda v: float(np.linalg.det(np.column_stack(v))),
    tex=r"\varepsilon(\mathbf{#1},\mathbf{#2},\mathbf{#3},\mathbf{#4})",
)
