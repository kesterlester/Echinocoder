"""
symcoder.operations.euclidean3
==============================
Standard operations for 3-D Euclidean space.

All vectors are assumed to be 3-component numpy arrays.

Operations
----------
mag   rank 1, scalar     |v|
dot   rank 2, scalar     v·w
eps   rank 3, pseudoscalar   v·(w×x)  (scalar triple product / Levi-Civita)
"""
import numpy as np
from symatom import Operation, ArgumentSymmetry

mag = Operation(
    "mag", rank=1, odd_parity=False,
    argument_symmetry=ArgumentSymmetry.SYMMETRIC,
    eval_fn=lambda v: float(np.linalg.norm(v[0])),
    tex=r"|\mathbf{#1}|",
    mass_dimension=1,                # |k v| = k |v|
)

dot = Operation(
    "dot", rank=2, odd_parity=False,
    argument_symmetry=ArgumentSymmetry.SYMMETRIC,
    eval_fn=lambda v: float(np.dot(v[0], v[1])),
    tex=r"\mathbf{#1} \cdot \mathbf{#2}",
    mass_dimension=2,                # (k u)·(k v) = k^2 u·v
)

eps = Operation(
    "eps", rank=3, odd_parity=True,
    argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
    eval_fn=lambda v: float(np.dot(v[0], np.cross(v[1], v[2]))),
    tex=r"\varepsilon(\mathbf{#1},\mathbf{#2},\mathbf{#3})",
    mass_dimension=3,                # (k u)·((k v) × (k w)) = k^3 u·(v×w)
)
