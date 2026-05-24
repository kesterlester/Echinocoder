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
)

dot = Operation(
    "dot", rank=2, odd_parity=False,
    argument_symmetry=ArgumentSymmetry.SYMMETRIC,
    eval_fn=lambda v: float(np.dot(v[0], v[1])),
    tex=r"\mathbf{#1} \cdot \mathbf{#2}",
)

eps = Operation(
    "eps", rank=3, odd_parity=True,
    argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
    eval_fn=lambda v: float(np.dot(v[0], np.cross(v[1], v[2]))),
    tex=r"\varepsilon(\mathbf{#1},\mathbf{#2},\mathbf{#3})",
)
