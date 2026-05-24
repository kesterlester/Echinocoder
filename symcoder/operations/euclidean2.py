"""
symcoder.operations.euclidean2
==============================
Standard operations for 2-D Euclidean space.

All vectors are assumed to be 2-component numpy arrays.

Operations
----------
mag   rank 1, scalar        |v|
dot   rank 2, scalar        v·w
eps   rank 2, pseudoscalar  v[0]*w[1] - v[1]*w[0]   (2-D "cross product" / det)

Note
----
eps here is distinct from euclidean3.eps: same symbolic shape (rank 2,
antisymmetric, odd parity) but a different Python object and a different
eval_fn.  Use ``is`` to test provenance::

    from symcoder.operations.euclidean2 import eps as eps2
    from symcoder.operations.euclidean3 import eps as eps3   # if you need it
    assert eps2 is not eps3
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
    "eps", rank=2, odd_parity=True,
    argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
    eval_fn=lambda v: float(v[0][0]*v[1][1] - v[0][1]*v[1][0]),
    tex=r"\varepsilon(\mathbf{#1},\mathbf{#2})",
)
