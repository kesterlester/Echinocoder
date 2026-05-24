"""
symcoder.operations — library of standard physics operations.

Available modules (nothing is pre-loaded; import what you need):

    euclidean2   2-D Euclidean space:      mag, dot, eps  (rank 2 antisymmetric)
    euclidean3   3-D Euclidean space:      mag, dot, eps  (rank 3 antisymmetric)
    minkowski    3+1-D Minkowski (+−−−):   dot, mag_sq, eps  (rank 4 antisymmetric)

Each module exposes module-level Operation singletons.  Repeated imports of
the same name return the identical Python object, so ``is`` reliably tests
whether two references came from the same library entry::

    from symcoder.operations.euclidean3 import eps
    from symcoder.operations.euclidean3 import eps as eps_again
    assert eps is eps_again          # True  — same singleton

    from symcoder.operations.euclidean2 import eps as eps2d
    assert eps is not eps2d          # True  — different operations

Typical usage::

    from symcoder.operations import euclidean3 as SO3
    plan = Plan(context=ctx, operations=(SO3.mag, SO3.dot, SO3.eps))

    # or cherry-pick
    from symcoder.operations.euclidean3 import dot, eps
"""
