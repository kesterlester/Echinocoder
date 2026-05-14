"""
symcoder.encoders.poly_encoder
================================
PolyEncoder: polynomial (compressed) embedding of a single atom orbit.

Accepts FLAVOURED_OPERATOR and REPRESENTATIVE_ATOM specs.  Returns can_encode=False
for EXPLICIT_ORBIT because a FlavouredOperator is needed to recover the sign-
correlation structure; if that need arises, reconstruct a FlavouredOperator from
the explicit atoms and use from_flavoured_operator() instead.

Priority is 1.0 — higher than SortEncoder — so it is preferred when both are
registered and both are capable.

Relationship to encode.py
--------------------------
The existing symcoder/encode.py contains _embed_compressed(), which encodes
*pair* orbits (two atoms evaluated together to form a complex z_k value).  The
logic for single-atom orbits is simpler: the Phase-1 sorted evaluation in encode()
is already the permutation-invariant representation, and the polynomial embedding
adds nothing for purely real values from symmetric operations.

For ANTISYMMETRIC operations the orbit is {a_k, -a_k}, and the natural
polynomial approach squares the values (like TYPE_NEG in _embed_compressed):

    w_k = a_k²    (n values, each in R≥0)
    embed w as a multiset via _zip_embed → n complex polynomial coefficients
    → 2n reals

This halves the output dimension compared to SortEncoder for ANTISYMMETRIC
single-atom orbits (since the sorted list has the redundant {-a_k} entries).

For SYMMETRIC operations the orbit values are all positive (or mixed sign but
without the {a,-a} structure), and polynomial embedding gives no compression
benefit over sorting.  In that case PolyEncoder may choose to delegate to
SortEncoder or simply return can_encode=False for SYMMETRIC operations.

Internal delegation
-------------------
This encoder holds a SortEncoder instance internally.  After obtaining the
sorted real evaluations it applies polynomial compression on top:

    self._sort_enc.encode(spec, event, plan).values → sorted float64 array
    # then apply polynomial compression

Implementation guide (stubs below raise NotImplementedError)
------------------------------------------------------------

assess()
~~~~~~~~
from symatom.atoms import ArgumentSymmetry
from symatom.rep import FlavouredOperator, Flavour

# 1. Obtain the FlavouredOperator:
if spec.form == OrbitSpecForm.FLAVOURED_OPERATOR:
    fo = spec.payload
elif spec.form == OrbitSpecForm.REPRESENTATIVE_ATOM:
    atom = spec.payload
    flavour = Flavour(tuple(
        sum(1 for lbl in atom.labels if lbl in set(g.labels))
        for g in plan.context.types
    ))
    fo = FlavouredOperator(atom.operation, flavour, plan.context)
else:
    return EncodingCapability(can_encode=False, output_dim=None,
                              method_name=None, priority=0.0)

# 2. Check the operation is evaluable (has eval_fn):
from symcoder.eval import EvaluableOperation
if not isinstance(fo.operation, EvaluableOperation):
    return EncodingCapability(can_encode=False, output_dim=None,
                              method_name=None, priority=0.0,
                              metadata={"reason": "operation has no eval_fn"})

# 3. Only ANTISYMMETRIC operations benefit from polynomial compression here.
#    For SYMMETRIC, delegate to SortEncoder (or return can_encode=False).
if fo.operation.argument_symmetry != ArgumentSymmetry.ANTISYMMETRIC:
    return EncodingCapability(can_encode=False, output_dim=None,
                              method_name=None, priority=0.0,
                              metadata={"reason": "no compression benefit for SYMMETRIC"})

# 4. Compute output_dim.
#    The orbit has fo.count_of_atoms_one_per_sign() atoms; after squaring and polynomial embedding
#    of the positive half the output is 2 * (fo.count_of_atoms_one_per_sign() // 2) reals.
#    (For ANTISYMMETRIC orbits fo.count_of_atoms_one_per_sign() is always even.)
n = fo.count_of_atoms_one_per_sign()
output_dim = n  # one complex coeff per squared root → 2n reals via _complex_to_reals
               # Revise once the exact embedding is decided.

return EncodingCapability(
    can_encode=True,
    output_dim=output_dim,
    method_name="poly_antisym",
    priority=_PRIORITY,
    metadata={"orbit_size": fo.count_of_atoms_one_per_sign(), "fo": fo},
)

encode()
~~~~~~~~
# 1. Obtain sorted evaluations from the sort encoder:
sorted_vals = self._sort_enc.encode(spec, event, plan).values
# sorted_vals is float64 array of length fo.count_of_atoms_one_per_sign(), symmetric about 0 for ANTISYM.

# 2. Take the positive half (second half of sorted array for ANTISYM orbits):
n = len(sorted_vals)
pos_vals = sorted_vals[n // 2:]   # the n//2 non-negative values

# 3. Square to form the sign-invariant multiset:
w = pos_vals ** 2

# 4. Polynomial-embed w as a multiset of reals using _zip_embed:
from symcoder.encode import _zip_embed, _complex_to_reals
# _zip_embed expects complex input; cast:
coeffs, _, _ = _zip_embed(w.astype(complex))
values = _complex_to_reals(coeffs)  # 2*(n//2) reals

return EncodingResult(
    values=values,
    metadata={"method": "poly_antisym", "orbit_size": n},
)

Note on future pair-orbit encoders
------------------------------------
The full _embed_compressed() logic in encode.py targets *pair* orbits (PairFlavour).
A future PolyPairEncoder should accept an OrbitSpec wrapping a PairFlavour and
delegate directly to _embed_compressed().  That encoder would live in a separate
file (e.g. poly_pair_encoder.py) and register alongside this one; the registry
manager would then pick between them based on the spec form.
"""
from __future__ import annotations

import numpy as np

from symatom.context import Plan

from ._base import (
    AtomOrbitEncoder,
    OrbitSpec,
    OrbitSpecForm,
    EncodingCapability,
    EncodingResult,
)
from .sort_encoder import SortEncoder

_PRIORITY = 1.0


class PolyEncoder(AtomOrbitEncoder):
    """
    Polynomial-compressed embedding for ANTISYMMETRIC single-atom orbits.

    Holds a SortEncoder internally and uses it to obtain Phase-1 sorted
    evaluations before applying compression.

    See module docstring for the full implementation guide.
    """

    def __init__(self) -> None:
        self._sort_enc = SortEncoder()

    def assess(self, spec: OrbitSpec, plan: Plan) -> EncodingCapability:
        # TODO: implement — see module docstring for step-by-step guide
        raise NotImplementedError(
            "PolyEncoder.assess() is a stub — see module docstring for the "
            "implementation guide (FlavouredOperator extraction, ANTISYM check, "
            "output_dim calculation)"
        )

    def encode(self, spec: OrbitSpec, event: dict, plan: Plan) -> EncodingResult:
        # TODO: implement — see module docstring for step-by-step guide
        raise NotImplementedError(
            "PolyEncoder.encode() is a stub — see module docstring for the "
            "implementation guide (delegate to _sort_enc, square, _zip_embed)"
        )
