"""
symcoder.encoders.poly_encoder
================================
PolyEncoderFactory / PolyEncoder: polynomial (compressed) embedding of a
single atom orbit.

Accepts FLAVOURED_OPERATOR and REPRESENTATIVE_ATOM specs.  Returns [] for
EXPLICIT_ORBIT because a FlavouredOperator is needed to recover the sign-
correlation structure; if that need arises, reconstruct a FlavouredOperator
from the explicit atoms and use from_flavoured_operator() instead.

Priority is 1.0 — higher than SortEncoderFactory — so it is preferred when
both are registered and both are capable.

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
benefit over sorting.  PolyEncoderFactory returns [] for SYMMETRIC operations.

Implementation guide (stubs below raise NotImplementedError)
------------------------------------------------------------

PolyEncoderFactory.assess()
~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    return []   # EXPLICIT_ORBIT not supported

# 2. Check the operation is evaluable (has eval_fn):
from symcoder.eval import EvaluableOperation
if not isinstance(fo.operation, EvaluableOperation):
    return []

# 3. Only ANTISYMMETRIC operations benefit from polynomial compression here.
if fo.operation.argument_symmetry != ArgumentSymmetry.ANTISYMMETRIC:
    return []

# 4. Compute output_dim.
#    The orbit has fo.count_of_atoms_one_per_sign() atoms; after squaring and
#    polynomial embedding of the positive half the output is
#    2 * (fo.count_of_atoms_one_per_sign() // 2) reals.
#    (For ANTISYMMETRIC orbits fo.count_of_atoms_one_per_sign() is always even.)
n = fo.count_of_atoms_one_per_sign()

return [PolyEncoder(fo=fo, n=n)]


PolyEncoder.encode()
~~~~~~~~~~~~~~~~~~~~
# 1. Obtain sorted evaluations from a SortEncoderFactory for the same orbit:
#    (The factory holds a SortEncoderFactory and ran it during __init__.)
sorted_vals = self._sort_encoder.encode(event).values
# sorted_vals: float64 array of length n, symmetric about 0 for ANTISYM.

# 2. Take the positive half (second half of sorted array for ANTISYM orbits):
pos_vals = sorted_vals[self._n // 2:]   # the n//2 non-negative values

# 3. Square to form the sign-invariant multiset:
w = pos_vals ** 2

# 4. Polynomial-embed w as a multiset of reals using _zip_embed:
from symcoder.encode import _zip_embed, _complex_to_reals
coeffs, _, _ = _zip_embed(w.astype(complex))
values = _complex_to_reals(coeffs)  # 2*(n//2) reals

return EncodingResult(
    values=values,
    metadata={"method": "poly_antisym", "orbit_size": self._n},
)

Note on future pair-orbit encoders
------------------------------------
The full _embed_compressed() logic in encode.py targets *pair* orbits (PairFlavour).
A future PolyPairEncoderFactory should accept an OrbitSpec wrapping a PairFlavour
and delegate directly to _embed_compressed().
"""
from __future__ import annotations

import numpy as np

from symatom.context import Plan

from ._base import (
    AtomOrbitEncoder,
    AtomOrbitEncoderFactory,
    OrbitSpec,
    OrbitSpecForm,
    EncodingResult,
)

_PRIORITY    = 1.0
_METHOD_NAME = "poly_antisym"


class PolyEncoder(AtomOrbitEncoder):
    """
    Polynomial-compressed embedding for ANTISYMMETRIC single-atom orbits.

    Constructed by PolyEncoderFactory.assess() with fo and n pre-computed.
    Internally holds a SortEncoder for the same orbit to obtain Phase-1
    sorted evaluations before applying compression.

    See module docstring for the full implementation guide.
    """

    def __init__(self, fo, n: int, sort_encoder: AtomOrbitEncoder) -> None:
        self._fo           = fo
        self._n            = n
        self._sort_encoder = sort_encoder  # pre-built SortEncoder for same orbit

    @property
    def output_dim(self) -> int:
        return self._n  # 2*(n//2) reals — revise once exact embedding is decided

    @property
    def priority(self) -> float:
        return _PRIORITY

    @property
    def method_name(self) -> str:
        return _METHOD_NAME

    def encode(self, event: dict) -> EncodingResult:
        # TODO: implement — see module docstring for step-by-step guide
        raise NotImplementedError(
            "PolyEncoder.encode() is a stub — see module docstring for the "
            "implementation guide (delegate to _sort_encoder, square, _zip_embed)"
        )


class PolyEncoderFactory(AtomOrbitEncoderFactory):
    """
    Factory that produces PolyEncoder instances for ANTISYMMETRIC single-atom orbits.

    Internally uses a SortEncoderFactory to pre-build the sort encoder that
    PolyEncoder will delegate to at encode time.

    See module docstring for the full implementation guide.
    """

    def __init__(self) -> None:
        from .sort_encoder import SortEncoderFactory
        self._sort_factory = SortEncoderFactory()

    def assess(self, spec: OrbitSpec, plan: Plan) -> list[AtomOrbitEncoder]:
        # TODO: implement — see module docstring for step-by-step guide
        raise NotImplementedError(
            "PolyEncoderFactory.assess() is a stub — see module docstring for the "
            "implementation guide (FlavouredOperator extraction, ANTISYM check, "
            "output_dim calculation)"
        )
