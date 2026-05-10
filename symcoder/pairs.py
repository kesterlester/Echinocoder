from __future__ import annotations
from symatom.context import Plan
from symatom.rep import PairFlavour
from .eval import evaluate


def eval_pair_orbit(pf: PairFlavour, plan: Plan, event: dict) -> list:
    """
    Return the complex numbers z_k = eval(u'_k, E) + i*eval(v'_k, E) for
    every atom-pair (u'_k, v'_k) in the G-orbit of the canonical
    representative of this PairFlavour.

    Returns pf.orbit_size(group_sizes) values.  Used as the brute-force
    reference and in group-theory tests; encode() uses eval_pair_orbit_positive
    instead.
    """
    return [
        complex(evaluate(u, event), evaluate(v, event))
        for u, v in plan.orbit_enumerator.orbit_elements(pf, plan.context)
    ]


def eval_single_orbit(fo, plan: Plan, event: dict) -> list:
    """
    Return eval(u, event) for every atom u in the orbit of a single
    FlavouredOperator fo.

    Uses fo.atoms() directly — no pair machinery needed.  The returned values
    are real scalars; for ANTISYMMETRIC operations both sign=+1 and sign=-1
    atoms are included (fo was created by repL with signed=True), so the list
    contains {a_k, -a_k} pairs.

    Used in Phase 1 of encode() to sort-encode the single-atom orbit before
    processing pair associations.
    """
    return [evaluate(atom, event) for atom in fo.atoms()]


def eval_pair_orbit_positive(pf: PairFlavour, plan: Plan, event: dict) -> list:
    """
    Return z_k = eval(u'_k, E) + i*eval(v'_k, E) for the (u.sign=+1, v.sign=+1)
    subset of the G-orbit.

    This always contains exactly pf.count(group_sizes) values regardless of
    which operations are ANTISYMMETRIC.  The proof: the Atom constructor sorts
    labels and absorbs the permutation parity into the stored sign, so for each
    base label combination the two u-atoms have opposite stored signs (+s and -s)
    — exactly one of the four (u,v)-sign combinations in an ANTISYM×ANTISYM
    orbit lands in the (+1,+1) cell.  For SYMMETRIC operations sign is always
    +1, so filtering changes nothing.

    Why this subset is sufficient for encoding (see _embed_compressed in
    encode.py): forming the symmetry-class-appropriate invariant multiset from
    these n values (e.g. {z_k, conj(z_k)} for SYM×ANTISYM) recovers the full
    permutation-invariant polynomial with its forced-zero / forced-real
    coefficient structure, yielding the same information as the full orbit
    embedding at a fraction of the cost.
    """
    return [
        complex(evaluate(u, event), evaluate(v, event))
        for u, v in plan.orbit_enumerator.orbit_elements(pf, plan.context)
        if u.sign == +1 and v.sign == +1
    ]
