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


def eval_single_orbit_compressed(fo, plan: Plan, event: dict) -> list:
    """
    Return Phase 1 values for fo, exploiting sign compression (5c) for
    ANTISYMMETRIC operations.

    ANTISYMMETRIC: the orbit under repL contains pairs {+u, -u} for every
    label combination, so eval values come in {y, -y} pairs.  Only the
    sign=+1 atoms are evaluated, and their absolute values are returned.
    This is permutation-invariant: label permutations may swap which atom
    carries sign=+1 vs -1, but |eval| is unchanged, and the multiset
    {|y_i|} is invariant.  Returns fo.count() // 2 values.

    SYMMETRIC / UNSTRUCTURED: all atoms have sign=+1 and the full eval
    list is returned (identical to eval_single_orbit).  Returns fo.count()
    values.

    The distinction is made entirely from fo.operation.argument_symmetry —
    never from the eval values themselves.
    """
    from symatom.atoms import ArgumentSymmetry
    if fo.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC:
        return [abs(evaluate(atom, event)) for atom in fo.atoms() if atom.sign == +1]
    return [evaluate(atom, event) for atom in fo.atoms()]


def _is_self_pair(pf) -> bool:
    """True if every atom-pair (u_k, v_k) in this PairFlavour pairs an atom with itself.

    Condition: same operation, same flavour counts, and full (maximum) overlap.
    Under this condition z_k = eval(u_k) + i*eval(u_k) = (1+i)*a_k, which is
    entirely determined by Phase 1 orbit storage — no new information is added.
    """
    return (pf.op_u == pf.op_v and
            pf.flavour_u.counts == pf.flavour_v.counts and
            tuple(pf.overlap) == tuple(pf.flavour_u.counts))


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
