from __future__ import annotations
from symatom.context import Plan
from symatom.rep import PairFlavour
from .eval import evaluate


def eval_pair_orbit(pf: PairFlavour, plan: Plan, event: dict) -> list:
    """
    Return the complex numbers z_k = eval(u'_k, E) + i*eval(v'_k, E) for
    every atom-pair (u'_k, v'_k) in the G-orbit of the canonical
    representative of this PairFlavour.

    Returns pf.orbit_size(type_sizes) values.  Used as the brute-force
    reference and in group-theory tests.
    """
    return [
        complex(evaluate(u, event), evaluate(v, event))
        for u, v in plan.orbit_enumerator.orbit_elements(pf, plan.context)
    ]


## def eval_single_orbit(fo, plan: Plan, event: dict) -> list:
##     """
##     Return eval(u, event) for every atom u in the orbit of a single
##     FlavouredOperator fo.
## 
##     Uses fo.atoms_one_per_sign() directly — no pair machinery needed.
## 
##     Used in Phase 1 of encode() to sort-encode the single-atom orbit before
##     processing pair associations.
##     """
##     return [evaluate(atom, event) for atom in fo.atoms_one_per_sign()]


## def eval_single_orbit_compressed(fo, plan: Plan, event: dict) -> list:
##     print("FIXME : eval_single_orbit_compressed is wrong at present:")
##     # TODO FIXME As now using repS rather than repL we might try to output an
##     # eval of sort([eps2(a,b)]) which is not perm invariant.  (Under reps we would
##     # have returned sort([+eps2(a,b), -eps2(a,b)]) which would be OK.
##     # Furthermore it is wrong that sign compression happens only when the
##     # operations are antisymmetric.  Consider, for example, that eps2(a,p)
##     # does not change sign under any permutation of vectors. So lots wrong below.
##     # TODO: What is really needed is a way of ouputting the phase 1 numbers in a perm invariant way when that is needed, but as they are when they are already fine.
##     """
##     Return Phase 1 values for fo, exploiting sign compression (5c) for
##     ANTISYMMETRIC operations.
## 
##     ANTISYMMETRIC: the orbit under repL contains pairs {+u, -u} for every
##     label combination, so eval values come in {y, -y} pairs.  Only the
##     sign=+1 atoms are evaluated, and their absolute values are returned.
##     This is permutation-invariant: label permutations may swap which atom
##     carries sign=+1 vs -1, but |eval| is unchanged, and the multiset
##     {|y_i|} is invariant.  Returns fo.count_of_atoms_one_per_sign() // 2 values.
## 
##     SYMMETRIC / UNSTRUCTURED: all atoms have sign=+1 and the full eval
##     list is returned (identical to eval_single_orbit).  Returns fo.count_of_atoms_one_per_sign()
##     values.
## 
##     The distinction is made entirely from fo.operation.argument_symmetry —
##     never from the eval values themselves.
##     """
##     from symatom.atoms import ArgumentSymmetry
##     if fo.operation.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC:
##         return [abs(evaluate(atom, event)) for atom in fo.atoms_one_per_sign() if atom.sign == +1]
##     return [evaluate(atom, event) for atom in fo.atoms_one_per_sign()]


def _is_self_pair(pf) -> bool:
    """True if every atom-pair (u_k, v_k) in this PairFlavour pairs an atom with itself.

    Condition: same operation, same flavour counts, and full (maximum) overlap.
    Under this condition z_k = eval(u_k) + i*eval(u_k) = (1+i)*a_k, which is
    entirely determined by Phase 1 orbit storage — no new information is added.
    """
    return (pf.op_u == pf.op_v and
            pf.flavour_u.counts == pf.flavour_v.counts and
            tuple(pf.overlap) == tuple(pf.flavour_u.counts))

