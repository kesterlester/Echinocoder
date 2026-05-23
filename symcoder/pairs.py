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


def _is_self_pair(pf) -> bool:
    """True if every atom-pair (u_k, v_k) in this PairFlavour pairs an atom with itself.

    Condition: same operation, same flavour counts, and full (maximum) overlap.
    Under this condition z_k = eval(u_k) + i*eval(u_k) = (1+i)*a_k, which is
    entirely determined by Phase 1 orbit storage — no new information is added.
    """
    return (pf.op_u == pf.op_v and
            pf.flavour_u.counts == pf.flavour_v.counts and
            tuple(pf.overlap) == tuple(pf.flavour_u.counts))

