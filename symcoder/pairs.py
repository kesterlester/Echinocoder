from __future__ import annotations
from symatom.context import Plan
from symatom.rep import PairFlavour
from .eval import evaluate


def eval_pair_orbit(pf: PairFlavour, plan: Plan, event: dict) -> list:
    """
    Return the complex numbers z_k = eval(u'_k, E) + i*eval(v'_k, E) for
    every atom-pair (u'_k, v'_k) in the G-orbit of the canonical
    representative of this PairFlavour.

    For ANTISYMMETRIC operations both sign variants appear in the orbit —
    the orbit_elements enumeration in symatom handles this.
    """
    return [
        complex(evaluate(u, event), evaluate(v, event))
        for u, v in plan.orbit_enumerator.orbit_elements(pf, plan.context)
    ]
