from __future__ import annotations
import math
from itertools import permutations as _iperms, product as _product
from .atoms import Atom
from .context import Plan


def _full_group_size(plan: Plan) -> int:
    """Size of the full symmetry group: product of factorials of group sizes."""
    return math.prod(math.factorial(g.size) for g in plan.context.groups)


def orbit(atom_tuple: tuple, plan: Plan) -> list:
    """
    Return the list of distinct canonical atom tuples in the orbit of
    atom_tuple under the full symmetry group.

    The input atom_tuple need not itself be in canonical form; the orbit
    is the same regardless of which representative is supplied.
    """
    context    = plan.context
    seen       = set()
    result     = []

    group_perms = [list(_iperms(g.labels)) for g in context.groups]

    for perm_combo in _product(*group_perms):
        relabeling = {}
        for group, perm in zip(context.groups, perm_combo):
            for orig, new in zip(group.labels, perm):
                relabeling[orig] = new

        relabeled = tuple(
            Atom(
                operation=a.operation,
                labels=tuple(relabeling.get(lbl, lbl) for lbl in a.labels),
                sign=a.sign,
            )
            for a in atom_tuple
        )

        canonical = plan.canonicalise(relabeled)

        if canonical not in seen:
            seen.add(canonical)
            result.append(canonical)

    return result


def stabiliser_size(atom_tuple: tuple, plan: Plan) -> int:
    """
    Return the size of the stabiliser of atom_tuple (the number of group
    elements that map it to the same canonical form).
    By the orbit-stabiliser theorem: |G| / |orbit|.
    """
    return _full_group_size(plan) // len(orbit(atom_tuple, plan))


def orbit_and_stabiliser_size(atom_tuple: tuple, plan: Plan) -> tuple:
    """Return (orbit_list, stabiliser_size) together to avoid computing twice."""
    orb       = orbit(atom_tuple, plan)
    stab_size = _full_group_size(plan) // len(orb)
    return orb, stab_size
