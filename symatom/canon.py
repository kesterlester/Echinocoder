from __future__ import annotations
from itertools import permutations as _iperms, product as _product
from typing import runtime_checkable, Protocol
from .atoms import Atom


def _atom_key(atom: Atom) -> tuple:
    """A comparable sort key for a single Atom."""
    return (atom.operation.name, atom.operation.rank, atom.labels, atom.sign)


def _atom_tuple_key(atom_tuple: tuple) -> tuple:
    """A comparable sort key for a tuple of Atoms."""
    return tuple(_atom_key(a) for a in atom_tuple)


# ---------------------------------------------------------------------------
# Canonicaliser protocol
# ---------------------------------------------------------------------------

@runtime_checkable
class Canonicaliser(Protocol):
    """
    Any object satisfying this protocol may be used as the canonicaliser in a
    Plan.  The only required method is canonicalise().
    """
    def canonicalise(self, atom_tuple: tuple, context) -> tuple:
        """
        Return the canonical form of atom_tuple under the full symmetry group
        determined by context.  Must satisfy the contract C1–C5 in SPEC.md.
        All sign changes are absorbed into the atoms' sign fields; no external
        sign is returned.
        """
        ...


# ---------------------------------------------------------------------------
# Concrete implementation
# ---------------------------------------------------------------------------

class SimpleCanonicaliser:
    """
    A brute-force canonicaliser.

    Canonical form = the lexicographically smallest atom tuple obtained by:
      1. Trying every combination of permutations of the labels in each
         VectorGroup (i.e. every element of the full symmetry group).
      2. For each relabeling: construct relabeled Atoms (the Atom constructor
         automatically sorts each atom's own argument list and adjusts sign
         for ANTISYMMETRIC operations).
      3. Sort the atoms within the tuple into a canonical order.
      4. Return the smallest such tuple.

    Correctness: O(n1! * n2! * ...) where ni are the group sizes.  Fine for
    the small groups (4–6 vectors) this package targets.  A faster
    implementation can be swapped in later without changing the interface.

    Note on C4: the sign consistency property C4 as written in SPEC.md rev 2
    is under review — see the discussion in tests/test_canon.py.  C1, C2, C3,
    and C5 are all satisfied by this implementation.
    """

    def canonicalise(self, atom_tuple: tuple, context) -> tuple:
        if not atom_tuple:
            return atom_tuple

        best_key = None
        best     = None

        group_perms = [list(_iperms(g.labels)) for g in context.groups]

        for perm_combo in _product(*group_perms):
            # Build the label relabeling for this group-permutation combination.
            relabeling = {}
            for group, perm in zip(context.groups, perm_combo):
                for orig, new in zip(group.labels, perm):
                    relabeling[orig] = new

            # Apply relabeling; the Atom constructor handles internal arg sorting.
            relabeled = [
                Atom(
                    operation=a.operation,
                    labels=tuple(relabeling.get(lbl, lbl) for lbl in a.labels),
                    sign=a.sign,
                )
                for a in atom_tuple
            ]

            # Sort atoms within the tuple for a fully canonical ordering.
            candidate     = tuple(sorted(relabeled, key=_atom_key))
            candidate_key = _atom_tuple_key(candidate)

            if best_key is None or candidate_key < best_key:
                best_key = candidate_key
                best     = candidate

        return best
