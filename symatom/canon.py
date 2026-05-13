from __future__ import annotations
from itertools import permutations as _iperms, product as _product
from typing import runtime_checkable, Protocol
from .atoms import Atom, ArgumentSymmetry


def _atom_key(atom: Atom) -> tuple:
    """A comparable sort key for a single Atom."""
    return (atom.operation.name, atom.operation.rank, atom.labels, -atom.sign)


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
         VectorType (i.e. every element of the full symmetry group).
      2. For each relabeling: construct relabeled Atoms (the Atom constructor
         automatically sorts each atom's own argument list and adjusts sign
         for ANTISYMMETRIC operations).
      3. Sort the atoms within the tuple into a canonical order.
      4. Return the smallest such tuple.

    Correctness: O(n1! * n2! * ...) where ni are the group sizes.  Fine for
    the small groups (4–6 vectors) this package targets.  A faster
    implementation can be swapped in later without changing the interface.

    """

    def canonicalise(self, atom_tuple: tuple, context) -> tuple:
        if not atom_tuple:
            return atom_tuple

        best_key = None
        best     = None

        type_perms = [list(_iperms(g.labels)) for g in context.types]

        for perm_combo in _product(*type_perms):
            # Build the label relabeling for this group-permutation combination.
            relabeling = {}
            for vt, perm in zip(context.types, perm_combo):
                for orig, new in zip(vt.labels, perm):
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


class DirectCanonicaliser:
    """
    Fast canonicaliser that exploits the fact that only the labels which
    actually appear in the atom tuple need to be permuted.

    Background and motivation
    --------------------------
    SimpleCanonicaliser tries every element of the full symmetry group
    G = S_{n1} x S_{n2} x ... (one symmetric group per VectorType).  For a
    group of size n, that is n! permutations.  At n=10 this is ~3.6 million;
    at n=30 it is ~2.7e32 — completely intractable.

    The atom tuple is the key observation: an atom tuple with k atoms of
    maximum rank r uses at most k*r DISTINCT labels from each group, and
    typically far fewer.  Call this number m (m << n for large n).  Labels
    that do not appear in any atom are "inactive" — permuting them has no
    effect on the atom tuple, so they need not be enumerated.  This reduces
    the search from n! to m! per group, which is manageable even for n=30
    provided m stays small (a pair of rank-2 and rank-3 atoms uses m <= 5).

    Key theorem: active labels should always map to the m smallest targets
    -----------------------------------------------------------------------
    Let the active labels of group g be {l1, ..., lm} (those appearing in
    the tuple), and let the full group have labels sorted as
    t1 < t2 < ... < tn.  The canonical targets are always {t1, ..., tm}.

    Proof by exchange argument: suppose the lex-min permutation pi* maps
    some active label to a target t_j with j > m (i.e. a label LARGER than
    the m-th smallest).  Then there exists some t_i (i <= m) that is not
    used as a target for any active label.  t_i must be assigned to some
    INACTIVE label.  Construct pi' identical to pi* except pi' swaps the
    assignments of t_i and t_j.  Because t_i < t_j, every atom that
    previously contained t_j now contains the smaller t_i.  After re-sorting
    the atom tuple, the key is lex-smaller or equal.  But pi* was lex-min,
    so pi' must give the same key — meaning there exists an equally lex-min
    permutation that uses t_i instead of t_j.  Applying the argument
    inductively over all j > m, there is always a lex-min permutation whose
    active labels map into {t1, ..., tm}.

    Algorithm (per canonicalise call)
    -----------------------------------
    For each VectorType g:
      1. Collect active_g = labels from g that appear in at least one atom.
         (Inactive labels are ignored entirely.)
      2. Set canonical_targets_g = sorted(g.labels)[:len(active_g)].
         These are the m lexicographically smallest labels in g; by the
         theorem above, the lex-min canonical form always uses exactly these
         labels for the active slots.
      3. Enumerate all len(active_g)! bijections from active_g to
         canonical_targets_g.  Each bijection is a partial relabeling dict.

    Take the Cartesian product of the bijection lists across all groups.
    For each combined relabeling:
      - Apply it to every atom (the Atom constructor re-sorts internal
        argument order and adjusts sign for ANTISYMMETRIC operations).
      - Sort the resulting atoms into a tuple using _atom_key.
      - Track the lex-min tuple seen so far.

    Return the lex-min tuple.

    Complexity
    -----------
    O(m1! x m2! x ... x |atom_tuple| x max_rank)
    where mi = number of distinct labels from group i appearing in the tuple.

    Typical speedup vs SimpleCanonicaliser:
      n=10, m=3:   3! / 10! =  6 / 3628800   — ~600,000x faster
      n=30, m=5:   5! / 30! = 120 / 2.7e32   — ~2e30x faster

    Worst case (atom tuple uses every label in every group): m=n, same
    as SimpleCanonicaliser.  In practice, for the orbit-enumeration workload
    this canonicaliser is designed for, m is bounded by the sum of atom ranks
    (typically 2-6 per group), making the speedup enormous.

    Correctness
    -----------
    Validated against SimpleCanonicaliser by the parametrised cross-comparison
    tests in symatom/tests/test_canon.py, which assert identical output on
    every test case.  SimpleCanonicaliser is retained permanently as the
    reference.
    """

    def canonicalise(self, atom_tuple: tuple, context) -> tuple:
        if not atom_tuple:
            return atom_tuple

        # Collect all labels appearing in the tuple (across all groups)
        all_active = {lbl for a in atom_tuple for lbl in a.labels}

        # For each group build the list of bijections: active -> canonical targets
        group_bijection_lists = []
        for g in context.types:
            active = [l for l in g.labels if l in all_active]
            if not active:
                group_bijection_lists.append([{}])
                continue
            # Canonical targets: the len(active) lex-smallest labels in g
            targets = sorted(g.labels)[:len(active)]
            # One bijection per permutation of targets (active order is fixed;
            # we permute which target each active label receives)
            group_bijection_lists.append(
                [dict(zip(active, perm)) for perm in _iperms(targets)]
            )

        best_key = None
        best     = None

        for combo in _product(*group_bijection_lists):
            # Merge per-group dicts into one relabeling
            relabeling = {}
            for mapping in combo:
                relabeling.update(mapping)

            # Apply relabeling; Atom constructor handles internal arg sort + sign
            relabeled = [
                Atom(
                    operation=a.operation,
                    labels=tuple(relabeling.get(lbl, lbl) for lbl in a.labels),
                    sign=a.sign,
                )
                for a in atom_tuple
            ]

            candidate     = tuple(sorted(relabeled, key=_atom_key))
            candidate_key = _atom_tuple_key(candidate)

            if best_key is None or candidate_key < best_key:
                best_key = candidate_key
                best     = candidate

        return best
