from __future__ import annotations
from abc import ABC, abstractmethod
from itertools import combinations, product as iproduct
from symatom.atoms import ArgumentSymmetry, Atom


class OrbitEnumerator(ABC):
    """
    Strategy interface for enumerating the G-orbit of an atom-pair given its
    PairFlavour.  Two implementations are provided:

    - BruteForceOrbitEnumerator  (correct by inspection; always available)
    - DirectOrbitEnumerator      (fast combinatorial construction)

    Pass the desired enumerator to Plan; symcoder calls
    plan.orbit_enumerator.orbit_elements(pf, context) without caring which
    implementation is active.
    """

    @abstractmethod
    def orbit_elements(self, pf, context) -> list:
        """
        Return all (Atom, Atom) pairs in the G-orbit of pf's canonical
        representative, given a context.
        """


class BruteForceOrbitEnumerator(OrbitEnumerator):
    """
    Cross-product + filter implementation.  Delegates to pf.orbit_elements(),
    which enumerates FlavouredOperator.atoms() for each side and keeps only
    pairs whose pair_flavour_of matches pf.

    Correct by inspection: every atom in FlavouredOperator.atoms() is a
    distinct concrete atom with the right flavour; pair_flavour_of is the
    canonical orbit-type test.  Keep this permanently as the reference
    against which DirectOrbitEnumerator is validated.
    """

    def orbit_elements(self, pf, context) -> list:
        return pf.orbit_elements(context)


class DirectOrbitEnumerator(OrbitEnumerator):
    """
    Combinatorial construction from the PairFlavour overlap structure.

    For large n (target: n ~ 30 per species) the cross-product approach in
    BruteForceOrbitEnumerator is O(atoms_u x atoms_v) >> O(orbit_size).
    This implementation iterates directly over the combinatorial structure
    encoded in pf.overlap, producing exactly orbit_size pairs with no
    filtering step anywhere.

    Core idea: label partitioning per group
    ---------------------------------------
    For each group g (with n labels), and a given PairFlavour specifying that
    atom u uses ku labels from this group and atom v uses kv labels with
    overlap s, we directly enumerate all valid ways to assign concrete labels
    from g.labels to u and v.  No pool of atoms is constructed, and no
    filtering step is applied.

    For each group the loop does this:

      1. Choose s shared labels from g.labels: combinations(g.labels, s).
         These go into both u and v.
      2. From the remaining n-s labels, choose ku-s u-only labels:
         combinations(rest1, ku-s).  These go into u only.
      3. From the remaining n-ku labels, choose kv-s v-only labels:
         combinations(rest2, kv-s).  These go into v only.

    The total count for one group is C(n,s) x C(n-s, ku-s) x C(n-ku, kv-s),
    which is exactly what pf.count() computes.  No atoms are ever generated
    and then discarded.

    Cross-product across groups
    ---------------------------
    After computing the list of (shared, u_only, v_only) triples for each
    group independently, itertools.product(*per_group) takes the Cartesian
    product across groups.  For each combined assignment, the labels for u
    and v are assembled by concatenation:

      labels_u = shared_g1 + u_only_g1 + shared_g2 + u_only_g2 + ...
      labels_v = shared_g1 + v_only_g1 + shared_g2 + v_only_g2 + ...

    Atom construction and sign handling
    ------------------------------------
    Atom(op, labels, sign=+1) is called with the concatenated (potentially
    unsorted) label tuple.  The Atom constructor internally sorts the labels
    and computes the permutation sign via perm_sign, so e.g.
    Atom(op, ('c', 'a'), +1) correctly produces the atom with sorted labels
    ('a', 'c') and sign=-1.  The label order passed in does not matter; the
    constructor handles it.

    For ANTISYMMETRIC operations, both sign variants are needed in the orbit.
    The code generates up to 4 sign combinations per label assignment:

      (u_pos, v_pos)   -- always
      (-u,    v_pos)   -- if op_u is ANTISYMMETRIC
      (u_pos, -v)      -- if op_v is ANTISYMMETRIC
      (-u,    -v)      -- if both are ANTISYMMETRIC

    There is no generate-and-prune anywhere.  The only combinatorial
    primitives used are combinations (for label selection within a group)
    and product (for cross-group assembly).  Every produced atom-pair is a
    valid orbit element by construction.

    The resulting set of atom-pairs is identical to BruteForceOrbitEnumerator;
    the parametrised cross-comparison tests in test_orbit_enum.py verify this.
    """

    def orbit_elements(self, pf, context) -> list:
        # For each group build all valid (shared, u_only, v_only) label partitions
        per_group = []
        for g, ku, kv, s in zip(
            context.groups,
            pf.flavour_u.counts,
            pf.flavour_v.counts,
            pf.overlap,
        ):
            group_choices = []
            for shared in combinations(g.labels, s):
                shared_set = set(shared)
                rest1 = [l for l in g.labels if l not in shared_set]
                for u_only in combinations(rest1, ku - s):
                    u_only_set = set(u_only)
                    rest2 = [l for l in rest1 if l not in u_only_set]
                    for v_only in combinations(rest2, kv - s):
                        group_choices.append((shared, u_only, v_only))
            per_group.append(group_choices)

        antisym_u = pf.op_u.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC
        antisym_v = pf.op_v.argument_symmetry == ArgumentSymmetry.ANTISYMMETRIC

        result = []
        for combo in iproduct(*per_group):
            labels_u = tuple(
                lbl for shared, u_only, _ in combo for lbl in shared + u_only
            )
            labels_v = tuple(
                lbl for shared, _, v_only in combo for lbl in shared + v_only
            )
            u_pos = Atom(pf.op_u, labels_u, sign=+1)
            v_pos = Atom(pf.op_v, labels_v, sign=+1)

            result.append((u_pos, v_pos))
            if antisym_u:
                result.append((Atom(pf.op_u, labels_u, sign=-1), v_pos))
            if antisym_v:
                result.append((u_pos, Atom(pf.op_v, labels_v, sign=-1)))
            if antisym_u and antisym_v:
                result.append((
                    Atom(pf.op_u, labels_u, sign=-1),
                    Atom(pf.op_v, labels_v, sign=-1),
                ))

        return result
