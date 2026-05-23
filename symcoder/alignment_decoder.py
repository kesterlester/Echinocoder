"""
symcoder.alignment_decoder
===========================
Alignment decoder: recovers the G-orbit of repS evaluations from Phase 1
and Phase 2 (OverlapBlock) decoded outputs.

Algorithm overview
------------------
The alignment table has |repS| rows (one canonical atom per FlavouredOperator)
and n_cols = |G| columns (one per group element g ∈ G).  Entry (r, c) =
eval(g_c · repS[r], E).

Given:
  - Phase 1 decoded values: the row multisets for each repS atom (with
    repetition, n_cols values per row).
  - Phase 2 decoded pair values: for cross-FO canonical pairs (u_can, v_can)
    where BOTH u_can and v_can are repS canonical atoms, the n_cols
    (u_val, v_val) joint values.

The decoder builds the columns by partitioning column indices and refining
the partition using pair constraints:

  1. Seed: assign the first atom's (sorted) row values to columns 0..n_cols-1.
  2. Refine: for each subsequent atom, for each group of columns sharing
     identical placed-atom values, use pair constraints to determine the
     distribution of new-atom values, then split the group by new-atom value.
  3. Collect: assemble column vectors as the output multiset.

Pair-relation constraint
------------------------
A decoded pair result for canonical pair (u_can, v_can) with u_can, v_can both
in repS provides, for each column g, the joint values
(eval(g·u_can, E), eval(g·v_can, E)) = (alignment[r_u, g], alignment[r_v, g]).

The decoder builds full_pairs = decoded.pairs × stabiliser_size_pair(u_can, v_can)
to expand the decoded orbit back to n_cols pairs.  Columns within the same
partition group all share the same value for u_can, so the pair constraint
directly gives the Counter of v_can values across those columns.

Multiple constraints on the same new atom are intersected via Counter &
Counter (minimum-count intersection), which resolves ambiguities without
special-casing.

Ambiguity and free fill
-----------------------
If no pair constraint links a new atom to any placed atom, values are assigned
from the remaining row-multiset pool in ascending sort order.  For a multiset
output this is always valid: any assignment within a group produces the same
set of column vectors.

Pairs with ONLY u_can in repS (v_can outside repS) are not usable as joint-
column constraints and are silently discarded.  The same applies to NULL_COMP
selections that own multiple orbits (their full_pairs list length exceeds n_cols
and is caught by the length sanity check).
"""
from __future__ import annotations

from collections import Counter, defaultdict

from symatom import repS as _repS_fn
from symcoder.decoded_types import AnnotatedMultisetOfRepSEvalVectors


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def decode_alignment(
    plan,
    phase1_results: dict,
    all_pair_decoded: list,
    atol: float = 0.0,
) -> AnnotatedMultisetOfRepSEvalVectors:
    """
    Recover the G-orbit of repS evaluations from Phase 1 and Phase 2 decoded outputs.

    Parameters
    ----------
    plan : Plan
        Encoding plan (holds context and operations).
    phase1_results : dict[(str, tuple), AnnotatedMultisetOfReals]
        Phase 1 decoded multisets keyed by (op_name, flavour_counts_tuple).
        Each entry's ``values`` list has orbit_size elements (one per distinct
        G-orbit element for that FlavouredOperator).
    all_pair_decoded : list[AnnotatedMultisetOfRealPairs]
        Flat list of all Phase 2 decoded pair multisets from all OverlapBlock
        decoders.  Self-pairs (NULL_SELF), self-pairing-block results, and
        multi-orbit NULL_COMP entries are filtered automatically.
    atol : float
        Absolute tolerance for matching reference values in pair lookups.
        Use 0.0 (default) when inputs are polished Phase 1 values (bitwise
        exact for integers; use 0.0 or a small epsilon for floats).

    Returns
    -------
    AnnotatedMultisetOfRepSEvalVectors
        Multiset of n_cols = |G| column vectors (fewer if the event has a
        non-trivial stabiliser on repS — but the multiset is always n_cols
        tuples in total, with identical vectors repeated for stabilised events).
        Each vector has length |repS|; position r holds eval(g_c · repS[r], E).
    """
    the_group = plan.context.the_group
    n_cols    = the_group.order()

    fo_list = _repS_fn(plan.context, plan.operations)
    if not fo_list or n_cols == 0:
        return AnnotatedMultisetOfRepSEvalVectors(vectors=[], atoms=[])

    repS_atoms    = [fo.canonical_representative() for fo in fo_list]
    repS_atom_set = set(repS_atoms)

    # ------------------------------------------------------------------ #
    # Step 1: Build row value lists.                                       #
    #                                                                      #
    # row_val_list[atom] is a sorted list of exactly n_cols values.        #
    # Phase 1 gives orbit_size distinct values; each is replicated         #
    # stab_u = n_cols / orbit_size times (orbit-stabiliser theorem).       #
    # ------------------------------------------------------------------ #
    row_val_list: dict = {}
    for fo, atom in zip(fo_list, repS_atoms):
        key    = (fo.operation.name, tuple(fo.flavour.counts))
        phase1 = phase1_results[key]
        orbit_size = len(phase1.values)
        if orbit_size == 0:
            row_val_list[atom] = []
            continue
        stab_u = n_cols // orbit_size          # = |stab_G(atom)|
        row_val_list[atom] = sorted(phase1.values * stab_u)

    # ------------------------------------------------------------------ #
    # Step 2: Build pair relations for cross-FO canonical pairs.           #
    #                                                                      #
    # pair_rel[(u_can, v_can)] = list of exactly n_cols (u_val, v_val)    #
    # tuples, one per group element.  Both (u,v) and (v,u) directions     #
    # are stored so the join loop only needs to look up (ref, new_atom).  #
    #                                                                      #
    # Filter: both atoms must be repS canonical representatives and        #
    # distinct.  Only then does (u_val, v_val) correspond to the SAME     #
    # column's entries for both rows.                                      #
    # ------------------------------------------------------------------ #
    pair_rel: dict = {}
    for decoded in all_pair_decoded:
        if not decoded.atom_pairs or not decoded.pairs:
            continue
        u0, v0 = decoded.atom_pairs[0]
        if u0 == v0 or u0 not in repS_atom_set or v0 not in repS_atom_set:
            continue
        stab_uv    = the_group.stabiliser_size_pair(u0, v0)
        full_pairs = list(decoded.pairs) * stab_uv
        if len(full_pairs) != n_cols:
            continue   # catches multi-orbit NULL_COMP or any other anomaly
        if (u0, v0) not in pair_rel:
            pair_rel[(u0, v0)] = full_pairs
        if (v0, u0) not in pair_rel:
            pair_rel[(v0, u0)] = [(y, x) for x, y in full_pairs]

    # ------------------------------------------------------------------ #
    # Step 3: Column-partition join.                                       #
    #                                                                      #
    # columns[j] maps atom → assigned value for column j.                 #
    # partition = list of groups; each group is a list of column indices   #
    # sharing identical values for all placed atoms so far.                #
    # ------------------------------------------------------------------ #

    # Process atoms in order of decreasing distinct-value count: the atom
    # with the most distinct values makes the finest initial partition,
    # maximising how much each subsequent pair constraint can contribute.
    ordered_atoms = sorted(
        repS_atoms,
        key=lambda a: -len(set(row_val_list.get(a, []))),
    )

    columns: list = [{} for _ in range(n_cols)]

    # Seed: assign sorted values of the first atom to columns 0..n_cols-1.
    first = ordered_atoms[0]
    for j, v in enumerate(row_val_list[first]):
        columns[j][first] = v

    placed: list    = [first]
    partition: list = _make_partition(columns, placed, n_cols)

    # Place remaining atoms one by one, refining the partition.
    for atom in ordered_atoms[1:]:
        new_partition: list  = []
        remaining:     Counter = Counter(row_val_list.get(atom, []))

        for group in partition:
            # All columns in this group share identical values for every
            # already-placed atom — use the first column as the representative.
            rep_col  = columns[group[0]]
            n_needed = len(group)

            # Gather v-value constraints from pair relations with placed atoms.
            candidates: Counter | None = None
            for ref in placed:
                key = (ref, atom)
                if key not in pair_rel:
                    continue
                ref_val = rep_col[ref]
                v_found = _lookup_v_vals(pair_rel[key], ref_val, atol)
                if not v_found:
                    continue
                c          = Counter(v_found)
                candidates = c if candidates is None else candidates & c

            # Determine the multiset of values to assign to this group.
            if candidates is None or sum(candidates.values()) == 0:
                # No usable constraint: free fill from remaining pool (sorted).
                group_vals = _take_n_smallest(remaining, n_needed)
            else:
                group_vals = sorted(candidates.elements())
                if len(group_vals) > n_needed:
                    group_vals = group_vals[:n_needed]
                elif len(group_vals) < n_needed:
                    # Constraint undershoots (collision or approximation error);
                    # pad from the remaining pool.
                    extra      = _take_n_smallest(remaining, n_needed - len(group_vals))
                    group_vals = sorted(group_vals + extra)

            # Deduct assigned values from the global remaining pool.
            for v in group_vals:
                remaining[v] -= 1

            # Assign values to columns and split group by new value.
            sub: dict = defaultdict(list)
            for col_idx, v in zip(group, group_vals):
                columns[col_idx][atom] = v
                sub[v].append(col_idx)
            new_partition.extend(sub.values())

        placed.append(atom)
        partition = new_partition

    # ------------------------------------------------------------------ #
    # Step 4: Assemble output multiset of column vectors.                  #
    # ------------------------------------------------------------------ #
    vectors = [
        tuple(col.get(a) for a in repS_atoms)
        for col in columns
    ]
    return AnnotatedMultisetOfRepSEvalVectors(vectors=vectors, atoms=repS_atoms)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _make_partition(columns: list, placed: list, n_cols: int) -> list:
    """Group column indices by their shared tuple of placed-atom values."""
    groups: dict = defaultdict(list)
    for j in range(n_cols):
        key = tuple(columns[j].get(a) for a in placed)
        groups[key].append(j)
    return list(groups.values())


def _lookup_v_vals(pairs: list, ref_val: float, atol: float) -> list:
    """Return second components of all pairs whose first component ≈ ref_val."""
    if atol == 0.0:
        return [v for u, v in pairs if u == ref_val]
    ref_f = float(ref_val)
    return [v for u, v in pairs if abs(float(u) - ref_f) <= atol]


def _take_n_smallest(pool: Counter, n: int) -> list:
    """Return (without removing) the n smallest values available in pool.

    Ties broken by sort order (stable).  Caller is responsible for updating
    pool after assigning the returned values.
    """
    if n <= 0:
        return []
    available = sorted(pool.elements())
    return available[:n]
