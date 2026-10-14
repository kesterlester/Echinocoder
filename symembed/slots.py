"""
Layer 1 (notes, Section 1): the slots of w and the signed-permutation action of G on them.

This module is pure Python (no torch, no sympy) and is deliberately an independent implementation of the definitions
in the notes, so that it can later be cross-checked against symatom.

Conventions
-----------
* Vectors are referred to by their position 0..n-1 along the input's n axis.
* A dot slot ``d`` has labels (i, j) with i <= j.  A determinant slot ``e`` has labels (i_1 < ... < i_k).
* For a group element g (a label permutation pi, optionally with parity), (g.w)_s is obtained by relabelling slot s
  by pi, restoring canonical order, and negating if s is an ``e`` and g contains parity:

      (g.w)[s] = sign[g][s] * w[index[g][s]]
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations, product
from typing import Iterator

from .spec import GroupSpec


@dataclass(frozen=True, order=True)
class Slot:
    """One entry of w: a dot product (kind ``"d"``) or a determinant (kind ``"e"``) of the vectors at ``positions``."""
    kind: str
    positions: tuple

    def name(self, labels: tuple) -> str:
        """Plain-text name, e.g. ``d_ab`` or ``e_abc`` (labels joined with commas if any is longer than one char)."""
        parts = [str(labels[i]) for i in self.positions]
        sep = "" if all(len(p) == 1 for p in parts) else ","
        return f"{self.kind}_{sep.join(parts)}"

    def latex(self, labels: tuple) -> str:
        parts = [str(labels[i]) for i in self.positions]
        sep = "" if all(len(p) == 1 for p in parts) else ","
        return f"{self.kind}_{{{sep.join(parts)}}}"

    def L_degree(self, k: int) -> int:
        """Degree as a polynomial in the components of the vectors."""
        return 2 if self.kind == "d" else k


@dataclass(frozen=True)
class GroupElement:
    """
    A label permutation (``perm[i]`` is the image of position i) together with a parity flag.
    """
    perm: tuple
    parity: bool

    def is_identity(self) -> bool:
        return not self.parity and all(i == p for i, p in enumerate(self.perm))

    def cycle_notation(self, labels: tuple) -> str:
        seen, cycles = set(), []
        for start in range(len(self.perm)):
            if start in seen:
                continue
            cycle, i = [], start
            while i not in seen:
                seen.add(i)
                cycle.append(i)
                i = self.perm[i]
            if len(cycle) > 1:
                cycles.append("(" + " ".join(str(labels[j]) for j in cycle) + ")")
        text = "".join(cycles) + ("p" if self.parity else "")
        return text or "e"


def _sort_with_sign(positions: tuple) -> tuple[tuple, int]:
    """Sort positions; return (sorted, sign of the sorting permutation). Positions are assumed distinct."""
    items = list(positions)
    sign = 1
    # Bubble sort: few elements, and each swap flips the sign.
    for i in range(len(items)):
        for j in range(len(items) - 1 - i):
            if items[j] > items[j + 1]:
                items[j], items[j + 1] = items[j + 1], items[j]
                sign = -sign
    return tuple(items), sign


def act(slot: Slot, perm: tuple, parity: bool) -> tuple[Slot, int]:
    """Image of a slot under a group element, as (canonical slot, sign)."""
    image = tuple(perm[i] for i in slot.positions)
    if slot.kind == "d":
        return Slot("d", tuple(sorted(image))), 1
    canonical, sign = _sort_with_sign(image)
    return Slot("e", canonical), (-sign if parity else sign)


def _generators(spec: GroupSpec) -> list[tuple[tuple, bool]]:
    """Adjacent transpositions within each species, plus parity if present. These generate G."""
    n = spec.n
    gens = []
    for block in spec.species_index_blocks:
        for a, b in zip(block, block[1:]):
            perm = list(range(n))
            perm[a], perm[b] = b, a
            gens.append((tuple(perm), False))
    if spec.parity:
        gens.append((tuple(range(n)), True))
    return gens


def build_slots(spec: GroupSpec, k: int) -> tuple:
    """
    The slots of w, grouped by G-orbit (signs ignored).

    Slots are first listed as: self-dots in position order, cross-dots (i < j) in lexicographic order, then
    determinants in lexicographic order. They are then stably sorted by orbit, orbits being ordered by the position
    of their first member in that listing. For L = (a, b, p), k = 2, G = S(2) x S(1) this gives

        d_aa, d_bb, d_pp, d_ab, d_ap, d_bp, e_ab, e_ap, e_bp

    as in the notes. Orbits are found from generators, without enumerating G.
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}")
    n = spec.n
    listing = ([Slot("d", (i, i)) for i in range(n)]
               + [Slot("d", (i, j)) for i, j in combinations(range(n), 2)]
               + [Slot("e", c) for c in combinations(range(n), k)])
    where = {s: idx for idx, s in enumerate(listing)}

    # Union-find over the listing, joined along the generators' action.
    parent = list(range(len(listing)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for perm, parity in _generators(spec):
        for idx, s in enumerate(listing):
            image, _ = act(s, perm, parity)
            ra, rb = find(idx), find(where[image])
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)

    first_member = {}
    for idx in range(len(listing)):
        first_member.setdefault(find(idx), idx)
    order = sorted(range(len(listing)), key=lambda idx: (first_member[find(idx)], idx))
    return tuple(listing[idx] for idx in order)


def slot_orbits(spec: GroupSpec, slots: tuple) -> tuple:
    """Consecutive runs of ``slots`` forming single G-orbits (as tuples of slot indices)."""
    where = {s: idx for idx, s in enumerate(slots)}
    parent = list(range(len(slots)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for perm, parity in _generators(spec):
        for idx, s in enumerate(slots):
            ra, rb = find(idx), find(where[act(s, perm, parity)[0]])
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    runs, current, current_root = [], [], None
    for idx in range(len(slots)):
        root = find(idx)
        if root != current_root and current:
            runs.append(tuple(current))
            current = []
        current_root = root
        current.append(idx)
    if current:
        runs.append(tuple(current))
    return tuple(runs)


def iter_group_elements(spec: GroupSpec) -> Iterator[GroupElement]:
    """All elements of G, identity first."""
    n = spec.n
    blocks = spec.species_index_blocks
    for parity in ((False, True) if spec.parity else (False,)):
        for block_images in product(*(permutations(b) for b in blocks)):
            perm = [0] * n
            for block, image in zip(blocks, block_images):
                for src, dst in zip(block, image):
                    perm[src] = dst
            yield GroupElement(tuple(perm), parity)


def signed_permutation_table(slots: tuple, elements: tuple) -> tuple[list, list]:
    """
    For each group element g and slot s: (g.w)[s] = sign[g][s] * w[index[g][s]].
    """
    where = {s: idx for idx, s in enumerate(slots)}
    index, sign = [], []
    for g in elements:
        row_index, row_sign = [], []
        for s in slots:
            image, sgn = act(s, g.perm, g.parity)
            row_index.append(where[image])
            row_sign.append(sgn)
        index.append(row_index)
        sign.append(row_sign)
    return index, sign
