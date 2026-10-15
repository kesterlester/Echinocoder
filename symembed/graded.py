"""
Graded packing (notes, Section 3.4): exploiting a dihedral subgroup <h, s> of G, s h s = h^-1, to kill most
coefficients of a complex Bursarial embedding of MS and to make the rest real.

Let h have order r and zeta = exp(2 pi i / r). Packed coordinates are complex linear functionals p of w with

    p(h.v) = zeta^j p(v)            when p is placed at power x^j   (grading),
    p(s.v) = conj(p(v))                                               (reality, if s is used).

Then q_{h.v}(x) = q_v(zeta x), and since h permutes MS, B(zeta x, y) = B(x, y): every coefficient E_{c,j} with
j != 0 (mod r) vanishes identically. With s, B equals its own conjugate, so the surviving coefficients are real.
The complex Bursarial embedding is injective on multisets of complex vectors (unique factorisation in C[x][y]) and the
packing is injective, so the surviving coefficients determine MS (Level 0). r = 1 (h the identity) is involution
packing.

Optionally the G-invariant linear functionals of w are output directly and left out of the packing. If h then acts on
every remaining coordinate by one scalar zeta^rho, the coordinates are placed consecutively and the grading is by c
instead ("scalar" mode): E_c = 0 unless r divides rho * c.

All coordinates are defined by explicit formulas (projections of slot functionals), so the packing is exact; it is
evaluated in floating point.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cyclotomic import Cyclotomic, field_for

TOL = 1e-9


# ----- linear algebra helpers ------------------------------------------------------------------------------------

def _real_rank(rows: list) -> int:
    if not rows:
        return 0
    A = np.array(rows)
    if np.iscomplexobj(A):
        A = np.hstack([A.real, A.imag])
    sv = np.linalg.svd(A, compute_uv=False)
    return int((sv > TOL * max(1.0, sv.max())).sum())


def _complex_rank(rows: list) -> int:
    if not rows:
        return 0
    sv = np.linalg.svd(np.array(rows), compute_uv=False)
    return int((sv > TOL * max(1.0, sv.max())).sum())


def _greedy(candidates, start=(), rank=_real_rank) -> list:
    """Candidates (in order) that increase the rank of ``start`` plus those already chosen."""
    chosen, current = [], list(start)
    base = rank(current)
    for c in candidates:
        if rank(current + [c]) > base:
            chosen.append(c)
            current.append(c)
            base += 1
    return chosen


def matrices(index: list, sign: list) -> list:
    """The signed-permutation matrices M_g with (g.w) = M_g w."""
    m = len(index[0])
    out = []
    for row_index, row_sign in zip(index, sign):
        M = np.zeros((m, m))
        M[np.arange(m), row_index] = row_sign
        out.append(M)
    return out


def _order(M) -> int:
    P, r, I = M.copy(), 1, np.eye(len(M))
    while not np.allclose(P, I):
        P, r = P @ M, r + 1
    return r


# ----- the packing -----------------------------------------------------------------------------------------------

@dataclass
class GradedPacking:
    """
    A graded packing of w. Position j of each q_g holds the coordinate plus[j] . v + i minus[j] . v (zero rows are
    padding). ``invariants`` rows are output directly.
    """
    h: int                      # index of h in the group-element list
    s: int | None               # index of s, or None (then coefficients are complex)
    r: int                      # order of h's action on w
    mode: str                   # "x": E_{c,j} = 0 unless r | j;  "scalar": E_c = 0 unless r | rho * c
    rho: int | None             # the common weight in scalar mode
    drop_trivial: bool
    invariants: np.ndarray      # (n_inv, m) real
    plus: np.ndarray            # (W, m) real
    minus: np.ndarray           # (W, m) real
    weights: list = field(default_factory=list)   # h-weight of each position, None for padding
    L: int = 4                                      # coefficients of the exact packing lie in Z[zeta_L]
    exact: np.ndarray | None = None                 # (W, m, f) exact coordinates (see build_packing)
    zeros: frozenset = frozenset()                  # (c, j) vanishing by character symmetry (symmetry_zeros)

    @property
    def real(self) -> bool:
        return self.s is not None

    @property
    def W(self) -> int:
        return self.plus.shape[0]

    def surviving(self, N: int) -> list[tuple[int, int]]:
        """(c, j) of the coefficients that are not forced to vanish, in order of c then j."""
        out = []
        if self.W == 0:
            return out
        # E_{c,j} is also zero, for every input, if j is not a sum of c occupied positions (one from each of c
        # factors of B): the product of c polynomials q_g has no other powers of x.
        occupied = [j for j, weight in enumerate(self.weights) if weight is not None]
        reachable = {0}
        for c in range(1, N + 1):
            reachable = {a + b for a in reachable for b in occupied}
            for j in range(c * (self.W - 1) + 1):
                if j not in reachable or (c, j) in self.zeros:
                    continue
                if self.mode == "x" and j % self.r == 0:
                    out.append((c, j))
                elif self.mode == "scalar" and (self.rho * c) % self.r == 0:
                    out.append((c, j))
        return out

    def output_count(self, N: int) -> int:
        return len(self.invariants) + (1 if self.real else 2) * len(self.surviving(N))


def _int_matrix(M) -> np.ndarray:
    return np.array(np.rint(M).astype(int), dtype=object)


def build_packing(Ms: list, h: int, s: int | None, drop_trivial: bool) -> GradedPacking:
    """
    Build the packing exactly, with coefficients in Z[zeta_L] (L = lcm(r, 4)); the floating-point coordinates used
    for evaluation are derived from the exact ones. An exact row is an object array of shape (m, f): entry t is the
    coefficient of w_t, as an integer vector in the power basis of Z[zeta_L]. Rank decisions (which candidates to
    keep) use floating point; the rows kept are exact.
    """
    m = Ms[0].shape[0]
    I = np.eye(m)
    Mh = Ms[h]
    Mss = None if s is None else Ms[s]
    r = _order(Mh)
    F = field_for(r)
    step = F.L // r                                                   # zeta_r = zeta_L^step
    powers = [np.linalg.matrix_power(Mh, t) for t in range(r)]
    int_powers = [_int_matrix(P) for P in powers]

    def embed(row) -> np.ndarray:                                     # integer row -> exact row
        U = np.zeros((m, F.f), dtype=object)
        U[:, 0] = [int(round(x)) for x in row]
        return U

    def field_map(A, U) -> np.ndarray:                                # apply a field-linear map to every entry
        return U @ A.T

    def times_matrix(U, M) -> np.ndarray:                             # the row (functional) U M
        return _int_matrix(M).T @ U

    def tau(U) -> np.ndarray:                                         # tau(p) = conj(p M_s)
        return field_map(F.conj_matrix, times_matrix(U, Mss))

    def greedy(exact_rows, start=(), rank=_real_rank) -> list:
        floats = [F.to_complex(U) for U in exact_rows]
        chosen, current = [], [np.asarray(x, dtype=complex) for x in start]
        base = rank(current)
        for U, x in zip(exact_rows, floats):
            if rank(current + [x]) > base:
                chosen.append(U)
                current.append(x)
                base += 1
        return chosen

    # Invariant functionals: rows of S = sum_g M_g (exact integers; S/N is the group average).
    invariants = _greedy(list(np.rint(sum(Ms)))) if drop_trivial else []
    i_matrix = F.mul_matrix(F.i)

    by_residue: dict[int, list] = {}
    for k in range(r):
        if k > r - k:
            continue
        if k == 0 or 2 * k == r:
            # Real weight space: the exact integer matrix r P_k = sum_t (+-1)^t M_h^t has the right row space.
            P = np.rint(sum(((-1) ** t if k else 1) * powers[t] for t in range(r)))
            if Mss is None:
                rows = [embed(x) for x in _greedy(list(P), start=invariants if k == 0 else ())]
                half = (len(rows) + 1) // 2
                even, odd = rows[:half], rows[half:]
            else:
                even = [embed(x) for x in _greedy(list(P @ (I + Mss)), start=invariants if k == 0 else ())]
                odd = [embed(x) for x in _greedy(list(P @ (I - Mss)))]
            n = max(len(even), len(odd))
            zero = np.zeros((m, F.f), dtype=object)
            even = even + [zero] * (n - len(even))
            odd = odd + [zero] * (n - len(odd))
            by_residue.setdefault(k, []).extend(a + field_map(i_matrix, b) for a, b in zip(even, odd))
        else:
            # Complex weight space: rows of r P_k = sum_t zeta_r^(-k t) M_h^t, exactly.
            rows = []
            for slot in range(m):
                U = np.zeros((m, F.f), dtype=object)
                for t in range(r):
                    U = U + np.outer(int_powers[t][slot], F.zeta_power(-k * t * step))
                rows.append(U)
            if Mss is None:
                basis = greedy(rows, rank=_complex_rank)
            else:
                candidates = [x for U in rows for x in (U + tau(U), field_map(i_matrix, U - tau(U)))]
                basis = greedy(candidates)
            for U in basis:                                               # balance residues k and r - k
                if len(by_residue.get(k, [])) <= len(by_residue.get(r - k, [])):
                    by_residue.setdefault(k, []).append(U)
                else:
                    by_residue.setdefault(r - k, []).append(field_map(F.conj_matrix, U))

    used = sorted(rho for rho, coords in by_residue.items() if coords)
    if not used:                                                      # everything is invariant: nothing to pack
        mode, rho, weights, exact = "none", None, [], np.zeros((0, m, F.f), dtype=object)
    elif len(used) == 1 and used[0] != 0:
        mode, rho = "scalar", used[0]
        coords = by_residue[rho]
        weights = [rho] * len(coords)
        exact = np.array(coords, dtype=object).reshape(len(coords), m, F.f)
    else:
        mode, rho = "x", None
        W = max(res + r * (len(coords) - 1) for res, coords in by_residue.items() if coords) + 1
        exact = np.zeros((W, m, F.f), dtype=object)
        weights = [None] * W
        for res, coords in by_residue.items():
            for t, U in enumerate(coords):
                exact[res + r * t] = U
                weights[res + r * t] = res
    Phi = F.to_complex(exact) if len(exact) else np.zeros((0, m), dtype=complex)
    inv = np.array(invariants).reshape(len(invariants), m)
    return GradedPacking(h, s, r, mode, rho, drop_trivial, inv, Phi.real.copy(), Phi.imag.copy(), weights,
                         L=F.L, exact=exact)


def check_packing_exact(Ms: list, packing: GradedPacking) -> None:
    """Exact check, in Z[zeta_L], of the grading and reality properties of every packed coordinate."""
    F = Cyclotomic(packing.L)
    step = F.L // packing.r
    for U, weight in zip(packing.exact, packing.weights):
        if weight is None:
            assert not np.any(U)
            continue
        lhs = _int_matrix(Ms[packing.h]).T @ U
        rhs = U @ F.mul_matrix(F.zeta_power(weight * step)).T
        assert np.array_equal(lhs, rhs), "grading fails exactly"
        if packing.s is not None:
            assert np.array_equal((_int_matrix(Ms[packing.s]).T @ U) @ F.conj_matrix.T, U), "reality fails exactly"


def check_packing(Ms: list, packing: GradedPacking) -> None:
    """Numerical check of the defining properties: grading, reality, injectivity."""
    m = Ms[0].shape[0]
    Phi = packing.plus + 1j * packing.minus
    zeta = np.exp(2j * np.pi / packing.r)
    for j, wt in enumerate(packing.weights):
        if wt is None:
            assert np.allclose(Phi[j], 0)
            continue
        assert np.allclose(Phi[j] @ Ms[packing.h], zeta ** wt * Phi[j], atol=1e-9), "grading fails"
        if packing.s is not None:
            assert np.allclose(Phi[j] @ Ms[packing.s], np.conj(Phi[j]), atol=1e-9), "reality fails"
    for v in packing.invariants:
        assert all(np.allclose(v @ M, v, atol=1e-9) for M in Ms), "not invariant"
    stacked = np.vstack([packing.plus, packing.minus, packing.invariants])
    assert np.linalg.matrix_rank(stacked, tol=1e-8) == m, "packing not injective"


def conjugacy_key(element, species_blocks) -> tuple:
    """Cycle type within each species, and parity: conjugate (h, s) pairs give identical output counts."""
    key = []
    for block in species_blocks:
        seen, lengths = set(), []
        for start in block:
            if start in seen:
                continue
            length, i = 0, start
            while i not in seen:
                seen.add(i)
                i = element.perm[i]
                length += 1
            lengths.append(length)
        key.append(tuple(sorted(lengths)))
    return tuple(key), element.parity


def choose_packing(index: list, sign: list, elements: tuple, species_blocks) -> GradedPacking:
    """
    Search (h, s, drop_trivial) for the fewest outputs: one h per conjugacy class, every involution s with
    s h s = h^-1 (and s = None). Ties go to the first found.
    """
    Ms = matrices(index, sign)
    N, m = len(Ms), Ms[0].shape[0]
    I = np.eye(m)
    involutions = [g for g in range(N) if np.allclose(Ms[g] @ Ms[g], I)]
    best, best_count, seen = None, None, set()
    for h, element in enumerate(elements):
        key = conjugacy_key(element, species_blocks)
        if key in seen:
            continue
        seen.add(key)
        Mh_inv = Ms[h].T                                               # signed permutations are orthogonal
        partners = [s for s in involutions if np.allclose(Ms[s] @ Ms[h] @ Ms[s], Mh_inv)] + [None]
        for s in partners:
            for drop in (False, True):
                packing = build_packing(Ms, h, s, drop)
                count = packing.output_count(N)
                if best_count is None or count < best_count:
                    best, best_count = packing, count
    check_packing(Ms, best)
    check_packing_exact(Ms, best)
    best.zeros = symmetry_zeros(best, diagonal_characters(Ms, best), N)
    return best


# ----- the evaluator ---------------------------------------------------------------------------------------------

import torch                                                             # noqa: E402

from .bursar import bursar_coefficients_complex, bursar_output_labels  # noqa: E402
from .module import InvariantEmbedding                                   # noqa: E402


def output_selection(packing: GradedPacking, N: int) -> list[tuple]:
    """
    The outputs after the invariants, in order: for each surviving (c, j), its real part and (if the coefficients
    are complex) its imaginary part. Each entry is (c, j, "re" | "im", index into the concatenated [re, im] array).
    """
    if packing.W == 0:
        return []
    position = {label: i for i, label in enumerate(bursar_output_labels(N, packing.W))}
    D_full = len(position)
    out = []
    for c, j in packing.surviving(N):
        out.append((c, j, "re", position[(c, j)]))
        if not packing.real:
            out.append((c, j, "im", D_full + position[(c, j)]))
    return out


class GradedEmbedding(InvariantEmbedding):
    """
    The graded-packing evaluator: G-invariant coordinates of w (if any), then the surviving coefficients of the
    complex Bursarial polynomial of the packed multiset {phi(g.w)}.
    """

    def __init__(self, n: int, k: int, slots: tuple, index: list, sign: list, dtype: torch.dtype,
                 packing: GradedPacking, keep: list | None = None, metric: tuple | None = None):
        super().__init__(n, k, slots, index, sign, dtype, metric=metric)
        selection = output_selection(packing, len(index))
        total = len(packing.invariants) + len(selection)
        keep = list(range(total)) if keep is None else list(keep)
        self.register_buffer("inv", torch.tensor(packing.invariants, dtype=dtype).reshape(-1, self.m))
        self.register_buffer("graded_plus", torch.tensor(packing.plus, dtype=dtype))
        self.register_buffer("graded_minus", torch.tensor(packing.minus, dtype=dtype))
        self.register_buffer("select", torch.tensor([t[3] for t in selection], dtype=torch.int64))
        self.register_buffer("graded_keep", torch.tensor(keep, dtype=torch.int64))
        self.graded_keeps_all = len(keep) == total
        self.output_shape = torch.tensor([len(keep)], dtype=torch.int64)

    def forward(self, x) -> torch.Tensor:
        x = self._prepare(x)
        batch_shape = x.shape[:-2]
        out = self.from_w(self.w(x).reshape(-1, self.m))
        return out.reshape(*batch_shape, out.shape[-1])

    def from_w(self, w: torch.Tensor) -> torch.Tensor:
        """The outputs for given slot values w, shape (B, m) -> (B, D)."""
        if self.graded_plus.shape[0] == 0:
            out = w @ self.inv.T
        else:
            gw = w[:, self.g_index] * self.g_sign
            re, im = bursar_coefficients_complex(gw @ self.graded_plus.T, gw @ self.graded_minus.T)
            out = torch.cat([re, im], -1)[:, self.select]
            if self.inv.shape[0]:
                out = torch.cat([w @ self.inv.T, out], -1)
        if not self.graded_keeps_all:
            out = out[:, self.graded_keep]
        return out


# ----- metadata --------------------------------------------------------------------------------------------------

import json                                                              # noqa: E402

from .meta import SYMEMBED_VERSION, CERTIFICATION_E_DROPPED, _echinocoder_commit, _latex_text, possible_L_degrees  # noqa: E402,E501
from .slots import slot_orbits                                           # noqa: E402

CERTIFICATION_GRADED = (
    "Level 0 (by construction): graded packing. Packed coordinates p satisfy p(h.v) = zeta^j p(v) at power x^j "
    "(zeta = exp(2 pi i/r)) and, when s is used, p(s.v) = conj p(v). As h and s permute MS, B(zeta x, y) = B(x, y) "
    "and conj B = B, so the omitted coefficients vanish identically and the kept ones are real (complex if no s). The "
    "complex Bursarial embedding is injective on multisets of complex vectors (unique factorisation in C[x][y]); the "
    "packing, together with the G-invariant coordinates output directly, is injective; so the outputs determine MS. MS "
    "determines L up to G and rotations (up to G and O(k) if parity is in G)."
)


def _form(row, names, digits: int = 4) -> str:
    """A real linear functional as text, e.g. '0.5 d_aa + 0.5 d_bb'."""
    terms = []
    for coeff, name in zip(row, names):
        if abs(coeff) < 1e-12:
            continue
        mag = abs(coeff)
        text = name if abs(mag - 1) < 1e-12 else f"{mag:.{digits}g} {name}"
        terms.append(("- " if coeff < 0 else "+ ") + text)
    out = " ".join(terms)
    return "0" if not out else (out[2:] if out.startswith("+ ") else "-" + out[2:])


class GradedMeta:
    """Metadata for embedder="bursar_graded" (same reading API as EmbeddingMeta: D, to_dict, to_json, latex)."""

    def __init__(self, group_spec, k, options, slots, elements, index, sign, packing: GradedPacking, e_dropped: bool,
                 relations: tuple = (), kept: list | None = None, linear_reduction_max_c: int = 1):
        self.group_spec, self.k, self.options, self.slots = group_spec, k, options, slots
        self.elements, self.index, self.sign, self.packing, self.e_dropped = elements, index, sign, packing, e_dropped
        self.relations = tuple(relations)          # LinearRelation, indices into the full output list
        total = len(packing.invariants) + len(output_selection(packing, len(elements)))
        self.kept = list(range(total)) if kept is None else list(kept)
        self.linear_reduction_max_c = linear_reduction_max_c
        self.echinocoder_commit = _echinocoder_commit()

    labels = property(lambda self: self.group_spec.labels)
    n = property(lambda self: self.group_spec.n)
    m = property(lambda self: len(self.slots))
    N = property(lambda self: len(self.elements))
    @property
    def certification_level(self) -> int:
        return 1 if self.relations else 0

    def slot_names(self) -> list[str]:
        return [s.name(self.labels) for s in self.slots]

    def selection(self) -> list[tuple]:
        return output_selection(self.packing, self.N)

    @property
    def D(self) -> int:
        return len(self.kept)

    def all_names(self) -> list[str]:
        """Names of the full output list (before any Level 1 reduction)."""
        names = [f"inv[{t}]" for t in range(len(self.packing.invariants))]
        names += [f"{'Re' if part == 're' else 'Im'} E_{{{c},{j}}}" for c, j, part, _ in self.selection()]
        return names

    def output_names(self) -> list[str]:
        names = self.all_names()
        return [names[i] for i in self.kept]

    def relation_summary(self, rel) -> str:
        names = self.all_names()
        if not rel.terms:
            return f"{names[rel.dropped]} = 0"
        pieces = []
        for t, coeff in rel.terms:
            text = coefficient_text(coeff, getattr(rel, "L", None))
            compound = any(ch in text[1:] for ch in "+-") or " " in text
            coefficient = "" if text == "1" else ("-" if text == "-1" else (f"({text}) " if compound else f"{text} "))
            pieces.append(f"{coefficient}{names[t]}")
        return f"{names[rel.dropped]} = " + " + ".join(pieces).replace("+ -", "- ")

    @property
    def certification(self) -> str:
        text = CERTIFICATION_GRADED + (" " + CERTIFICATION_E_DROPPED if self.e_dropped else "")
        if self.relations:
            text = ("Level 1 (explicit identity): each dropped output is an exact linear combination of kept outputs, "
                    "with rational coefficients (c = 1) or coefficients in Q(zeta_L) (c >= 2), verified in exact "
                    "integer arithmetic. " + text)
        return text

    def element_name(self, g) -> str | None:
        return None if g is None else self.elements[g].cycle_notation(self.labels)

    def coordinates(self) -> list[dict]:
        names = self.slot_names()
        return [{"position": j, "weight": wt, "re": _form(p, names), "im": _form(q, names)}
                for j, (wt, p, q) in enumerate(zip(self.packing.weights, self.packing.plus, self.packing.minus))]

    def reductions_log(self) -> list[dict]:
        log = []
        if self.e_dropped:
            log.append({"type": "omit_determinants_with_parity", "level": 0, "reason": CERTIFICATION_E_DROPPED})
        if self.options.reductions == "safe":
            log.append({"type": "linear_scope", "level": 1,
                        "reason": f"exact linear search: the c = 1 block (integer coefficients), and 2 <= c <= "
                                  f"{self.linear_reduction_max_c} over Q(zeta_{self.packing.L}) by restriction of "
                                  f"scalars; outputs with larger c are all kept"})
        names = self.all_names()
        for rel in self.relations:
            log.append({"type": "linear", "level": 1, "dropped": names[rel.dropped],
                        "terms": [{"name": names[t], "coefficient": coefficient_text(f, getattr(rel, "L", None))}
                                  for t, f in rel.terms],
                        "summary": self.relation_summary(rel), "verified": "exact (integer arithmetic)"})
        return log

    def reconstruct_base(self, Y):
        """The full output list (before Level 1 reduction), recomputed from the kept outputs Y."""
        total = len(self.all_names())
        columns = [None] * total
        for out, b in enumerate(self.kept):
            columns[b] = Y[..., out]
        for rel in self.relations:
            value = 0 * Y[..., 0]
            for t, coeff in rel.terms:
                value = value + columns[t] * coefficient_value(coeff, getattr(rel, "L", None))
            columns[rel.dropped] = value
        stack = torch.stack if isinstance(Y, torch.Tensor) else np.stack
        return stack(columns, -1)

    def possible_degrees(self) -> list[list[int]]:
        inv = [sorted({self.slots[t].L_degree(self.k) for t in range(self.m) if abs(v[t]) > 1e-12})
               for v in self.packing.invariants]
        if self.packing.W == 0:
            return [inv[i] for i in self.kept]
        kinds = [{self.slots[t].kind for t in range(self.m) if abs(p[t]) > 1e-12 or abs(q[t]) > 1e-12}
                 for p, q in zip(self.packing.plus, self.packing.minus)]
        full = possible_L_degrees(self.slots, self.k, self.N, kinds)
        position = {lab: i for i, lab in enumerate(bursar_output_labels(self.N, self.packing.W))}
        alldeg = inv + [full[position[(c, j)]] for c, j, _, _ in self.selection()]
        return [alldeg[i] for i in self.kept]

    def to_dict(self) -> dict:
        p = self.packing
        degrees = self.possible_degrees()
        return {
            "symembed_version": SYMEMBED_VERSION, "echinocoder_commit": self.echinocoder_commit,
            "group_spec": self.group_spec.to_dict(), "k": self.k, "options": self.options.to_dict(),
            "n": self.n, "m": self.m, "multiset_size": self.N, "D": self.D, "base_D": len(self.all_names()),
            "certification_level": self.certification_level, "certification": self.certification,
            "reductions": self.reductions_log(),
            "slots": [{"name": s.name(self.labels), "kind": s.kind, "L_degree": s.L_degree(self.k)} for s in self.slots],
            "slot_orbits": [list(run) for run in slot_orbits(self.group_spec, self.slots)],
            "graded": {"h": self.element_name(p.h), "s": self.element_name(p.s), "r": p.r, "mode": p.mode,
                       "rho": p.rho, "real_coefficients": p.real, "drop_trivial": p.drop_trivial,
                       "invariants": [_form(v, self.slot_names()) for v in p.invariants],
                       "coordinates": self.coordinates(),
                       "meaning": "q_g(x) = sum_j (re_j + i im_j)(g.w) x^j; outputs: invariants, then surviving "
                                  "coefficients of prod_g (y + q_g(x)) (real parts, and imaginary parts if complex)"},
            "outputs": [{"index": i, "name": name, "possible_L_degrees": deg}
                        for i, (name, deg) in enumerate(zip(self.output_names(), degrees))],
        }

    def to_json(self, **kwargs) -> str:
        return json.dumps(self.to_dict(), **kwargs)

    def latex(self, max_rows: int = 200) -> str:
        p, L = self.packing, self.labels
        spec = " \\times ".join(f"S({len(s)})" for s in self.group_spec.species)
        if self.group_spec.parity:
            spec += " \\times P"

        def tex(form):
            import re
            return re.sub(r"([de])_([A-Za-z0-9,]+)", r"\1_{\2}", form)

        s_name = self.element_name(p.s) or "\\text{none}"
        parts = [f"\\noindent Graded packing of $n = {self.n}$ vectors in $k = {self.k}$ dimensions, $G = {spec}$, "
                 f"multiset size $N = {self.N}$, $m = {self.m}$, $D = {self.D}$ outputs. "
                 f"$h = {self.element_name(p.h)}$ (order $r = {p.r}$), $s = {s_name}$, "
                 f"mode {p.mode}. Certification: {_latex_text(self.certification)}\n"]
        if len(p.invariants):
            rows = [f"{t} & ${tex(_form(v, self.slot_names()))}$ \\\\" for t, v in enumerate(p.invariants)]
            parts.append("\\begin{longtable}{rl}\n$t$ & invariant output \\\\ \\hline\n" + "\n".join(rows)
                         + "\n\\end{longtable}\n")
        rows = [f"{c['position']} & {'' if c['weight'] is None else c['weight']} & ${tex(c['re'])}$ & ${tex(c['im'])}$ \\\\"
                for c in self.coordinates()]
        parts.append("\\begin{longtable}{rrll}\n$j$ & weight & real part & imaginary part \\\\ \\hline\n"
                     + "\n".join(rows) + "\n\\end{longtable}\n")
        names, degrees = self.output_names(), self.possible_degrees()
        rows = [f"{i} & ${name}$ & {', '.join(map(str, deg))} \\\\"
                for i, (name, deg) in enumerate(zip(names[:max_rows], degrees))]
        if len(names) > max_rows:
            rows.append(f"\\multicolumn{{3}}{{l}}{{\\dots\\ {len(names) - max_rows} more outputs}} \\\\")
        parts.append("\\begin{longtable}{rll}\noutput & & possible degrees in $L$ \\\\ \\hline\n" + "\n".join(rows)
                     + "\n\\end{longtable}\n")
        return "\n".join(parts)


# ----- exact Level 1 reduction of the c = 1 block ----------------------------------------------------------------

from .reductions import LinearRelation, find_linear_relations          # noqa: E402


def c1_relations(Ms: list, packing: GradedPacking, N: int) -> list:
    """
    Exact linear relations involving the c = 1 outputs.

    E_{1,j} = sum_g p_j(g.w) = (p_j S) w with S = sum_g M_g, a linear function of w with integer coefficients
    (kept powers j are multiples of r, which hold weight-0 coordinates, built from integer rows). The invariant
    outputs are the integer rows of S. Any c = 1 output that is an exact rational combination of the invariant outputs
    and earlier c = 1 outputs is dropped (Level 1). Indices refer to the full output list: invariants, then the
    coefficients of ``output_selection``.
    """
    S = np.rint(sum(Ms))
    n_inv = len(packing.invariants)
    rows, owners = [np.rint(v) for v in packing.invariants], list(range(n_inv))
    for t, (c, j, part, _) in enumerate(output_selection(packing, N)):
        if c != 1:
            continue
        row = (packing.plus[j] if part == "re" else packing.minus[j]) @ S
        if not np.allclose(row, np.rint(row)):
            raise AssertionError("c = 1 output does not have integer coefficients")
        rows.append(np.rint(row))
        owners.append(n_inv + t)
    if len(rows) == n_inv:
        return []
    matrix = np.array([[int(x) for x in row] for row in rows], dtype=object)
    _, relations = find_linear_relations(matrix)
    return [LinearRelation(owners[rel.dropped], tuple((owners[t], f) for t, f in rel.terms)) for rel in relations]


# ----- exact Level 1 reduction for c >= 2, over Z[zeta_L] --------------------------------------------------------

from fractions import Fraction                                          # noqa: E402

from .reductions import choose_c_max, exact_expansion_field             # noqa: E402


def packed_field_forms(index: list, sign: list, packing: GradedPacking) -> list:
    """forms[g][j]: coefficient j of q_g as (slot t, a in Z[zeta_L]) pairs, from the exact packing."""
    f = packing.exact.shape[2]
    forms = []
    for row_index, row_sign in zip(index, sign):
        row = []
        for U in packing.exact:
            terms: dict[int, np.ndarray] = {}
            for s in range(U.shape[0]):
                if np.any(U[s]):
                    t = row_index[s]
                    terms[t] = terms.get(t, np.zeros(f, dtype=object)) + row_sign[s] * U[s]
            row.append([(t, a) for t, a in sorted(terms.items()) if np.any(a)])
        forms.append(row)
    return forms


@dataclass(frozen=True)
class FieldRelation:
    """Output ``dropped`` = sum over terms (t, lam) of lam * output t, lam in Q(zeta_L) (power basis, Fractions)."""
    dropped: int
    terms: tuple
    L: int


def coefficient_value(coeff, L: int | None = None) -> float:
    """Numerical value of a relation coefficient: a Fraction, or a power-basis tuple of Fractions in Q(zeta_L)."""
    if isinstance(coeff, Fraction):
        return coeff.numerator / coeff.denominator
    value = sum((c.numerator / c.denominator) * np.exp(2j * np.pi * a / L) for a, c in enumerate(coeff))
    return float(np.real(value))


def coefficient_text(coeff, L: int | None = None) -> str:
    if isinstance(coeff, Fraction):
        return str(coeff)
    import sympy
    z = sympy.exp(2 * sympy.pi * sympy.I / L)
    expr = sum(sympy.Rational(c.numerator, c.denominator) * z ** a for a, c in enumerate(coeff))
    return str(sympy.nsimplify(sympy.re(sympy.expand_complex(expr))))


def field_relations(index: list, sign: list, packing: GradedPacking, N: int, kept: list, budget: int):
    """
    Exact linear relations among kept outputs with 2 <= c <= c_max, with coefficients in Q(zeta_L).

    The outputs are expanded exactly in Z[zeta_L] (c_max chosen from ``budget``). First, as an exact check of the
    construction, every omitted coefficient is verified to vanish and (with s) every kept one to equal its conjugate.
    Then, by "restriction of scalars", each kept output P contributes the rows P, zeta P, ..., zeta^(f-1) P (as integer
    vectors). Q-linear relations among these rows are exactly Q(zeta_L)-linear relations among the outputs, so the
    integer relation search applies unchanged: an output is dropped when its first row depends on rows of earlier
    kept outputs. Returns (relations, c_max).
    """
    F = Cyclotomic(packing.L)
    m = packing.exact.shape[1]
    c_max = min(choose_c_max(N, m, max(1, budget // F.f), packing.W), N)
    if c_max < 1 or packing.W == 0:
        return [], c_max
    expansions = exact_expansion_field(packed_field_forms(index, sign, packing), m, c_max, F)
    conj = F.conj_matrix
    minus_i = F.mul_matrix(F.conj(F.i)) if F.L % 4 == 0 else None
    zeta = F.mul_matrix(F.zeta_power(1))
    surviving = set(packing.surviving(N))
    for c, (_, coeffs) in enumerate(expansions, start=1):
        for j in range(coeffs.shape[0]):
            V = coeffs[j].astype(object)
            if (c, j) not in surviving:
                assert not np.any(V), f"E_{{{c},{j}}} should vanish identically but does not"
            elif packing.real:
                assert np.array_equal(V @ conj.T, V), f"E_{{{c},{j}}} should be real but is not"
    n_inv = len(packing.invariants)
    selection = output_selection(packing, N)
    relations = []
    for c in range(2, c_max + 1):
        owners, blocks = [], []
        for idx in kept:
            if idx < n_inv:
                continue
            cc, j, part, _ = selection[idx - n_inv]
            if cc != c:
                continue
            V = expansions[c - 1][1][j].astype(object)
            if not packing.real:
                V = V + V @ conj.T if part == "re" else (V - V @ conj.T) @ minus_i.T     # 2 Re, 2 Im
            rows, P = [], V
            for _ in range(F.f):
                rows.append(P.reshape(-1))
                P = P @ zeta.T
            owners.append(idx)
            blocks.append(rows)
        if len(owners) < 2:
            continue
        matrix = np.array([row for rows in blocks for row in rows], dtype=object)
        _, rels = find_linear_relations(matrix)
        for rel in rels:
            o, a0 = divmod(rel.dropped, F.f)
            if a0:
                continue                                              # a multiple of a dropped output
            lam: dict[int, list] = {}
            for t, frac in rel.terms:
                ot, a = divmod(t, F.f)
                lam.setdefault(ot, [Fraction(0)] * F.f)[a] += frac
            relations.append(FieldRelation(owners[o], tuple((owners[ot], tuple(v)) for ot, v in lam.items()), F.L))
    return relations, c_max


# ----- exact zeros from sign (character) symmetries of individual coordinates ------------------------------------

def diagonal_characters(Ms: list, packing: GradedPacking) -> dict:
    """
    For each group element t that multiplies some packed coordinate by a non-trivial scalar, the map
    {position a: e}, meaning phi_a(t.v) = zeta_L^e phi_a(v) for all v. Found exactly in Z[zeta_L].
    """
    F = Cyclotomic(packing.L)
    scaled = {}
    for a, U in enumerate(packing.exact):
        if packing.weights[a] is None:
            continue
        candidates = {e: U @ F.mul_matrix(F.zeta_power(e)).T for e in range(F.L)}
        for t, M in enumerate(Ms):
            V = _int_matrix(M).T @ U
            for e, Z in candidates.items():
                if np.array_equal(V, Z):
                    scaled.setdefault(t, {})[a] = e
                    break
    return {t: chars for t, chars in scaled.items() if any(e for e in chars.values())}


def _patterns(deficits: list, d: int, max_parts: int, start: int = 0):
    """Multisets (as lists) of entries of ``deficits`` (positive, ascending) summing to d with at most max_parts parts."""
    if d == 0:
        yield []
        return
    if max_parts == 0:
        return
    for i in range(start, len(deficits)):
        x = deficits[i]
        if x > d:
            break
        for rest in _patterns(deficits, d - x, max_parts - 1, i):
            yield [x] + rest


def symmetry_zeros(packing: GradedPacking, characters: dict, N: int, max_deficit: int = 12) -> frozenset:
    """
    Coefficients E_{c,j} that vanish for every input because every position pattern contributing to them is killed by
    a character symmetry (see the notes, Section 3.4). Only coefficients within ``max_deficit`` of the bottom or top
    end are examined; there the patterns are few (most positions equal the end position).
    """
    if packing.W == 0 or not characters:
        return frozenset()
    occupied = [j for j, wt in enumerate(packing.weights) if wt is not None]
    lo, hi = occupied[0], occupied[-1]
    L = packing.L

    def killed(multiplicity: dict) -> bool:
        for chars in characters.values():
            if all(a in chars for a in multiplicity) and sum(n * chars[a] for a, n in multiplicity.items()) % L:
                return True
        return False

    def all_patterns_killed(c: int, end: int, d: int, deficit_of) -> bool:
        """True if E_{c,j} (at deficit d from ``end``) has at least one pattern and every pattern is killed."""
        by_deficit = {deficit_of(a): a for a in occupied if a != end}
        deficits = sorted(x for x in by_deficit if 0 < x <= d)
        found = False
        for pattern in _patterns(deficits, d, c):
            found = True
            multiplicity = {end: c - len(pattern)} if c > len(pattern) else {}
            for x in pattern:
                multiplicity[by_deficit[x]] = multiplicity.get(by_deficit[x], 0) + 1
            if not killed(multiplicity):
                return False
        return found

    zeros = set()
    for c, j in packing.surviving(N):
        bottom, top = j - c * lo, c * hi - j
        if (bottom <= max_deficit and all_patterns_killed(c, lo, bottom, lambda a: a - lo)) or \
                (top <= max_deficit and all_patterns_killed(c, hi, top, lambda a: hi - a)):
            zeros.add((c, j))
    return frozenset(zeros)


# ----- exact symbolic expansion of graded outputs (small cases) --------------------------------------------------

def _graded_expanded_outputs(meta, max_group_order: int = 24) -> list:
    """
    Every (kept) output of a graded embedding, fully expanded as a polynomial in the d's and e's, with exact
    coefficients (radicals where possible). Feasible only for small groups.
    """
    import sympy
    if meta.N > max_group_order:
        raise ValueError(f"multiset size {meta.N} exceeds max_group_order = {max_group_order}")
    p = meta.packing
    w = [sympy.Symbol(nm) for nm in meta.slot_names()]
    out_all = [sum((int(round(c)) * wi for c, wi in zip(row, w) if abs(c) > 1e-12), sympy.Integer(0))
               for row in p.invariants]
    if p.W:
        F = Cyclotomic(p.L)
        expansions = exact_expansion_field(packed_field_forms(meta.index, meta.sign, p), meta.m, meta.N, F)
        conj = F.conj_matrix
        minus_i = F.mul_matrix(F.conj(F.i))
        cache: dict = {}

        def value(vec) -> sympy.Expr:
            key = tuple(int(x) for x in vec)
            if key not in cache:
                cache[key] = F.to_sympy(np.array(key, dtype=object))
            return cache[key]

        for c, j, part, _ in output_selection(p, meta.N):
            monos, coeffs = expansions[c - 1]
            V = coeffs[j].astype(object)
            scale = sympy.Integer(1)
            if not p.real:
                V = V + V @ conj.T if part == "re" else (V - V @ conj.T) @ minus_i.T
                scale = sympy.Rational(1, 2)
            expr = sympy.Integer(0)
            for mono, vec in zip(monos, V):
                if np.any(vec):
                    expr += value(vec) * sympy.Mul(*[w[s] for s in mono])
            out_all.append(sympy.expand(scale * expr))
    return [out_all[i] for i in meta.kept]


GradedMeta.expanded_outputs = _graded_expanded_outputs
