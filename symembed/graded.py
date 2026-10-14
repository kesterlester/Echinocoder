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
        for c in range(1, N + 1):
            for j in range(c * (self.W - 1) + 1):
                if self.mode == "x" and j % self.r == 0:
                    out.append((c, j))
                elif self.mode == "scalar" and (self.rho * c) % self.r == 0:
                    out.append((c, j))
        return out

    def output_count(self, N: int) -> int:
        return len(self.invariants) + (1 if self.real else 2) * len(self.surviving(N))


def build_packing(Ms: list, h: int, s: int | None, drop_trivial: bool) -> GradedPacking:
    m = Ms[0].shape[0]
    I = np.eye(m)
    Mh = Ms[h]
    Mss = None if s is None else Ms[s]
    r = _order(Mh)
    zeta = np.exp(2j * np.pi / r)
    powers = [np.linalg.matrix_power(Mh, t) for t in range(r)]
    invariants = _greedy(list(sum(Ms) / len(Ms))) if drop_trivial else []

    by_residue: dict[int, list] = {}
    for k in range(r):
        if k > r - k:
            continue
        P = sum(zeta ** (-k * t) * powers[t] for t in range(r)) / r        # projector onto weight k (rows)
        if k == 0 or 2 * k == r:
            P = P.real                                                    # real (rational) weight space
            if Mss is None:
                rows = _greedy(list(P), start=invariants if k == 0 else ())
                half = (len(rows) + 1) // 2
                even, odd = rows[:half], rows[half:]
            else:
                even = _greedy(list(P @ (I + Mss) / 2), start=invariants if k == 0 else ())
                odd = _greedy(list(P @ (I - Mss) / 2))
            n = max(len(even), len(odd))
            even = even + [np.zeros(m)] * (n - len(even))
            odd = odd + [np.zeros(m)] * (n - len(odd))
            by_residue.setdefault(k, []).extend(a + 1j * b for a, b in zip(even, odd))
        else:
            rows = list(P)
            if Mss is None:
                basis = _greedy(rows, rank=_complex_rank)
            else:
                tau = [np.conj(u @ Mss) for u in rows]                     # tau(p) = conj(p M_s)
                candidates = [x for u, t in zip(rows, tau) for x in (u + t, 1j * (u - t))]
                basis = _greedy(candidates)
            for p in basis:                                               # balance residues k and r - k
                if len(by_residue.get(k, [])) <= len(by_residue.get(r - k, [])):
                    by_residue.setdefault(k, []).append(p)
                else:
                    by_residue.setdefault(r - k, []).append(np.conj(p))

    used = sorted(rho for rho, coords in by_residue.items() if coords)
    if not used:                                                      # everything is invariant: nothing to pack
        mode, rho, weights, Phi = "none", None, [], np.zeros((0, m), dtype=complex)
    elif len(used) == 1 and used[0] != 0:
        mode, rho = "scalar", used[0]
        coords = by_residue[rho]
        weights = [rho] * len(coords)
        Phi = np.array(coords)
    else:
        mode, rho = "x", None
        W = max(res + r * (len(coords) - 1) for res, coords in by_residue.items() if coords) + 1
        Phi = np.zeros((W, m), dtype=complex)
        weights = [None] * W
        for res, coords in by_residue.items():
            for t, p in enumerate(coords):
                Phi[res + r * t] = p
                weights[res + r * t] = res
    inv = np.array(invariants).reshape(len(invariants), m)
    return GradedPacking(h, s, r, mode, rho, drop_trivial, inv, Phi.real.copy(), Phi.imag.copy(), weights)


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
                 packing: GradedPacking):
        super().__init__(n, k, slots, index, sign, dtype)
        selection = output_selection(packing, len(index))
        self.register_buffer("inv", torch.tensor(packing.invariants, dtype=dtype).reshape(-1, self.m))
        self.register_buffer("graded_plus", torch.tensor(packing.plus, dtype=dtype))
        self.register_buffer("graded_minus", torch.tensor(packing.minus, dtype=dtype))
        self.register_buffer("select", torch.tensor([t[3] for t in selection], dtype=torch.int64))
        self.output_shape = torch.tensor([self.inv.shape[0] + len(selection)], dtype=torch.int64)

    def forward(self, x) -> torch.Tensor:
        x = self._prepare(x)
        batch_shape = x.shape[:-2]
        w = self.w(x).reshape(-1, self.m)
        if self.graded_plus.shape[0] == 0:
            out = w @ self.inv.T
        else:
            gw = w[:, self.g_index] * self.g_sign
            re, im = bursar_coefficients_complex(gw @ self.graded_plus.T, gw @ self.graded_minus.T)
            out = torch.cat([re, im], -1)[:, self.select]
            if self.inv.shape[0]:
                out = torch.cat([w @ self.inv.T, out], -1)
        return out.reshape(*batch_shape, out.shape[-1])


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

    def __init__(self, group_spec, k, options, slots, elements, index, sign, packing: GradedPacking, e_dropped: bool):
        self.group_spec, self.k, self.options, self.slots = group_spec, k, options, slots
        self.elements, self.index, self.sign, self.packing, self.e_dropped = elements, index, sign, packing, e_dropped
        self.relations = ()
        self.echinocoder_commit = _echinocoder_commit()

    labels = property(lambda self: self.group_spec.labels)
    n = property(lambda self: self.group_spec.n)
    m = property(lambda self: len(self.slots))
    N = property(lambda self: len(self.elements))
    certification_level = 0

    def slot_names(self) -> list[str]:
        return [s.name(self.labels) for s in self.slots]

    def selection(self) -> list[tuple]:
        return output_selection(self.packing, self.N)

    @property
    def D(self) -> int:
        return len(self.packing.invariants) + len(self.selection())

    def output_names(self) -> list[str]:
        names = [f"inv[{t}]" for t in range(len(self.packing.invariants))]
        names += [f"{'Re' if part == 're' else 'Im'} E_{{{c},{j}}}" for c, j, part, _ in self.selection()]
        return names

    @property
    def certification(self) -> str:
        return CERTIFICATION_GRADED + (" " + CERTIFICATION_E_DROPPED if self.e_dropped else "")

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
                        "reason": "not applied: the exact (integer) linear-dependency search does not yet support "
                                  "graded packings, whose coefficients involve exp(2 pi i / r)"})
        return log

    def reconstruct_base(self, Y):
        return Y

    def possible_degrees(self) -> list[list[int]]:
        inv = [sorted({self.slots[t].L_degree(self.k) for t in range(self.m) if abs(v[t]) > 1e-12})
               for v in self.packing.invariants]
        if self.packing.W == 0:
            return inv
        kinds = [{self.slots[t].kind for t in range(self.m) if abs(p[t]) > 1e-12 or abs(q[t]) > 1e-12}
                 for p, q in zip(self.packing.plus, self.packing.minus)]
        full = possible_L_degrees(self.slots, self.k, self.N, kinds)
        position = {lab: i for i, lab in enumerate(bursar_output_labels(self.N, self.packing.W))}
        return inv + [full[position[(c, j)]] for c, j, _, _ in self.selection()]

    def to_dict(self) -> dict:
        p = self.packing
        degrees = self.possible_degrees()
        return {
            "symembed_version": SYMEMBED_VERSION, "echinocoder_commit": self.echinocoder_commit,
            "group_spec": self.group_spec.to_dict(), "k": self.k, "options": self.options.to_dict(),
            "n": self.n, "m": self.m, "multiset_size": self.N, "D": self.D,
            "certification_level": 0, "certification": self.certification, "reductions": self.reductions_log(),
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
