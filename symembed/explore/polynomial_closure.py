"""
Exploration: which Bursarial outputs are exact polynomials in earlier ones?

For each case, the outputs E_{c,j} are visited in order (c, then j). An output is kept unless it equals an exact
polynomial (with rational coefficients) in the outputs already kept. Membership is tested in one of two settings:

* ``"w"``: formally, as polynomials in the slots of w (the d's and e's treated as independent variables);
* ``"L"``: as polynomials in the components of the vectors, so that the Gram and Plucker relations among the d's
  and e's hold automatically.

Every identity found is verified exactly (see ``symembed.reductions.find_linear_relations``), so the kept set is a
Level 1 certified separating set: every Bursarial output is an explicit polynomial in it.

Run from the repository root:

    /Users/lester/github/Echinocoder/venv/bin/python -m symembed.explore.polynomial_closure
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path

import numpy as np
import sympy
from sympy.polys.domains import ZZ
from sympy.polys.rings import ring

from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding
from symembed.explore.simple_cases import CASES, OUT_DIR, quotient_dimension


def slot_values(meta, setting: str):
    """The ring, and each slot of w as an element of it."""
    if setting == "w":
        names = meta.slot_names()
        R, *gens = ring(",".join(names), ZZ)
        return R, list(gens)
    k, labels = meta.k, meta.labels
    names = [f"{label}{axis}" for label in labels for axis in range(k)]
    R, *gens = ring(",".join(names), ZZ)
    vec = [gens[i * k:(i + 1) * k] for i in range(len(labels))]
    values = []
    for s in meta.slots:
        if s.kind == "d":
            i, j = s.positions
            values.append(sum((vec[i][t] * vec[j][t] for t in range(k)), R.zero))
        else:
            M = sympy.Matrix([[sympy.Symbol(str(vec[i][t])) for t in range(k)] for i in s.positions])
            values.append(R(sympy.expand(M.det())))
    return R, values


def bursar_polys(meta, setting: str) -> dict:
    """Every Bursarial output E_{c,j} as an exact polynomial in the chosen setting."""
    R, w = slot_values(meta, setting)
    m, N = meta.m, meta.N
    E = [[R.one]]                                   # E[c][j]
    for g in range(N):
        q = [meta.sign[g][jp] * w[meta.index[g][jp]] for jp in range(m)]
        new = [E[0]] + [list(row) for row in E[1:]] + [[R.zero] * (len(E) * (m - 1) + 1)]
        for c in range(1, len(E) + 1):
            prev = E[c - 1]
            for jp, qv in enumerate(q):
                if not qv:
                    continue
                for j, pv in enumerate(prev):
                    if pv:
                        new[c][j + jp] += qv * pv
        E = new
    return {(c, j): E[c][j] for c in range(1, N + 1) for j in range(c * (m - 1) + 1)}


@dataclass
class ClosureResult:
    setting: str
    kept: list                 # (c, j) of kept outputs, in order
    identities: dict           # (c, j) -> list of (Fraction, tuple of (c, j) factors)
    unverified: list           # (c, j) of outputs kept only because a candidate identity failed exact verification


PRIMES = (2147483629, 2147483587, 2147483579, 2147483563)     # products of residues fit in int64


def _rational_reconstruct_big(a: int, modulus: int) -> Fraction | None:
    """The fraction r/s with |r|, |s| <= sqrt(modulus/2) congruent to a mod modulus, if one exists."""
    a %= modulus
    bound = isqrt(modulus // 2)
    r0, r1, s0, s1 = modulus, a, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound or gcd(r1, abs(s1)) != 1:
        return None
    return Fraction(r1, s1)


def _crt(residues, moduli) -> tuple[int, int]:
    x, M = 0, 1
    for r, m in zip(residues, moduli):
        # x = x (mod M), x = r (mod m)
        t = ((r - x) * pow(M, -1, m)) % m
        x, M = x + M * t, M * m
    return x, M


class _IncrementalSpan:
    """Row echelon form modulo a prime, tracking each basis row as a combination of the rows added."""

    def __init__(self, ncols: int, nrows_max: int, prime: int):
        self.ncols, self.nrows_max, self.prime = ncols, nrows_max, prime
        self.pivots: list[int] = []
        self.rows: list[np.ndarray] = []
        self.combs: list[np.ndarray] = []

    def _reduce(self, row: np.ndarray):
        P = self.prime
        r = row % P
        comb = np.zeros(self.nrows_max, dtype=np.int64)       # r = row - sum comb[i] * added_i
        for pivot, b, bc in zip(self.pivots, self.rows, self.combs):
            f = int(r[pivot])
            if f:
                r = (r - f * b) % P
                comb = (comb + f * bc) % P
        return r, comb

    def add(self, index: int, row: np.ndarray) -> bool:
        r, comb = self._reduce(row)
        nz = np.nonzero(r)[0]
        if nz.size == 0:
            return False
        P = self.prime
        pivot = int(nz[0])
        inv = pow(int(r[pivot]), P - 2, P)
        own = np.zeros(self.nrows_max, dtype=np.int64)
        own[index] = 1
        self.pivots.append(pivot)
        self.rows.append((r * inv) % P)
        self.combs.append(((own - comb) * inv) % P)
        return True

    def express(self, row: np.ndarray):
        """Coefficients (mod the prime) expressing row in terms of rows added, or None if independent."""
        r, comb = self._reduce(row)
        return None if np.any(r) else comb


def _exact_check(f, terms, polys_of_terms) -> bool:
    denominator = 1
    for frac in terms:
        denominator = denominator * frac.denominator // gcd(denominator, frac.denominator)
    total = f * denominator
    for frac, p in zip(terms, polys_of_terms):
        total -= p * (frac.numerator * (denominator // frac.denominator))
    return not total


def polynomial_closure(meta, setting: str) -> ClosureResult:
    """
    Visit outputs in order of c; keep an output unless it lies in the algebra generated by the kept outputs.

    For each c, a basis of the span of all products of kept outputs of total w-degree c (or at most c, when the
    outputs are not homogeneous as functions) is grown incrementally, as products (generator) x (basis element of
    lower degree). Identities are found modulo a prime and verified exactly before an output is dropped.
    """
    polys = bursar_polys(meta, setting)
    graded = setting == "w" or meta.k == 2       # outputs are then homogeneous, of degree c (in w or in L)
    N = meta.N
    kept, identities, unverified = [], {}, []
    generators: dict[int, list] = {}             # c -> list of (combo, poly) for kept outputs of degree c
    basis: dict[int, list] = {0: []}             # c -> list of (combo, poly), independent
    for c in range(1, N + 1):
        outputs = [((cc, j), f) for (cc, j), f in polys.items() if cc == c]
        candidates = list(basis[c - 1]) if not graded else []
        seen = {combo for combo, _ in candidates}
        for cg in range(1, c):
            lower = basis[c - cg] if graded else basis[c - cg]
            for g_combo, g in generators.get(cg, []):
                for a_combo, a in lower:
                    combo = tuple(sorted(g_combo + a_combo))
                    if combo not in seen:
                        seen.add(combo)
                        candidates.append((combo, g * a))
        rows = [p for _, p in candidates] + [f for _, f in outputs]
        monomials = sorted({mono for p in rows for mono in p.keys()})
        col = {mono: t for t, mono in enumerate(monomials)}

        def dense(p, prime):
            v = np.zeros(len(monomials), dtype=np.int64)
            for mono, coeff in p.items():
                v[col[mono]] = int(coeff) % prime
            return v

        # One echelon form per prime. A candidate is a member if independent modulo the first prime; if another
        # prime disagrees (an unlucky prime), identities through it fail exact verification and are not used.
        spans = [_IncrementalSpan(len(monomials), len(rows), prime) for prime in PRIMES]
        members: list = []                        # index in addition order -> (combo, poly), or None
        for combo, p in candidates:
            independent = [span.add(len(members), dense(p, span.prime)) for span in spans]
            members.append((combo, p) if independent[0] else None)
        basis[c] = [m for m in members if m is not None]
        for (cc, j), f in outputs:
            if not f:
                identities[(cc, j)] = []
                continue
            combs = [span.express(dense(f, span.prime)) for span in spans]
            if all(comb is not None for comb in combs):
                support = sorted(set().union(*(set(np.nonzero(comb)[0].tolist()) for comb in combs)))
                terms, used, ok = [], [], True
                for t in support:
                    value, modulus = _crt([int(comb[t]) for comb in combs], PRIMES)
                    frac = _rational_reconstruct_big(value, modulus)
                    if frac is None or members[t] is None:
                        ok = False
                        break
                    if frac:
                        terms.append(frac)
                        used.append(members[t])
                if ok and _exact_check(f, terms, [p for _, p in used]):
                    identities[(cc, j)] = [(frac, combo) for frac, (combo, _) in zip(terms, used)]
                    continue
                unverified.append((cc, j))
            kept.append((cc, j))
            generators.setdefault(c, []).append((((cc, j),), f))
            for span in spans:
                span.add(len(members), dense(f, span.prime))
            members.append((((cc, j),), f))
            basis[c].append((((cc, j),), f))
    return ClosureResult(setting, kept, identities, unverified)


def _degree(p) -> int:
    """Total degree of a ring element (in w or in the vector components)."""
    return max((sum(mono) for mono in p.keys()), default=0)


def _express_exactly(f, candidates):
    """Exact rational coefficients expressing f in the span of the candidate polys, or None."""
    rows = [p for _, p in candidates] + [f]
    monomials = sorted({mono for p in rows for mono in p.keys()})
    col = {mono: t for t, mono in enumerate(monomials)}

    def dense(p, prime):
        v = np.zeros(len(monomials), dtype=np.int64)
        for mono, coeff in p.items():
            v[col[mono]] = int(coeff) % prime
        return v

    spans = [_IncrementalSpan(len(monomials), len(rows), prime) for prime in PRIMES]
    members = []
    for combo, p in candidates:
        independent = [span.add(len(members), dense(p, span.prime)) for span in spans]
        members.append((combo, p) if independent[0] else None)
    combs = [span.express(dense(f, span.prime)) for span in spans]
    if any(comb is None for comb in combs):
        return None
    terms, used = [], []
    for t in sorted(set().union(*(set(np.nonzero(comb)[0].tolist()) for comb in combs))):
        value, modulus = _crt([int(comb[t]) for comb in combs], PRIMES)
        frac = _rational_reconstruct_big(value, modulus)
        if frac is None or members[t] is None:
            return None
        if frac:
            terms.append(frac)
            used.append(members[t])
    if not _exact_check(f, terms, [p for _, p in used]):
        return None
    return [(frac, combo) for frac, (combo, _) in zip(terms, used)]


def prune_backwards(meta, result: ClosureResult) -> ClosureResult:
    """
    Retry each kept output, latest first: drop it if it is an exact polynomial in the other kept outputs. Products
    are matched by total degree (in w, or in the vector components for setting "L"), so for example a^2 + b^2 can
    be expressed through (a + b)^2 even though e has lower degree in w than d.
    """
    polys = bursar_polys(meta, result.setting)
    kept = list(result.kept)
    identities = dict(result.identities)
    for target in reversed(list(kept)):
        others = [o for o in kept if o != target]
        f = polys[target]
        bound = _degree(f)
        degrees = {o: _degree(polys[o]) for o in others}
        candidates = []

        def extend(start, combo, deg, product):
            for i in range(start, len(others)):
                o = others[i]
                d = deg + degrees[o]
                if d <= bound:
                    new_combo, new_product = combo + (o,), (polys[o] if product is None else product * polys[o])
                    candidates.append((new_combo, new_product))
                    extend(i, new_combo, d, new_product)

        extend(0, (), 0, None)
        terms = _express_exactly(f, candidates) if candidates else None
        if terms is not None:
            kept.remove(target)
            identities[target] = terms
    return ClosureResult(result.setting, kept, identities, result.unverified)


def describe(term_list) -> str:
    if not term_list:
        return "0"
    pieces = []
    for frac, combo in term_list:
        factors = " ".join(f"E{c},{j}" for c, j in combo)
        coefficient = "" if frac == 1 else ("-" if frac == -1 else f"({frac}) ")
        pieces.append(f"{coefficient}{factors}")
    return " + ".join(pieces).replace("+ -", "- ")


def main() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "polynomial_closure.txt"
    lines = ["symembed: polynomial closure of the Bursarial outputs (generated by",
             "symembed/explore/polynomial_closure.py; do not edit by hand).",
             "",
             "An output is KEPT unless it is an exact polynomial in outputs already kept (visited in order c, j).",
             "Setting 'w': formal polynomials in the d's and e's. Setting 'L': polynomials in vector components,",
             "where the Gram/Plucker relations hold. Every identity was verified exactly. 'Ec,j' means E_{c,j}."]
    for title, k, spec in CASES:
        print(f"closure: {title}", file=sys.stderr)
        _, meta = build_invariant_embedding(k, spec, EmbeddingOptions())
        plain_title = title.replace("$", "").replace("\\times", "x")
        lines += ["", "=" * 100,
                  f"{plain_title}   (n = {spec.n}, k = {k}, |G| = {meta.N}, m = {meta.m}, "
                  f"Bursarial outputs = {len(meta.base_output_labels)}, quotient dimension = "
                  f"{quotient_dimension(spec.n, k)})",
                  "slots: " + ", ".join(meta.slot_names()), "=" * 100]
        expanded = dict(zip(meta.base_output_labels, meta.expanded_base_outputs()))
        for setting in ("w", "L"):
            forward = polynomial_closure(meta, setting)
            result = prune_backwards(meta, forward)
            by_c = {}
            for c, _ in result.kept:
                by_c[c] = by_c.get(c, 0) + 1
            lines += ["", f"--- setting {setting}: {len(result.kept)} kept after backward pruning "
                          f"({len(forward.kept)} after the forward pass); kept per c: {by_c}"]
            for c, j in result.kept:
                lines.append(f"  KEPT  E{c},{j} = {sympy.sstr(sympy.expand(expanded[(c, j)]))}")
            if setting == "L":
                for (c, j), terms in result.identities.items():
                    text = describe(terms)
                    if len(text) > 300:
                        text = text[:300] + f" ... ({len(terms)} terms)"
                    lines.append(f"  drop  E{c},{j} = {text}")
    out.write_text("\n".join(lines) + "\n")
    return out


if __name__ == "__main__":
    print(main())
