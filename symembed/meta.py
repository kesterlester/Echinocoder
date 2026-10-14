"""
Metadata describing a built embedding: what each output is, symbolically, and on what its injectivity rests.

The metadata is for documentation and auditing only; evaluation never uses it. The symbolic description never
expands below the dots d and determinants e (its leaves), and uses nested definitions (q_g, B, E_c) rather than
expanding everything.
"""
from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .bursar import bursar_output_labels
from .reductions import LinearRelation
from .slots import GroupElement, Slot, slot_orbits
from .spec import EmbeddingOptions, GroupSpec

SYMEMBED_VERSION = "0.1.0"

CERTIFICATION_LEVEL_0_BURSAR = (
    "Level 0 (by construction): Bursarial embedding. The polynomial B(x, y) = prod_g (y + q_g(x)) "
    "determines the multiset {q_g} by unique factorisation in F[x][y], hence determines MS. MS determines L up to "
    "G and rotations (up to G and O(k) if parity is in G)."
)

CERTIFICATION_E_DROPPED = (
    "Determinants omitted because parity is in G (Level 0 for the smaller w): over the reals the Gram matrix "
    "determines L up to O(k), which is all that is required when parity is in G. Parity then acts trivially on w, "
    "so the multiset is taken over the permutation part of G only."
)

CERTIFICATION_INVOLUTION = (
    "Involution packing (Level 0): g0 acts on the packed coordinates as complex conjugation, and g -> g0 g permutes "
    "the group, so the multiset of packed polynomials is closed under conjugation and B has real coefficients. The "
    "complex Bursarial embedding is injective on multisets of complex vectors (unique factorisation in C[x][y]) and "
    "the packing is a bijection, so the real coefficients determine MS."
)

CERTIFICATION_LEVEL_1_LINEAR = (
    "Level 1 (explicit identity): each dropped output is an exact rational linear combination of kept outputs, as "
    "polynomials in w, verified in exact integer arithmetic. Injectivity is inherited from the Level 0 embedding."
)


_LATEX_ESCAPES = {"\\": "\\textbackslash{}", "_": "\\_", "{": "\\{", "}": "\\}", "^": "\\^{}",
                  "&": "\\&", "%": "\\%", "#": "\\#", "$": "\\$"}


def _latex_text(text: str) -> str:
    """Escape plain text for use in LaTeX text mode (single pass, so replacements are not re-escaped)."""
    return re.sub(r"[\\_{}^&%#$]", lambda match: _LATEX_ESCAPES[match.group()], text)


def _latex_form(text: str) -> str:
    """A packed-coordinate form like 'd_aa + d_bb' as LaTeX: d_aa -> d_{aa}, '*' removed."""
    return re.sub(r"([de])_([A-Za-z0-9,]+)", r"\1_{\2}", text).replace("*", "")


def _echinocoder_commit() -> str | None:
    try:
        result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=Path(__file__).parent,
                                capture_output=True, text=True, timeout=5)
        return result.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def possible_L_degrees(slots: tuple, k: int, N: int, position_kinds: list | None = None) -> list[list[int]]:
    """
    For each Bursarial output E_{c,j}, the L-degrees of the monomials that can appear in it (before any cancellation).

    A monomial of E_{c,j} is a product of c values, one from each of c distinct q_g, whose powers of x sum to j. Each
    d has L-degree 2 and each e has L-degree k, so the degree is 2c + (k-2)t where t counts the e factors. We track
    which (j, t) are reachable. By default position j of q_g holds a single slot (``slots[j]``); ``position_kinds``
    (a set of "d"/"e" per position) covers packed coordinates, which may combine slots of both kinds.
    """
    if position_kinds is None:
        position_kinds = [{s.kind} for s in slots]
    m = len(position_kinds)
    L = N * (m - 1) + 1
    reach = np.zeros((L, N + 1), dtype=bool)        # reach[j, t] for the current c
    for s, kinds in enumerate(position_kinds):
        for kind in kinds:
            reach[s, int(kind == "e")] = True
    result = []
    for c in range(1, N + 1):
        if c > 1:
            new = np.zeros_like(reach)
            for s, kinds in enumerate(position_kinds):
                for kind in kinds:
                    shift_t = int(kind == "e")
                    new[s:, shift_t:] |= reach[:L - s, :N + 1 - shift_t]
            reach = new
        for j in range(c * (m - 1) + 1):
            ts = np.nonzero(reach[j])[0]
            result.append(sorted({2 * c + (k - 2) * int(t) for t in ts}))
    return result


@dataclass
class EmbeddingMeta:
    """
    Everything needed to document, audit and rebuild an embedding.

    The embedding is the Bursarial embedding of {g.w} (the *base* outputs E_{c,j}), of which the subset ``kept`` is
    output, in that order. Dropped base outputs are recorded in ``relations`` as exact linear combinations of kept ones.

    Use ``to_dict``/``to_json`` for machine-readable output, ``latex`` for a human-readable description,
    ``expanded_output`` to see an output fully expanded in the d's and e's (small cases only), and
    ``reconstruct_base`` to recompute every base output from the kept ones.
    """
    group_spec: GroupSpec
    k: int
    options: EmbeddingOptions
    slots: tuple
    elements: tuple
    index: list
    sign: list
    kept: list | None = None                 # base indices of the outputs, in output order; None means all
    relations: tuple = ()                    # LinearRelation per dropped base output, indices are base indices
    e_dropped: bool = False
    linear_reduction_max_c: int = 0
    packing: object = None                   # involution.Packing when embedder == "bursar_involution"
    involution: int | None = None            # index in ``elements`` of the involution g0 used for the packing
    echinocoder_commit: str | None = field(default_factory=_echinocoder_commit)
    _expanded_base: list | None = field(default=None, repr=False)

    # ----- sizes and naming ------------------------------------------------------------------------------------

    @property
    def n(self) -> int:
        return self.group_spec.n

    @property
    def m(self) -> int:
        return len(self.slots)

    @property
    def N(self) -> int:
        """The number of members of the multiset (|G|, or |G|/2 if parity acts trivially on w)."""
        return len(self.elements)

    @property
    def width(self) -> int:
        """Number of coefficients of each q_g: m, or the number of packed coordinates."""
        return self.m if self.packing is None else self.packing.r

    def position_kinds(self) -> list:
        """For each position of q_g, the kinds ("d", "e") of the slots it involves."""
        if self.packing is None:
            return [{s.kind} for s in self.slots]
        return [{self.slots[t].kind for t in range(self.m) if p[t] or q[t]}
                for p, q in zip(self.packing.plus, self.packing.minus)]

    def packed_coordinates(self) -> list[tuple[str, str]]:
        """Each packed coordinate as (real part, imaginary part) in the slot names; empty if not packed."""
        return [] if self.packing is None else self.packing.describe(self.slot_names())

    @property
    def base_output_labels(self) -> list[tuple[int, int]]:
        """(c, j) for each base output: base output i is E_{c,j}."""
        return bursar_output_labels(self.N, self.width)

    @property
    def kept_base_indices(self) -> list[int]:
        return list(range(len(self.base_output_labels))) if self.kept is None else list(self.kept)

    @property
    def output_labels(self) -> list[tuple[int, int]]:
        """(c, j) for each output: output i is E_{c,j}."""
        base = self.base_output_labels
        return [base[i] for i in self.kept_base_indices]

    @property
    def D(self) -> int:
        return len(self.kept_base_indices)

    @property
    def labels(self) -> tuple:
        return self.group_spec.labels

    def slot_names(self) -> list[str]:
        return [s.name(self.labels) for s in self.slots]

    def base_name(self, i: int) -> str:
        c, j = self.base_output_labels[i]
        return f"E_{{{c},{j}}}"

    @property
    def certification(self) -> str:
        parts = []
        if self.relations:
            parts.append(CERTIFICATION_LEVEL_1_LINEAR)
        parts.append(CERTIFICATION_LEVEL_0_BURSAR)
        if self.packing is not None:
            parts.append(CERTIFICATION_INVOLUTION)
        if self.e_dropped:
            parts.append(CERTIFICATION_E_DROPPED)
        return " ".join(parts)

    @property
    def certification_level(self) -> int:
        return 1 if self.relations else 0

    # ----- reductions ------------------------------------------------------------------------------------------

    def relation_summary(self, relation: LinearRelation) -> str:
        """Plain text, e.g. ``E_{1,1} = E_{1,0}``."""
        lhs = self.base_name(relation.dropped)
        if not relation.terms:
            return f"{lhs} = 0"
        pieces = []
        for t, frac in relation.terms:
            coefficient = "" if frac == 1 else ("-" if frac == -1 else f"{frac} ")
            pieces.append(f"{coefficient}{self.base_name(t)}")
        return f"{lhs} = " + " + ".join(pieces).replace("+ -", "- ")

    def reductions_log(self) -> list[dict]:
        log = []
        if self.e_dropped:
            log.append({"type": "omit_determinants_with_parity", "level": 0, "reason": CERTIFICATION_E_DROPPED})
        if self.options.reductions == "safe":
            log.append({"type": "linear_scope", "level": 1,
                        "reason": f"exact linear dependencies checked for c <= {self.linear_reduction_max_c}; "
                                  f"outputs with larger c are all kept"})
        for rel in self.relations:
            log.append({"type": "linear", "level": 1, "dropped": self.base_name(rel.dropped),
                        "dropped_base_index": rel.dropped,
                        "terms": [{"base_index": t, "name": self.base_name(t), "coefficient": str(frac)}
                                  for t, frac in rel.terms],
                        "summary": self.relation_summary(rel),
                        "verified": "exact (integer arithmetic)"})
        return log

    def reconstruct_base(self, Y):
        """
        All base outputs E_{c,j}, recomputed from the kept outputs ``Y`` (shape (..., D)) using the recorded
        relations. Works for torch tensors and NumPy arrays.
        """
        base_D = len(self.base_output_labels)
        columns = [None] * base_D
        for out, b in enumerate(self.kept_base_indices):
            columns[b] = Y[..., out]
        for rel in sorted(self.relations, key=lambda r: r.dropped):
            value = 0 * Y[..., 0]
            for t, frac in rel.terms:
                value = value + columns[t] * (frac.numerator / frac.denominator)
            columns[rel.dropped] = value
        stack = __import__("torch").stack if hasattr(Y, "detach") else np.stack
        return stack(columns, -1)

    # ----- machine readable ------------------------------------------------------------------------------------

    def to_dict(self) -> dict:
        degrees = possible_L_degrees(self.slots, self.k, self.N, self.position_kinds())
        return {
            "symembed_version": SYMEMBED_VERSION,
            "echinocoder_commit": self.echinocoder_commit,
            "group_spec": self.group_spec.to_dict(),
            "k": self.k,
            "options": self.options.to_dict(),
            "n": self.n,
            "m": self.m,
            "multiset_size": self.N,
            "D": self.D,
            "base_D": len(self.base_output_labels),
            "certification_level": self.certification_level,
            "certification": self.certification,
            "reductions": self.reductions_log(),
            "slots": [{"name": s.name(self.labels), "kind": s.kind,
                       "labels": [self.labels[i] for i in s.positions], "L_degree": s.L_degree(self.k)}
                      for s in self.slots],
            "slot_orbits": [list(run) for run in slot_orbits(self.group_spec, self.slots)],
            "group_elements": [{"cycles": g.cycle_notation(self.labels),
                                "perm": [self.labels[p] for p in g.perm], "parity": g.parity}
                               for g in self.elements],
            "signed_permutation": {"index": self.index, "sign": self.sign,
                                   "meaning": "(g.w)[s] = sign[g][s] * w[index[g][s]]"},
            "packing": None if self.packing is None else {
                "involution": self.elements[self.involution].cycle_notation(self.labels),
                "trace": sum(self.sign[self.involution][s] for s in range(self.m)
                             if self.index[self.involution][s] == s),
                "coordinates": [{"re": re_, "im": im_} for re_, im_ in self.packed_coordinates()],
                "meaning": "z_j = re_j(g.w) + i im_j(g.w); q_g(x) = sum_j z_j x^j; outputs are real parts",
            },
            "definitions": {
                "q_g(x)": ("sum_{j=0}^{m-1} (g.w)_j x^j" if self.packing is None
                           else "sum_j z_j(g.w) x^j, with z_j the packed coordinates"),
                "B(x,y)": "prod_{g} (y + q_g(x)) = sum_{c=0}^{N} y^(N-c) E_c(x), product over the group elements",
                "E_c(x)": "sum_j E_{c,j} x^j",
            },
            "outputs": [{"index": i, "name": self.base_name(b), "base_index": b, "c": c, "j": j, "w_degree": c,
                         "possible_L_degrees": degrees[b]}
                        for i, (b, (c, j)) in enumerate(zip(self.kept_base_indices, self.output_labels))],
        }

    def to_json(self, **kwargs) -> str:
        return json.dumps(self.to_dict(), **kwargs)

    # ----- symbolic expansion (small cases only) ---------------------------------------------------------------

    def slot_symbols(self):
        import sympy
        return [sympy.Symbol(s.name(self.labels)) for s in self.slots]

    def expanded_base_outputs(self, max_group_order: int = 24) -> list:
        """
        Every base output fully expanded as a polynomial in the d's and e's (sympy expressions). Feasible only for
        small groups; raises ValueError if the multiset has more than max_group_order members.
        """
        if self._expanded_base is not None:
            return self._expanded_base
        if self.N > max_group_order:
            raise ValueError(f"multiset size {self.N} exceeds max_group_order = {max_group_order}; full expansion "
                             f"is likely infeasible. Raise max_group_order to force it.")
        import sympy
        x = sympy.Symbol("x")
        w = self.slot_symbols()
        domain = sympy.ZZ.poly_ring(*w)

        def poly(expr):
            return sympy.Poly(expr, x, domain=domain)

        # Complex polynomials as (real part, imaginary part).
        E = [(poly(1), poly(0))]
        for g in range(self.N):
            gw = [self.sign[g][s] * w[self.index[g][s]] for s in range(self.m)]
            if self.packing is None:
                q = (poly(sum(gw[j] * x ** j for j in range(self.m))), poly(0))
            else:
                q = (poly(sum(sum(c * v for c, v in zip(row, gw)) * x ** j
                              for j, row in enumerate(self.packing.plus))),
                     poly(sum(sum(c * v for c, v in zip(row, gw)) * x ** j
                              for j, row in enumerate(self.packing.minus))))

            def times(a, b):
                return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

            E = ([E[0]] + [(E[c][0] + times(q, E[c - 1])[0], E[c][1] + times(q, E[c - 1])[1])
                           for c in range(1, len(E))] + [times(q, E[-1])])
        for c, j in self.base_output_labels:
            if E[c][1].coeff_monomial(x ** j) != 0:
                raise AssertionError(f"E_{{{c},{j}}} has a non-zero imaginary part")
        self._expanded_base = [E[c][0].coeff_monomial(x ** j) for c, j in self.base_output_labels]
        return self._expanded_base

    def expanded_outputs(self, max_group_order: int = 24) -> list:
        """Every (kept) output fully expanded in the d's and e's. Small cases only."""
        base = self.expanded_base_outputs(max_group_order)
        return [base[b] for b in self.kept_base_indices]

    def expanded_output(self, i: int, max_group_order: int = 24):
        return self.expanded_outputs(max_group_order)[i]

    # ----- LaTeX -----------------------------------------------------------------------------------------------

    def latex(self, max_rows: int = 200, expanded: bool = False) -> str:
        """
        A LaTeX fragment (needs amsmath and longtable) describing the embedding: slots of w, the q_g, the
        definitions, a table of outputs and a table of reductions. Tables longer than ``max_rows`` are truncated,
        with a note. With ``expanded=True`` (small cases only) each output is also shown expanded in the d's and e's.
        """
        import sympy
        L = self.labels
        parts = []
        spec = " \\times ".join(f"S({len(s)})" for s in self.group_spec.species)
        if self.group_spec.parity:
            spec += " \\times P"
        parts.append(f"\\noindent Embedding of $n = {self.n}$ vectors in $k = {self.k}$ dimensions, "
                     f"$G = {spec}$, multiset size $N = {self.N}$, $m = {self.m}$, $D = {self.D}$ outputs "
                     f"(of {len(self.base_output_labels)} Bursarial coefficients). "
                     f"Certification: {_latex_text(self.certification)}\n")

        orbits = slot_orbits(self.group_spec, self.slots)
        orbit_end = {run[-1] for run in orbits}
        rows = []
        for j, s in enumerate(self.slots):
            rows.append(f"{j} & ${s.latex(L)}$ & {s.L_degree(self.k)} \\\\" + (" \\hline" if j in orbit_end else ""))
        parts.append("\\begin{longtable}{rcc}\nslot $j$ & $w_j$ & degree in $L$ \\\\ \\hline\n"
                     + "\n".join(rows) + "\n\\end{longtable}\n")

        parts.append("\\[ q_g(x) = \\sum_{j=0}^{m-1} (g\\cdot w)_j\\, x^j, \\qquad "
                     "B(x,y) = \\prod_{g} \\bigl(y + q_g(x)\\bigr) = \\sum_{c=0}^{N} y^{N-c} E_c(x), "
                     "\\qquad E_c(x) = \\sum_j E_{c,j}\\, x^j. \\]\n")

        if self.packing is not None:
            rows = [f"{j} & ${_latex_form(re_)}$ & ${_latex_form(im_)}$ \\\\"
                    for j, (re_, im_) in enumerate(self.packed_coordinates())]
            g0 = self.elements[self.involution].cycle_notation(L)
            parts.append(f"Involution packing with $g_0 = {g0}$: each $q_g$ has packed coefficients "
                         "$z_j(g\\cdot w) = \\mathrm{Re}_j + i\\,\\mathrm{Im}_j$, on which $g_0$ acts as complex "
                         "conjugation, so the outputs (real parts) are the whole of $B$.\n")
            parts.append("\\begin{longtable}{rll}\n$j$ & real part & imaginary part \\\\ \\hline\n"
                         + "\n".join(rows) + "\n\\end{longtable}\n")

        def signed(g, j):
            s = self.slots[self.index[g][j]].latex(L)
            return ("-" if self.sign[g][j] < 0 else "") + s

        rows = []
        for g in range(min(self.N, max_rows)):
            cells = ", ".join(signed(g, j) for j in range(self.m))
            rows.append(f"${self.elements[g].cycle_notation(L)}$ & $({cells})$ \\\\")
        if self.N > max_rows:
            rows.append(f"\\multicolumn{{2}}{{l}}{{\\dots\\ {self.N - max_rows} more group elements}} \\\\")
        parts.append("\\begin{longtable}{ll}\n$g$ & $g\\cdot w$ (coefficients of $q_g$) \\\\ \\hline\n"
                     + "\n".join(rows) + "\n\\end{longtable}\n")

        degrees = possible_L_degrees(self.slots, self.k, self.N, self.position_kinds())
        expanded_exprs = self.expanded_outputs() if expanded else None
        rows = []
        for i, (b, (c, j)) in enumerate(list(zip(self.kept_base_indices, self.output_labels))[:max_rows]):
            deg = ", ".join(str(d) for d in degrees[b]) if degrees[b] else "--"
            row = f"{i} & $E_{{{c},{j}}}$ & {c} & {deg}"
            if expanded:
                row += f" & ${sympy.latex(expanded_exprs[i])}$"
            rows.append(row + " \\\\")
        if self.D > max_rows:
            rows.append(f"\\multicolumn{{4}}{{l}}{{\\dots\\ {self.D - max_rows} more outputs}} \\\\")
        header = "output & & degree in $w$ & possible degrees in $L$"
        colspec = "rccc"
        if expanded:
            header += " & expanded"
            colspec += "l"
        parts.append(f"\\begin{{longtable}}{{{colspec}}}\n{header} \\\\ \\hline\n"
                     + "\n".join(rows) + "\n\\end{longtable}\n")

        if self.relations:
            rows = []
            for rel in self.relations[:max_rows]:
                rhs = self.relation_summary(rel).split(" = ", 1)[1]
                rhs = re.sub(r"(\d+)/(\d+) ", r"\\tfrac{\1}{\2}\\,", rhs)
                rows.append(f"${self.base_name(rel.dropped)}$ & $= {rhs}$ \\\\")
            if len(self.relations) > max_rows:
                rows.append(f"\\multicolumn{{2}}{{l}}{{\\dots\\ {len(self.relations) - max_rows} more}} \\\\")
            parts.append("\\begin{longtable}{ll}\ndropped & exact identity (Level 1) \\\\ \\hline\n"
                         + "\n".join(rows) + "\n\\end{longtable}\n")
        return "\n".join(parts)
