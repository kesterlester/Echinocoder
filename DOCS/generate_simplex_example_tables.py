#!/usr/bin/env python
"""Generate the LaTeX tables for the Algorithm-2 worked example in
simplicial_complex_embedders_v3.tex, by running EncDec.py on the example
multiset.  This keeps the document in lock-step with the code: to change the
worked example, edit M below and re-run (the Makefile does this automatically).

Fragments are written to DOCS/generated/ and \\input by the .tex source:

  sx2_input.tex            the multiset M displayed as column vectors
  sx2_minima.tex           inline list of row minima  m_0, m_1, m_2
  sx2_balance.tex          the balance multiset display
  sx2_stepC.tex            three mini Abel-summation tables (one per row)
  sx2_merged_deltas.tex    the merged gap list, sorted descending
  sx2_stepD_main.tex       the big 9-row table (default tie-break order)
  sx2_stepD_variant.tex    the same table under a rotated tie-break order
  sx2_sum_check.tex        the sum-equality display for the verification box
  sx2_threeway.tex         the three-way lin-comb equality for the same box
  sx2_count_matrices.tex   the Eji_LinComb count-matrix table
  sx2_stepF.tex            the final assembled-output display

Greying convention (matches the hand-drawn v2 tables): a cell is greyed iff
its value depends on the tie-break order, determined by comparing the default
and variant orderings cell-by-cell:
  - the (delta, V) pair is greyed iff V differs between the two orderings;
  - the (Delta, L) and (alpha, Lhat) pairs are greyed iff L differs.
Both tables receive the same grey pattern.
"""

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import EncDec  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent / "generated"

# ── The example multiset ──────────────────────────────────────────────────────
# Rows are vectors j = 0..n-1; columns are coordinates i = 0..k-1.
# This is the same example as in v1/v2 of the document.
M = np.asarray([[4, 2, 3],
                [-3, 5, 1],
                [8, 9, 2],
                [2, 7, 2]])


# ── Eji-sum rendering ─────────────────────────────────────────────────────────

def eji_term(coeff: int, j: int, i: int) -> str:
    """One term c*e^j_i, omitting a coefficient of 1."""
    c = "" if coeff == 1 else str(coeff)
    return rf"{c}\eji{{{j}}}{{{i}}}"


def eji_lines(basis_vec: np.ndarray) -> list[str]:
    """Render a basis matrix as one line of Eji terms per coordinate row i.

    basis_vec has shape (n, k) indexed [j][i].  Terms are grouped by i
    (matching the document's convention), ascending j within each group.
    Each returned line starts with its own leading '+' except the first.
    """
    bv = np.asarray(basis_vec).astype(int)
    n, k = bv.shape
    lines = []
    for i in range(k):
        terms = [eji_term(int(bv[j][i]), j, i) for j in range(n) if bv[j][i] != 0]
        if terms:
            lines.append("+".join(terms))
    return lines


def eji_sum_cell(basis_vec, prefactor_index: int | None = None) -> str:
    """Render a basis matrix as a (possibly multi-line) maths cell.

    If prefactor_index is given and > 1, wrap as \\tfrac{1}{idx}( ... ).
    Continuation lines are separated by \\newline in the p{} column style
    used throughout the document.
    """
    lines = eji_lines(basis_vec)
    if not lines:
        return "$0$"
    pre = ""
    close = ""
    if prefactor_index is not None and prefactor_index > 1:
        pre = rf"\tfrac{{1}}{{{prefactor_index}}} ("
        close = ")"
    if len(lines) == 1:
        return rf"${pre}{lines[0]}{close}$"
    parts = [rf"${pre}{lines[0]}$"]
    for line in lines[1:-1]:
        parts.append(rf"${{}}+{line}$")
    parts.append(rf"${{}}+{lines[-1]}{close}$")
    return "\\newline\n    ".join(parts)


def grey(cell: str) -> str:
    r"""Wrap a cell's content in \textcolor{gray}{...} (maths-safe)."""
    return rf"\textcolor{{gray}}{{{cell}}}"


# ── Abel summation bookkeeping (wraps EncDec) ────────────────────────────────

def per_row_first_diffs(data: np.ndarray):
    """First Abel summation, one call per coordinate row (as in EncDec's
    simplex_2_preprocess_steps Step 2).  Returns (diffs_per_row, offsets_per_row)
    where each diffs entry is a LinComb of (gap, prefix) pairs."""
    n, k = data.shape
    lin_comb_0 = [EncDec.array_to_lin_comb(data, fixed_axes={1: i}) for i in range(k)]
    diffs, offsets = [], []
    for i in range(k):
        d, o = EncDec.barycentric_subdivide(
            lin_comb_0[i], return_offset_separately=True, preserve_scale=False)
        diffs.append(d)
        offsets.append(o)
    return diffs, offsets


def merged_sorted_pairs(diffs_per_row, rotate_tie_groups: bool = False):
    """Merge per-row (gap, prefix) pairs and sort by gap descending (stable,
    preserving row order within ties — EncDec's behaviour).

    If rotate_tie_groups, each group of equal gaps is rotated right by one,
    producing the document's alternative tie-break ordering."""
    merged = []
    for d in diffs_per_row:
        merged.extend(zip(d.coeffs, d.basis_vecs))
    merged.sort(key=lambda x: -x[0])
    if rotate_tie_groups:
        out, s = [], 0
        while s < len(merged):
            e = s
            while e + 1 < len(merged) and merged[e + 1][0] == merged[s][0]:
                e += 1
            group = merged[s:e + 1]
            out.extend([group[-1]] + group[:-1])  # rotate right by one
            s = e + 1
        merged = out
    return merged


def second_abel_rows(pairs):
    """Compute the big-table rows from the sorted (delta, V) pairs.

    Returns a list of dicts with keys s, delta, V, Delta, L, alpha
    (Lhat is L with prefactor 1/(s+1), rendered at table time).
    Cross-checks the cumulative sums against EncDec's own second
    barycentric_subdivide output."""
    deltas = [p[0] for p in pairs]
    vs = [p[1] for p in pairs]
    m = len(pairs)
    rows = []
    L = np.zeros_like(np.asarray(vs[0]))
    for s in range(m):
        L = L + np.asarray(vs[s])
        nxt = deltas[s + 1] if s + 1 < m else 0
        Delta = deltas[s] - nxt
        rows.append(dict(s=s, delta=int(deltas[s]), V=np.asarray(vs[s]).astype(int),
                         Delta=int(Delta), L=L.astype(int), alpha=int((s + 1) * Delta)))

    # Cross-check against EncDec (preserve_scale=False gives (Delta, L) directly).
    lc = EncDec.LinComb()
    for c, b in pairs:
        lc += EncDec.MonoLinComb(c, b)
    enc = EncDec.barycentric_subdivide(lc, return_offset_separately=False,
                                       preserve_scale=False)
    assert len(enc) == m
    for s in range(m):
        assert enc.coeffs[s] == rows[s]["Delta"], f"Delta mismatch at s={s}"
        assert np.array_equal(np.asarray(enc.basis_vecs[s]).astype(int), rows[s]["L"]), \
            f"L mismatch at s={s}"
    return rows


# ── Table emitters ────────────────────────────────────────────────────────────

TABLE_PREAMBLE = r"""\begin{center}
\small
\setlength{\tabcolsep}{4pt}
\renewcommand{\arraystretch}{1.35}
"""


def emit_stepC(data: np.ndarray) -> str:
    """Three mini Abel-summation tables, one per coordinate row, each shaped
    like a small version of the Step-D table: input pair (x, single Eji),
    output pair (gap delta, cumulative prefix V), then the offset row."""
    n, k = data.shape
    out = []
    for i in range(k):
        vals = [(int(data[j][i]), j) for j in range(n)]
        vals.sort(key=lambda x: -x[0])  # descending; stable in j for ties
        body = []
        prefix = np.zeros((n, k), dtype=int)
        for r in range(n - 1):
            x, j = vals[r]
            prefix[j][i] = 1
            delta = x - vals[r + 1][0]
            body.append(
                rf"{r} & ${x}$ & $\eji{{{j}}}{{{i}}}$ & {delta}"
                f"\n  & {eji_sum_cell(prefix)}"
                "\n  \\\\[2pt]"
            )
        m_i, j_last = vals[-1]
        all_ones = "+".join(rf"\eji{{{j}}}{{{i}}}" for j in range(n))
        offset_row = (
            rf"\multicolumn{{5}}{{l}}{{offset $\;=\;x^{{({i})}}_{{({n-1})}}\cdot"
            rf"\textstyle\sum_j\eji{{j}}{{{i}}}\;=\;{m_i}\,({all_ones})$}} \\"
        )
        table = (
            TABLE_PREAMBLE
            + r"\begin{tabular}{c|cl|cl}" + "\n"
            + r"\toprule" + "\n"
            + rf"\multicolumn{{5}}{{l}}{{\textbf{{Row $i={i}$}}}} \\" + "\n"
            + r"\midrule" + "\n"
            + rf"  & & & $x^{{({i})}}_{{(r)}}-x^{{({i})}}_{{(r+1)}}$ & $V^{{({i})}}_{{r-1}}+\eji{{\sigma_{i}(r)}}{{{i}}}$ \\" + "\n"
            + r"\midrule" + "\n"
            + rf"$r$ & $x^{{({i})}}_{{(r)}}$ & $\eji{{\sigma_{i}(r)}}{{{i}}}$ & $\delta^{{({i})}}_r$ & $V^{{({i})}}_r$ \\" + "\n"
            + r"\midrule" + "\n"
            + "\n".join(body) + "\n"
            + r"\midrule" + "\n"
            + offset_row + "\n"
            + r"\bottomrule" + "\n"
            + r"\end{tabular}" + "\n"
            + r"\end{center}"
        )
        out.append(table)
    return "\n\n".join(out) + "\n"


BIG_TABLE_HEADER = r"""\begin{tabular}{c|cp{3cm}|cp{3cm}|cp{3cm}}
\toprule
    &
    &
    & $\delta_{(s)}-\delta_{(s+1)}$
    & $L_{(s-1)}+V_{(s)}$
    & $(s+1)\Delta_{(s)}$
    & $L_{(s)}/(s+1)$\\
\midrule
$s$ & $\delta_{(s)}$
    & $V_{(s)}$
    & $\Delta_{(s)}$
    & $L_{(s)}$
    & $\alpha_s$
    & ${\hat L}_{(s)}$ \\
\midrule
"""


def emit_stepD_table(rows, grey_V: set[int], grey_L: set[int]) -> str:
    """The 9-row Step-D table.  grey_V / grey_L are the sets of s indices whose
    (delta, V) / (Delta, L, alpha, Lhat) cells are tie-break-dependent."""
    body = []
    for row in rows:
        s = row["s"]
        gv, gl = s in grey_V, s in grey_L
        d_cell = rf"$ {row['delta']} $"
        v_cell = eji_sum_cell(row["V"])
        D_cell = rf"$ {row['Delta']} $"
        L_cell = eji_sum_cell(row["L"])
        a_cell = rf"$ {row['alpha']} $"
        H_cell = eji_sum_cell(row["L"], prefactor_index=s + 1)
        if gv and gl:
            # whole row grey: rowcolor with the s-index cell kept white
            line = (rf"\rowcolor{{gray!15}}\cellcolor{{white}}{s} & {grey(d_cell)} & {grey(v_cell)} & {grey(D_cell)}"
                    f"\n  & {grey(L_cell)}"
                    f"\n    & {grey(a_cell)}"
                    f"\n  & {grey(H_cell)}")
        else:
            def maybe(cell, g):
                return rf"\cellcolor{{gray!15}}{grey(cell)}" if g else cell
            line = (rf"{s} & {maybe(d_cell, gv)} & {maybe(v_cell, gv)} & {maybe(D_cell, gl)}"
                    f"\n  & {maybe(L_cell, gl)}"
                    f"\n    & {maybe(a_cell, gl)}"
                    f"\n  & {maybe(H_cell, gl)}")
        body.append(line + "\n   \\\\[3pt]")

    # Totals block (matches the hand-drawn layout).
    sum_delta = sum(r["delta"] for r in rows)
    sum_Delta = sum(r["Delta"] for r in rows)
    sum_alpha = sum(r["alpha"] for r in rows)
    assert sum_alpha == sum_delta, "alpha-sum must equal delta-sum"
    ne = "" if sum_Delta == sum_delta else rf"\ne {sum_delta}"

    total = np.zeros_like(rows[0]["L"])
    for r in rows:
        total = total + r["delta"] * r["V"]
    check_DL = sum(r["Delta"] * r["L"] for r in rows)
    assert np.array_equal(total, check_DL), "sum dV must equal sum DL"
    total_cell = eji_sum_cell(total)

    totals = (
        "\\midrule\n"
        rf"&\multicolumn{{2}}{{p{{3.5cm}}|}}{{$\sum_s \delta_{{(s)}}={sum_delta}$}}" "\n"
        rf"  & \multicolumn{{2}}{{p{{3.5cm}}|}}{{$\sum_s \Delta_{{(s)}}={sum_Delta}{ne}$}}" "\n"
        rf"  & \multicolumn{{2}}{{p{{3.5cm}}}}{{$\sum_s \alpha_{{(s)}}={sum_alpha}$}}" "\n"
        "  \\\\\n"
        "     \\midrule\n"
        r"& \multicolumn{2}{p{3.5cm}|}{" "\n"
        r"   $\sum_s \delta_{(s)} V_{(s)}=$}" "\n"
        r"  &\multicolumn{2}{p{3.5cm}|}{$" "\n"
        r"  \sum_s \Delta_{(s)} L_{(s)}=$}" "\n"
        r"  &\multicolumn{2}{p{3.5cm}}{$" "\n"
        r"  \sum_s \alpha_{(s)} {\hat L}_{(s)}=$}" "\n"
        "   \\\\\n"
        "&\n&\n   " + total_cell + "\n  &\n  &\n  " + total_cell
        + "\n  &\n  & " + total_cell + "\n  \\\\\n"
    )

    return (TABLE_PREAMBLE + BIG_TABLE_HEADER + "\n".join(body) + "\n"
            + totals + r"\bottomrule" + "\n" + r"\end{tabular}" + "\n"
            + r"\end{center}" + "\n")


def emit_count_matrices(rows, n: int, k: int) -> str:
    """The Eji_LinComb count-matrix table for the alpha != 0 terms.
    Count matrices are displayed transposed (rows = coordinates i,
    columns = vectors j), matching Eji_LinComb's display convention."""
    body = []
    for row in rows:
        if row["alpha"] == 0:
            continue
        s = row["s"]
        mat = row["L"].T  # (k, n): rows i, cols j
        mat_tex = r"\\".join("&".join(str(int(x)) for x in mat[i]) for i in range(k))
        lhat = eji_sum_cell(row["L"], prefactor_index=s + 1).replace(
            "\\newline\n    ", "\\newline\n    $\\phantom{\\tfrac{1}{%d}(}$" % (s + 1))
        body.append(
            rf"{s} & {row['alpha']}" "\n"
            rf"  & $\begin{{pmatrix}}{mat_tex}\end{{pmatrix}}$, idx {s + 1}" "\n"
            rf"  & {lhat} \\[8pt]"
        )
    return (
        "\\begin{center}\n\\small\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\renewcommand{\\arraystretch}{1.2}\n"
        "\\begin{tabular}{cl l p{6.6cm}}\n\\toprule\n"
        "$s$ & $\\alpha_s$ & Eji\\_LinComb structure\n"
        "    &  ${\\hat L}_{(s)}=L_{(s)}/(s{+}1)$ \\\\\n\\midrule\n"
        + "\n".join(body)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{center}\n"
    )


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(exist_ok=True)
    n, k = M.shape

    # Input display.
    cols = ";\\;\n    ".join(
        r"\begin{pmatrix}" + r"\\".join(str(int(M[j][i])) for i in range(k)) + r"\end{pmatrix}"
        for j in range(n))
    write("sx2_input.tex",
          "\\[\n  M \\;=\\; \\left\\{\\;\n    " + cols + "\n  \\;\\right\\}\n"
          rf"  \quad (n={n} \text{{ vectors in }} \R^{k})" + "\n\\]\n")

    # Row minima and balance.
    minima = [int(M[:, i].min()) for i in range(k)]
    write("sx2_minima.tex",
          ", ".join(rf"$m_{i} = {minima[i]}$ (row~{i})" for i in range(k)) + "%\n")
    balance = M - np.array(minima)
    bal_cols = ";\\;\n    ".join(
        r"\begin{pmatrix}" + r"\\".join(str(int(balance[j][i])) for i in range(k)) + r"\end{pmatrix}"
        for j in range(n))
    write("sx2_balance.tex",
          "\\[\n  \\text{balance} = \\left\\{\\;\n    " + bal_cols + "\n  \\;\\right\\}\n\\]\n")

    # Step C mini tables.
    write("sx2_stepC.tex", emit_stepC(M))

    # First Abel summation (per row), merge, and both orderings of the second.
    diffs, offsets = per_row_first_diffs(M)
    for i in range(k):  # offsets must reproduce the Step-B anchors
        off = np.asarray(offsets[i].basis_vec).astype(int)
        assert offsets[i].coeff == minima[i]
        assert np.array_equal(off[:, i], np.ones(n, dtype=int))

    pairs_default = merged_sorted_pairs(diffs)
    pairs_variant = merged_sorted_pairs(diffs, rotate_tie_groups=True)
    m = len(pairs_default)
    assert m == k * (n - 1)

    write("sx2_merged_deltas.tex",
          "\\[\n  " + ",\\;".join(str(int(p[0])) for p in pairs_default) + "\n\\]\n")

    rows_default = second_abel_rows(pairs_default)
    rows_variant = second_abel_rows(pairs_variant)

    # Grey pattern: cells whose values are tie-break-dependent.
    grey_V = {s for s in range(m)
              if not np.array_equal(rows_default[s]["V"], rows_variant[s]["V"])}
    grey_L = {s for s in range(m)
              if not np.array_equal(rows_default[s]["L"], rows_variant[s]["L"])}

    write("sx2_stepD_main.tex", emit_stepD_table(rows_default, grey_V, grey_L))
    write("sx2_stepD_variant.tex", emit_stepD_table(rows_variant, grey_V, grey_L))

    # Verification-box fragments.
    sum_delta = sum(r["delta"] for r in rows_default)
    write("sx2_sum_check.tex",
          "\\[\n"
          rf"  \sum_{{s=0}}^{{{m - 1}}}\delta_{{(s)}} \;=\; "
          + "+".join(str(r["delta"]) for r in rows_default)
          + rf" \;=\; {sum_delta}" + "\n"
          rf"  \;=\; " + "+".join(str(r["alpha"]) for r in rows_default)
          + rf" \;=\; \sum_{{s=0}}^{{{m - 1}}}\alpha_s." + "\n\\]\n")

    total = np.zeros_like(rows_default[0]["L"])
    for r in rows_default:
        total = total + r["delta"] * r["V"]
    assert np.array_equal(total, balance), "totals must reproduce the balance matrix"
    tot_lines = eji_lines(total)
    align = (r"\begin{align*}" "\n"
             r"  \sum_s \delta_{(s)}\,V_{(s)}" "\n"
             r"  \;=\; \sum_s \Delta_{(s)}\,L_{(s)}" "\n"
             r"  \;=\; \sum_s \alpha_s\,{\hat L}_{(s)}" "\n"
             rf"  &= {tot_lines[0]} \\" + "\n"
             + "\n".join(rf"  &\quad{{}}+{line} \\" for line in tot_lines[1:-1]) + "\n"
             + rf"  &\quad{{}}+{tot_lines[-1]}," + "\n"
             r"\end{align*}" "\n")
    mat = total.T  # rows i, cols j
    sm = r"\\".join("&".join(str(int(x)) for x in mat[i]) for i in range(k))
    align += ("i.e.\\ coefficient matrix\n"
              rf"$\left(\begin{{smallmatrix}}{sm}\end{{smallmatrix}}\right)$" + "\n"
              rf"(rows $=i=0,\ldots,{k - 1}$; columns $=j=0,\ldots,{n - 1}$), "
              "matching \\texttt{balance} exactly.\n")
    write("sx2_threeway.tex", align)

    # Count-matrix table.
    write("sx2_count_matrices.tex", emit_count_matrices(rows_default, n, k))

    # Step F display.
    nhash = 2 * m + 1
    out_dim = nhash + k
    alpha_terms = " + ".join(rf"{r['alpha']}\,p_{{{r['s']}}}"
                             for r in rows_default if r["alpha"] != 0)
    minima_str = ",\\;".join(str(v) for v in minima)
    write("sx2_stepF.tex",
          "\\[\n  f(M) \\;=\\;\n"
          rf"  \bigl[\;\underbrace{{{minima_str}}}_{{{k}\text{{ row min.}}}},\;" + "\n"
          rf"          {alpha_terms}\;\bigr]" + "\n"
          rf"  \;\in\; \R^{{{out_dim}}}." + "\n\\]\n")

    print(f"Wrote {len(list(OUT_DIR.glob('sx2_*.tex')))} fragments to {OUT_DIR}")


def write(name: str, content: str):
    header = ("% AUTO-GENERATED by generate_simplex_example_tables.py — do not edit.\n"
              "% Regenerate with: make tables   (or run the script directly).\n")
    (OUT_DIR / name).write_text(header + content)


if __name__ == "__main__":
    main()
