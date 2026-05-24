#!/usr/bin/env python3
"""
symcoder/demo_roundtrip.py
===========================
Round-trip encode → decode demonstration for a small plan.

Context  : electrons = (a, b,c),  muons = (p,q)
Operations: mag, dot (user-defined); euclidean3.dot, euclidean3.eps (library)
Group    : G = S_3 × S_2  (order 12)

Outputs simultaneously to:
  stdout               — plain text
  demo_roundtrip.tex   — LaTeX source; compile with  pdflatex demo_roundtrip.tex

Run from repo root:
  venv/bin/python symcoder/demo_roundtrip.py

Note: Uses Phase 2 with use_complementarity_drop=False so that the alignment
decoder can recover the full repS evaluation table.  The script is intentionally
hacky — it accesses internal attributes of encoder objects to narrate what is
happening, but adds NO new production code.

Use one of 
    use_comp_drop=True
    use_comp_drop=False
at (***) in the code.

"""
from __future__ import annotations

import os, sys, textwrap, hashlib
import numpy as np

# ── locate repo root and activate imports ──────────────────────────────────
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan, repS as _repS_fn
import symcoder.operations.euclidean2
import symcoder.operations.euclidean3
#import symcoder.operations as ops
#from symcoder.operations.euclidean3 import eps as eps3
from symcoder import evaluate, decode_alignment
from symcoder.encoders import (
    OrbitEncoderFactory, SortEncoderFactory, HalfSortEncoderFactory,
    standard_row_pair_factories, OverlapBlockEncoderFactory, Phase2EncoderFactory,
)
from symcoder.encoders.sort_encoder import SortEncoder, HalfSortEncoder

_TEX_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_roundtrip.tex")


# ─────────────────────────────────────────────────────────────────────────────
# TeX macro registry
#
# Operations that carry a tex= payload get a \newcommand whose name is derived
# deterministically from the *content* of the operation (name + tex string),
# not just its name.  This means two operations that share a name but differ in
# their tex= payload (e.g. a user-defined dot and euclidean3.dot) receive
# distinct command names and can coexist in the same .tex file.
#
# Within a session the registry is keyed by id(op) so that different Python
# objects never collide, even if they happen to produce the same command name
# (i.e. same name and same tex= string → same rendered output → same command,
# which is harmless).
# ─────────────────────────────────────────────────────────────────────────────

_TEX_CMD_REGISTRY: dict[int, str] = {}   # id(op) → r"\macroname"


def _tex_cmd_for(op) -> str:
    """Derive a deterministic lowercase \\macroname from the operation's identity.

    The seed is (op.name, op.tex) so that two operations with the same name but
    different tex= payloads produce different command names.  The SHA-256 digest
    is mapped nibble-by-nibble to letters a–p to stay all-alphabetic.
    """
    seed = f"{op.name}\x00{op.tex or ''}"
    digest = hashlib.sha256(seed.encode()).hexdigest()
    letters = "".join(chr(ord("a") + int(c, 16)) for c in digest[:14])
    return "\\" + letters


def _register_operations(operations) -> None:
    """Assign a deterministic TeX command name to every operation that has tex set."""
    for op in operations:
        if op.tex is not None and id(op) not in _TEX_CMD_REGISTRY:
            _TEX_CMD_REGISTRY[id(op)] = _tex_cmd_for(op)


# ─────────────────────────────────────────────────────────────────────────────
# Dual output: write to stdout and TeX file simultaneously
# ─────────────────────────────────────────────────────────────────────────────

class DualOut:
    """Context-manager that writes to stdout and a .tex file."""

    def __init__(self, tex_path: str):
        self._tex = open(tex_path, "w", encoding="utf-8")

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self._tex.write("\n\\end{document}\n")
        self._tex.close()

    # ── plain text helpers ────────────────────────────────────────────────

    def _p(self, text=""):
        print(text)

    def _t(self, tex=""):
        self._tex.write(tex + "\n")

    # ── section / subsection ─────────────────────────────────────────────

    def section(self, title: str):
        bar = "═" * 70
        self._p(f"\n{bar}\n  {title}\n{bar}")
        self._t(f"\n\\section{{{_tex_escape(title)}}}")

    def subsection(self, title: str, tex: str | None = None):
        self._p(f"\n  ── {title} ──")
        tex_title = tex if tex is not None else _tex_escape(title)
        self._t(f"\n\\subsection{{{tex_title}}}")

    def line(self, text: str = "", tex: str | None = None):
        """Print text to stdout, tex (or text) to .tex file."""
        self._p(text)
        self._t((tex if tex is not None else _tex_escape(text)))

    def blank(self):
        self._p()
        self._t()

    def verbatim_line(self, text: str):
        """Output literal text to both (not escaped for TeX)."""
        self._p(text)
        self._t(f"\\texttt{{{_tex_escape(text)}}}" + r"\\")

    def math_line(self, text: str, tex_math: str):
        """stdout: plain text version.  TeX: math-mode version."""
        self._p(f"    {text}")
        self._t(f"\\[ {tex_math} \\]")

    def kv(self, key: str, val_text: str, tex: str | None = None):
        """Key=value line."""
        self._p(f"    {key}: {val_text}")
        vt = tex if tex is not None else _tex_escape(val_text)
        self._t(f"\\textbf{{{_tex_escape(key)}}}: {vt}" + r"\\")

    # ── table helpers ─────────────────────────────────────────────────────

    def begin_table(self, col_spec: str, caption: str = ""):
        self._t(r"\begin{table}[h]\centering")
        if caption:
            self._t(f"\\caption{{{_tex_escape(caption)}}}")
        self._t(f"\\begin{{tabular}}{{{col_spec}}}")
        self._t(r"\toprule")

    def table_row(self, cells_text: list, cells_tex: list | None = None):
        """Print one row.  cells_text → stdout (tab-separated), cells_tex → .tex."""
        self._p("    " + "  |  ".join(str(c) for c in cells_text))
        ct = cells_tex if cells_tex is not None else [_tex_escape(str(c)) for c in cells_text]
        self._t(" & ".join(str(c) for c in ct) + r" \\")

    def table_mid(self):
        self._t(r"\midrule")

    def end_table(self):
        self._t(r"\bottomrule\end{tabular}\end{table}")

    def write_tex_header(self, title: str, operations=()):
        """
        Emit the LaTeX preamble and \\begin{document}.

        If any operation carries a tex= payload, _register_operations is called
        first to assign random command names, then a \\newcommand is emitted for
        each such operation before \\begin{document}.
        """
        _register_operations(operations)
        self._t(textwrap.dedent(r"""
            \documentclass[a4paper,11pt]{article}
            \usepackage{booktabs}
            \usepackage{amsmath,amssymb}
            \usepackage[margin=2cm]{geometry}
            \usepackage{xcolor}
            \usepackage{array}
            \newcommand{\vc}[1]{\ensuremath{\mathbf{#1}}}
        """))
        for op in operations:
            if op.tex is not None:
                cmd = _TEX_CMD_REGISTRY[id(op)]
                self._t(f"\\newcommand{{{cmd}}}[{op.rank}]{{{op.tex}}}")
        self._t(f"\\title{{symcoder round-trip demo\\\\\\large {_tex_escape(title)}}}")
        self._t(r"\date{}\begin{document}\maketitle")


# ─────────────────────────────────────────────────────────────────────────────
# TeX pretty-printing for operators and atoms
# ─────────────────────────────────────────────────────────────────────────────


def _tex_escape(s: str) -> str:
    """Minimal TeX escaping for plain-text strings."""
    return (s.replace("\\", r"\textbackslash ")
             .replace("_", r"\_")
             .replace("&", r"\&")
             .replace("%", r"\%")
             .replace("$", r"\$")
             .replace("#", r"\#")
             .replace("{", r"\{")
             .replace("}", r"\}"))


def _atom_text(atom) -> str:
    """Plain-text representation of an atom, e.g. 'dot(a,b)' or '-eps3(a,b,p)'."""
    sign = "" if atom.sign > 0 else "-"
    return f"{sign}{atom.operation.name}({','.join(atom.labels)})"


def _atom_tex(atom) -> str:
    """TeX math representation of an atom (no surrounding $)."""
    sign = "" if atom.sign > 0 else "-"
    op = atom.operation
    # Prefer the registered \newcommand when the operation carries a tex= payload.
    if op.tex is not None and id(op) in _TEX_CMD_REGISTRY:
        cmd  = _TEX_CMD_REGISTRY[id(op)]
        args = "".join(f"{{{lbl}}}" for lbl in atom.labels)
        return sign + cmd + args
    return _tex_escape(_atom_text(atom))


def _op_tex_sample(op, sample_labels) -> str:
    """TeX snippet for an operation applied to placeholder labels (e.g. for titles)."""
    if op.tex is not None and id(op) in _TEX_CMD_REGISTRY:
        cmd  = _TEX_CMD_REGISTRY[id(op)]
        args = "".join(f"{{{lbl}}}" for lbl in sample_labels[:op.rank])
        return cmd + args
    return _tex_escape(op.name)


def _fmtf(v: float, width: int = 8) -> str:
    """Format a float for plain-text output."""
    return f"{v:+{width}.4f}"


def _tex_num(v: float) -> str:
    return f"{v:+.4f}"


# ─────────────────────────────────────────────────────────────────────────────
# Main demonstration
# ─────────────────────────────────────────────────────────────────────────────

def run():
    # Mymag is defined here to show that user-defined operations work
    # just as well as library ones.  Others comes from the standard library.
    mymag  = Operation("Mymag",  rank=1, odd_parity=False,
                              argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                              eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))),
                              tex=r"|\vc{#1}|")
    # eps is the library singleton from symcoder.operations.euclidean3.
    # Its tex= payload uses \mathbf{} rather than \vc{}, which is fine —
    # each operation's \newcommand is emitted independently.

    # ctx and plan are built before the output context manager so that
    # write_tex_header can register all operation macros from plan.operations.
    electrons = VectorType("electrons", ("a", "b"))
    muons     = VectorType("muons",     ("p",))
    ctx       = Context((electrons, muons))
    plan      = Plan(context=ctx, operations=(mymag,
              symcoder.operations.euclidean3.dot, 
              symcoder.operations.euclidean2.dot, 
              symcoder.operations.euclidean3.eps,
              ))

    with DualOut(_TEX_FILE) as out:
        ctx_title = "electrons=(a,b), muons=(p)"
        out.write_tex_header(ctx_title, operations=plan.operations)

        # ── Setup ─────────────────────────────────────────────────────────
        out.section("Setup")

        out.kv("Particle types",
               "electrons=(a,b),  muons=(p)",
               r"electrons $= \{a, b\}$, \quad muons $= \{p\}$")
        out.kv("Group G",
               f"S_electrons × S_muons = S_2 × S_1  (order {ctx.the_group.order()})",
               rf"$G = S_{{electrons}} \times S_{{muons}} = S_2 \times S_1$,"
               rf" order $= {ctx.the_group.order()}$")

        out.blank()
        out.line("Operations:\n")
        _pl = list("xyzw")
        for op in plan.operations:
            pl   = _pl[:op.rank]
            args = ",".join(pl)
            sym  = op.argument_symmetry.name.lower()
            args_vc = ",".join(rf"\vc{{{p}}}" for p in pl)
            out.line(
                f"  {op.name}({args})  (rank {op.rank}, {sym})",
                tex=rf"${_op_tex_sample(op, pl)}$ \quad (rank {op.rank}, {sym})\\"
            )

        # Fixed event (3-D vectors; eps3 needs ≥3 dimensions)
        event = {
            "a": np.array([ 1.2,  0.3,   0.5]),
            "b": np.array([ 0.4,  1.1,   0.2]),
            "c": np.array([ -0.5, 1.2,  -0.7]),
            "p": np.array([ 0.0,  0.0,  +1.0]),
            "q": np.array([ 0.0,  0.0,  -1.0]),
        }

        out.blank()
        out.line("Event vectors (3-D):", tex=r"\textbf{Event vectors (3-D):}\\")
        for lbl, vec in event.items():
            vec_text = f"[{', '.join(f'{v:.2f}' for v in vec)}]"
            vec_tex  = r",\ ".join(f"{v:.2f}" for v in vec)
            out.kv(f"  {lbl}", vec_text, tex=rf"$\vc{{{lbl}}} = ({vec_tex})$\\")

        # Factories
        # Toggle this flag to see the effect of complementarity drop:
        #   False → every PairFlavour is explicitly encoded (larger output, easier to read)
        #   True  → the largest non-self PairFlavour is dropped and reconstructed at decode
        use_comp_drop = False # (***) - see reference in docstring at top of file.

        orbit_factory   = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
        phase2_factory  = Phase2EncoderFactory([
            OverlapBlockEncoderFactory(
                standard_row_pair_factories(), use_complementarity_drop=use_comp_drop
            )
        ])

        # Build encoders and encode
        orbit_enc  = orbit_factory.build(plan)
        phase2_enc = phase2_factory.build(plan)
        phase1_vals = orbit_enc.encode(event).values
        phase2_vals = phase2_enc.encode(event).values

        fo_list    = _repS_fn(ctx, plan.operations)
        repS_atoms = [fo.canonical_representative() for fo in fo_list]

        # ── Phase 1 encoding ──────────────────────────────────────────────
        out.section("Phase 1 Encoding  (one orbit per FlavouredOperator in repS)")
        out.line(f"  {len(orbit_enc._selections)} FlavouredOperators in repS")
        out.blank()

        phase1_results = {}
        cursor = 0
        for fo, enc in orbit_enc._selections:
            key    = (fo.operation, tuple(fo.flavour.counts))
            chunk  = phase1_vals[cursor : cursor + enc.output_dim]
            decoded = enc.decode(chunk)
            phase1_results[key] = decoded

            # Orbit atoms
            if isinstance(enc, HalfSortEncoder):
                orbit_atoms = enc._representatives
                method = "HalfSortEncoder  (stores |eval| of sign=+1 reps)"
                method_tex = r"\texttt{HalfSortEncoder} — stores $|\text{eval}|$ of sign$=+1$ representatives"
            else:
                orbit_atoms = enc._orbit
                method = "SortEncoder  (stores eval of full orbit)"
                method_tex = r"\texttt{SortEncoder} — stores eval of full orbit"

            _sl = list("xyzw")[:fo.operation.rank]
            out.subsection(
                f"FlavouredOperator: {fo.operation.name}  flavour={tuple(fo.flavour.counts)}",
                tex=f"FlavouredOperator: ${_op_tex_sample(fo.operation, _sl)}$  flavour $= {tuple(fo.flavour.counts)}$"
            )
            out.kv("  Canonical representative", _atom_text(fo.canonical_representative()),
                   tex=f"$\\displaystyle{_atom_tex(fo.canonical_representative())}$")
            out.kv("  Encoder", method, tex=method_tex + r"\\")
            out.kv("  Orbit atoms", "  ".join(_atom_text(a) for a in orbit_atoms),
                   tex=r",\ ".join(f"${_atom_tex(a)}$" for a in orbit_atoms))

            eval_vals = [evaluate(a, event) for a in orbit_atoms]
            out.kv("  Eval values",
                   "  ".join(_fmtf(v) for v in eval_vals),
                   tex=r",\ ".join(_tex_num(v) for v in eval_vals))

            out.kv(f"  Encoded [{chunk.size} reals]",
                   "  ".join(_fmtf(v) for v in chunk),
                   tex=r",\ ".join(_tex_num(v) for v in chunk))

            cursor += enc.output_dim

        # ── Phase 2 encoding ──────────────────────────────────────────────
        _drop_label = "enabled" if use_comp_drop else "disabled"
        out.section(f"Phase 2 Encoding  (OverlapBlocks, complementarity-drop {_drop_label})")

        phase2_results = {}   # pf → decoded AnnotatedMultisetOfRealPairs list
        p2cursor = 0
        for spec, block_enc in phase2_enc._block_encoders:
            pf0 = block_enc._selections[0].pf
            block_title = (f"{pf0.op_u.name}[{tuple(pf0.flavour_u.counts)}] "
                           f"× {pf0.op_v.name}[{tuple(pf0.flavour_v.counts)}]")
            _qs = list("xyzw")  # placeholder labels for block title
            block_title_tex = (f"${_op_tex_sample(pf0.op_u, _qs)}"
                               f" \\times "
                               f"{_op_tex_sample(pf0.op_v, _qs)}$"
                               f"  (block)")

            out.subsection(f"Block: {block_title}", tex=f"Block: {block_title_tex}")

            block_slice = phase2_vals[p2cursor : p2cursor + block_enc.output_dim]
            decoded_list = block_enc.decode(block_slice, phase1_results)

            sel_cursor = 0
            for sel, decoded in zip(block_enc._selections, decoded_list):
                pf   = sel.pf
                kind = "NULL_SELF" if sel.encoder.output_dim == 0 and sel.encoder.method_name == "null_self" \
                       else ("NULL_COMP" if sel.is_comp_drop else "ASSOC")
                enc_dim  = sel.encoder.output_dim
                method   = sel.encoder.method_name
                # NULL_COMP slots contribute 0 bytes to block_slice (not stored);
                # enc_dim here is the *original* encoder's dim (informative: shows space saved).
                stored_dim = 0 if sel.is_comp_drop else enc_dim

                out.line(f"    PairFlavour  overlap={tuple(pf.overlap)}"
                         f"  type={method}  kind={kind}  dim={enc_dim}")
                out._t(f"\\textbf{{PairFlavour}} overlap$={tuple(pf.overlap)}$,"
                       f" type=\\texttt{{{_tex_escape(method)}}}, kind=\\texttt{{{_tex_escape(kind)}}},"
                       f" encoded dim$={enc_dim}$\\\\")

                sel_vals = block_slice[sel_cursor : sel_cursor + stored_dim]
                if stored_dim > 0:
                    out.kv("      Encoded values",
                           "  ".join(_fmtf(v) for v in sel_vals),
                           tex=r",\ ".join(_tex_num(v) for v in sel_vals))
                elif sel.is_comp_drop:
                    out.line("      [not stored — reconstructed from complement]",
                             tex=r"\textit{not stored --- reconstructed from complement}\\")


                if decoded.pairs:
                    out.line("      Decoded pairs  (u_val, v_val):",
                             tex=r"\quad Decoded pairs $(u, v)$:\\")
                    for ap, pair in zip(decoded.atom_pairs[:8], decoded.pairs[:8]):
                        u_atom, v_atom = ap
                        u_val,  v_val  = pair
                        out.line(f"        ({_atom_text(u_atom)}, {_atom_text(v_atom)})"
                                 f"  →  ({_fmtf(u_val)}, {_fmtf(v_val)})",
                                 tex=(f"$({_atom_tex(u_atom)},\\ {_atom_tex(v_atom)})$"
                                      f"$\\to ({_tex_num(u_val)},\\ {_tex_num(v_val)})$\\\\"))
                    if len(decoded.pairs) > 8:
                        out.line(f"        ... ({len(decoded.pairs)-8} more pairs)")

                sel_cursor += stored_dim

            p2cursor += block_enc.output_dim

        # Rebuild all_pair_decoded list for alignment decoder
        all_pair_decoded = []
        p2cursor = 0
        for _, block_enc in phase2_enc._block_encoders:
            block_slice  = phase2_vals[p2cursor : p2cursor + block_enc.output_dim]
            decoded_list = block_enc.decode(block_slice, phase1_results)
            all_pair_decoded.extend(decoded_list)
            p2cursor += block_enc.output_dim

        # ── Alignment decoding ────────────────────────────────────────────
        out.section("Alignment Decoder  (recovers the G-orbit of repS evaluations)")

        decoded_align = decode_alignment(plan, phase1_results, all_pair_decoded, atol=1e-10)

        n_cols = ctx.the_group.order()
        n_reps = len(repS_atoms)
        out.kv("  |repS|", str(n_reps), tex=rf"$|repS| = {n_reps}$\\")
        out.kv("  |G|",    str(n_cols), tex=rf"$|G| = {n_cols}$\\")
        out.kv("  Decoded vectors", str(len(decoded_align.vectors)),
               tex=rf"decoded column vectors: ${len(decoded_align.vectors)}$\\")

        # Ground truth — TheGroup is sole authority on group action
        the_group = ctx.the_group
        gt_vectors = [
            tuple(evaluate(g.apply(atom), event) for atom in repS_atoms)
            for g in the_group.all_group_elements()
        ]

        # Show the alignment table
        out.blank()
        out.line("  Decoded alignment table  (rows=repS atoms, cols=group elements):")
        out._t(r"\textbf{Decoded alignment table} (rows$=\text{repS}$ atoms, columns$=g\in G$):\\")

        # Column headers
        col_spec = "l" + "r" * n_cols
        col_hdrs_text = ["atom"] + [f"col {c}" for c in range(n_cols)]
        col_hdrs_tex  = ["Atom"]  + [f"$g_{{{c}}}$" for c in range(n_cols)]
        out.begin_table(col_spec, caption="Decoded alignment table")
        out.table_row(col_hdrs_text, col_hdrs_tex)
        out.table_mid()

        # One decoded vector per column  (decoded_align.vectors is a list of n_cols tuples)
        # Transpose: rows are atoms, columns are group elements
        if decoded_align.vectors:
            for r, atom in enumerate(repS_atoms):
                row_vals = [vec[r] for vec in decoded_align.vectors]
                cells_text = [_atom_text(atom)] + [_fmtf(v, 7) for v in row_vals]
                cells_tex  = [f"${_atom_tex(atom)}$"] + [_tex_num(v) for v in row_vals]
                out.table_row(cells_text, cells_tex)

        out.end_table()

        # Ground truth table
        out.blank()
        out.line("  Ground truth  (direct evaluation of g·repS on event):")
        out._t(r"\textbf{Ground truth} (direct evaluation of $g \cdot \text{repS}$ on event):\\")
        out.begin_table(col_spec, caption="Ground truth alignment table")
        out.table_row(col_hdrs_text, col_hdrs_tex)
        out.table_mid()
        for r, atom in enumerate(repS_atoms):
            gt_row = [vec[r] for vec in gt_vectors]
            cells_text = [_atom_text(atom)] + [_fmtf(v, 7) for v in gt_row]
            cells_tex  = [f"${_atom_tex(atom)}$"] + [_tex_num(v) for v in gt_row]
            out.table_row(cells_text, cells_tex)
        out.end_table()

        # Comparison
        out.blank()
        out.line("  Multiset comparison  (decoded vs ground truth):")
        out._t(r"\textbf{Multiset comparison:}\\")

        dec_sorted = sorted(decoded_align.vectors)
        gt_sorted  = sorted(gt_vectors)
        atol = 1e-8
        ok = (len(dec_sorted) == len(gt_sorted) and
              all(max(abs(a - b) for a, b in zip(da, ga)) <= atol
                  for da, ga in zip(dec_sorted, gt_sorted)))

        status_text = "✓ MATCH  (multisets equal within atol=1e-8)" if ok \
                 else "✗ MISMATCH — decoded multiset differs from ground truth"
        status_tex  = r"\textcolor{green!60!black}{\textbf{MATCH}} (multisets equal, atol$=10^{-8}$)" if ok \
                 else r"\textcolor{red}{\textbf{MISMATCH}} — decoded multiset differs from ground truth"
        out.line(f"  {status_text}", tex=status_tex)

        # Per-atom residuals
        out.blank()
        out.line("  Max |error| per repS row (decoded vs ground truth):")
        out._t(r"\textbf{Max $|$error$|$ per repS row:}\\")
        for r, atom in enumerate(repS_atoms):
            dec_row = sorted(vec[r] for vec in decoded_align.vectors)
            gt_row  = sorted(vec[r] for vec in gt_vectors)
            max_err = max(abs(d - g) for d, g in zip(dec_row, gt_row))
            out.kv(f"    {_atom_text(atom)}",
                   f"max_err = {max_err:.2e}",
                   tex=f"${_atom_tex(atom)}$: max err $= {max_err:.2e}$\\\\")

        # ── Summary ───────────────────────────────────────────────────────
        out.section("Summary")
        phase1_dim = orbit_enc.output_dim
        phase2_dim = phase2_enc.output_dim
        out.kv("  Phase 1 encoded length", str(phase1_dim),
               tex=rf"Phase 1 encoded length: ${phase1_dim}$ reals\\")
        out.kv("  Phase 2 encoded length", str(phase2_dim),
               tex=rf"Phase 2 encoded length: ${phase2_dim}$ reals\\")
        out.kv("  Total encoded length",   str(phase1_dim + phase2_dim),
               tex=rf"Total encoded length: ${phase1_dim + phase2_dim}$ reals\\")
        out.kv("  repS size",              str(n_reps),
               tex=rf"$|\text{{repS}}| = {n_reps}$ atoms\\")
        out.kv("  |G|",                    str(n_cols),
               tex=rf"$|G| = {n_cols}$ group elements\\")
        out.kv("  Decoded vectors",        str(len(decoded_align.vectors)),
               tex=rf"decoded alignment vectors: ${len(decoded_align.vectors)}$\\")
        out.line(f"\n  TeX output written to: {_TEX_FILE}")
        out.line(  "  Compile with:  pdflatex demo_roundtrip.tex")


if __name__ == "__main__":
    run()
