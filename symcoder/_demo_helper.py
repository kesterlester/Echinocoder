"""
symcoder._demo_helper
=====================
Shared infrastructure for the demo scripts in this directory
(``demo_roundtrip.py``, ``demo_shadow_event.py``).

Private (underscore-prefixed) — end users have no business importing
from here.  It exists only to keep the two demo scripts from drifting
apart whenever the demo narrative is improved.

Public API
----------
* ``DualOut(tex_path)``         — context-managed dual writer (stdout + .tex)
* ``run_demo(plan, ctx, event, tex_path, *, title_main, subtitle,
             intro_callback=None, use_comp_drop=False)`` — render the full
  encode → decode → explain narrative for one ``(plan, event)`` pair.

The TeX preamble, operation-macro registry, atom/operator pretty-printing,
phase-1 / phase-2 encoding explanations, alignment-decoder narration, and
summary section all live here.
"""
from __future__ import annotations

import hashlib
import os
import textwrap
from typing import Callable, Optional

import numpy as np

from symatom import repS as _repS_fn
from symcoder import evaluate, decode_alignment
from symcoder.encoders import (
    OrbitEncoderFactory, SortEncoderFactory, HalfSortEncoderFactory,
    standard_row_pair_factories, OverlapBlockEncoderFactory, Phase2EncoderFactory,
)
from symcoder.encoders.sort_encoder import SortEncoder, HalfSortEncoder


# ─────────────────────────────────────────────────────────────────────────────
# TeX macro registry
#
# Operations that carry a tex= payload get a \newcommand whose name is derived
# deterministically from the *content* of the operation (name + tex string),
# not just its name.  This means two operations that share a name but differ
# in their tex= payload (e.g. a user-defined dot and euclidean3.dot) receive
# distinct command names and can coexist in the same .tex file.
# ─────────────────────────────────────────────────────────────────────────────

_TEX_CMD_REGISTRY: dict[int, str] = {}


def _tex_cmd_for(op) -> str:
    """Derive a deterministic lowercase \\macroname from the operation's identity.

    Seed is (op.name, op.tex) so two operations with the same name but different
    tex= payloads produce different command names.  Digest letters are mapped
    a–p so the macro is all-alphabetic.
    """
    seed = f"{op.name}\x00{op.tex or ''}"
    digest = hashlib.sha256(seed.encode()).hexdigest()
    letters = "".join(chr(ord("a") + int(c, 16)) for c in digest[:14])
    return "\\" + letters


def _register_operations(operations) -> None:
    for op in operations:
        if op.tex is not None and id(op) not in _TEX_CMD_REGISTRY:
            _TEX_CMD_REGISTRY[id(op)] = _tex_cmd_for(op)


# ─────────────────────────────────────────────────────────────────────────────
# TeX string utilities (module-private, but used by run_demo and any
# caller-supplied intro_callback)
# ─────────────────────────────────────────────────────────────────────────────


def tex_escape(s: str) -> str:
    """Minimal TeX escaping for plain-text strings."""
    return (s.replace("\\", r"\textbackslash ")
             .replace("_", r"\_")
             .replace("&", r"\&")
             .replace("%", r"\%")
             .replace("$", r"\$")
             .replace("#", r"\#")
             .replace("{", r"\{")
             .replace("}", r"\}"))


def atom_text(atom) -> str:
    sign = "" if atom.sign > 0 else "-"
    return f"{sign}{atom.operation.name}({','.join(atom.labels)})"


def atom_tex(atom) -> str:
    sign = "" if atom.sign > 0 else "-"
    op = atom.operation
    if op.tex is not None and id(op) in _TEX_CMD_REGISTRY:
        cmd  = _TEX_CMD_REGISTRY[id(op)]
        args = "".join(f"{{{lbl}}}" for lbl in atom.labels)
        return sign + cmd + args
    return tex_escape(atom_text(atom))


def op_tex_sample(op, sample_labels) -> str:
    if op.tex is not None and id(op) in _TEX_CMD_REGISTRY:
        cmd  = _TEX_CMD_REGISTRY[id(op)]
        args = "".join(f"{{{lbl}}}" for lbl in sample_labels[:op.rank])
        return cmd + args
    return tex_escape(op.name)


def fmtf(v: float, width: int = 8) -> str:
    return f"{v:+{width}.4f}"


def tex_num(v: float) -> str:
    return f"{v:+.4f}"


# ─────────────────────────────────────────────────────────────────────────────
# Dual output: stdout and a .tex file simultaneously
# ─────────────────────────────────────────────────────────────────────────────


class DualOut:
    """Context-managed writer that emits plain text to stdout and TeX to a file.

    The TeX file is opened on construction, closed (with ``\\end{document}``
    appended) on context-manager exit.
    """

    def __init__(self, tex_path: str):
        self._tex = open(tex_path, "w", encoding="utf-8")

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self._tex.write("\n\\end{document}\n")
        self._tex.close()

    # ── plain text and TeX writers ────────────────────────────────────────

    def _p(self, text: str = "") -> None:
        print(text)

    def _t(self, tex: str = "") -> None:
        self._tex.write(tex + "\n")

    # ── section / subsection helpers ─────────────────────────────────────

    def section(self, title: str) -> None:
        bar = "═" * 70
        self._p(f"\n{bar}\n  {title}\n{bar}")
        self._t(f"\n\\section{{{tex_escape(title)}}}")

    def subsection(self, title: str, tex: Optional[str] = None) -> None:
        self._p(f"\n  ── {title} ──")
        tex_title = tex if tex is not None else tex_escape(title)
        self._t(f"\n\\subsection{{{tex_title}}}")

    def line(self, text: str = "", tex: Optional[str] = None) -> None:
        self._p(text)
        self._t((tex if tex is not None else tex_escape(text)))

    def blank(self) -> None:
        self._p()
        self._t()

    def verbatim_line(self, text: str) -> None:
        self._p(text)
        self._t(f"\\texttt{{{tex_escape(text)}}}" + r"\\")

    def math_line(self, text: str, tex_math: str) -> None:
        self._p(f"    {text}")
        self._t(f"\\[ {tex_math} \\]")

    def kv(self, key: str, val_text: str, tex: Optional[str] = None) -> None:
        self._p(f"    {key}: {val_text}")
        vt = tex if tex is not None else tex_escape(val_text)
        self._t(f"\\textbf{{{tex_escape(key)}}}: {vt}" + r"\\")

    # ── table helpers ────────────────────────────────────────────────────

    def begin_table(self, col_spec: str, caption: str = "") -> None:
        self._t(r"\begin{table}[h]\centering")
        if caption:
            self._t(f"\\caption{{{tex_escape(caption)}}}")
        self._t(f"\\begin{{tabular}}{{{col_spec}}}")
        self._t(r"\toprule")

    def table_row(self, cells_text: list, cells_tex: Optional[list] = None) -> None:
        self._p("    " + "  |  ".join(str(c) for c in cells_text))
        ct = cells_tex if cells_tex is not None else [tex_escape(str(c)) for c in cells_text]
        self._t(" & ".join(str(c) for c in ct) + r" \\")

    def table_mid(self) -> None:
        self._t(r"\midrule")

    def end_table(self) -> None:
        self._t(r"\bottomrule\end{tabular}\end{table}")

    # ── document preamble ────────────────────────────────────────────────

    def write_tex_header(self, *, title_main: str, subtitle: str, operations=()) -> None:
        """Emit the LaTeX preamble and \\begin{document}.

        ``title_main`` is the main demo title (e.g. ``symcoder round-trip
        demo`` or ``symcoder parity-blindness demo``).  ``subtitle`` is the
        ``\\large`` line that follows (e.g. ``electrons=(a,b,c), muons=(p,q)``
        or ``shadow event #0,  parity = +1``).

        Any operation that carries a ``tex=`` payload gets a ``\\newcommand``
        emitted before ``\\begin{document}`` so it can be used downstream.
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
        self._t(f"\\title{{{tex_escape(title_main)}\\\\\\large {tex_escape(subtitle)}}}")
        self._t(r"\date{}\begin{document}\maketitle")


# ─────────────────────────────────────────────────────────────────────────────
# Shared demo body
# ─────────────────────────────────────────────────────────────────────────────


# Sentinel for "caller didn't pass a factory; build the standard one".
# Distinct from explicit ``None`` (which means "skip that phase entirely").
_DEFAULT_FACTORY = object()


def run_demo(
    plan,
    ctx,
    event: dict,
    tex_path: str,
    *,
    title_main: str,
    subtitle: str,
    intro_callback: Optional[Callable[[DualOut], None]] = None,
    use_comp_drop: bool = False,
    orbit_factory=_DEFAULT_FACTORY,
    phase2_factory=_DEFAULT_FACTORY,
    phase3_factory=None,
) -> None:
    """Render an encode → decode → explain narrative to ``tex_path``.

    The narrative is *phase-conditional*: each of Phase 1 (orbit), Phase 2
    (overlap-block pair), and Phase 3 (whole-table) is narrated only when
    its factory argument is present.  The Alignment-Decoder section is
    additionally gated on Phase 1 AND Phase 2 being present, since the
    alignment decoder is fed from both.

    Parameters
    ----------
    plan, ctx
        As built by the caller.  ``plan.operations`` and ``ctx.types`` are
        the sources of truth for everything the demo describes.
    event
        ``{label: ndarray}`` dict mapping particle labels to numerical
        vectors.  Must cover every label in every ``ctx.types[*].labels``.
    tex_path
        Destination ``.tex`` filename.  Will be overwritten.
    title_main, subtitle
        Document title parts (main heading + ``\\large`` subtitle line).
    intro_callback
        Optional callback ``cb(out)``.  If given, it is invoked *after* the
        TeX header (so it can call ``out.section`` / ``out.kv`` / ``out.line``
        etc.) and *before* the standard "Setup" section.  Used by the
        parity-blindness demo to emit an "About this document" preamble.
    use_comp_drop
        Whether the Phase 2 OverlapBlock factory should drop the
        complementarity-implied selection.  Only affects the *default*
        Phase 2 factory; a caller-supplied ``phase2_factory`` ignores this
        flag.  ``False`` produces the most explicit/readable Phase 2
        output; ``True`` produces the most compact.
    orbit_factory
        Phase 1 sub-factory chain.  Default (sentinel): build the
        standard ``OrbitEncoderFactory([HalfSortEncoderFactory(),
        SortEncoderFactory()])``.  Pass ``None`` explicitly to *skip*
        Phase 1 narration entirely.
    phase2_factory
        Phase 2 factory.  Default (sentinel): build the standard
        ``Phase2EncoderFactory([OverlapBlockEncoderFactory(...)])``.
        Pass ``None`` explicitly to skip Phase 2 narration.
    phase3_factory
        Phase 3 factory (e.g.\ a stacking ``Phase3EncoderFactory`` or
        a single sub-factory with ``.build(plan)``).  Default ``None``,
        meaning Phase 3 is not narrated.  Pass a factory instance to
        include the Phase 3 section.
    """
    # Resolve factory sentinels.  An explicit ``None`` from the caller is
    # treated as "skip"; the sentinel means "build the standard one here".
    if orbit_factory is _DEFAULT_FACTORY:
        orbit_factory = OrbitEncoderFactory(
            [HalfSortEncoderFactory(), SortEncoderFactory()]
        )
    if phase2_factory is _DEFAULT_FACTORY:
        phase2_factory = Phase2EncoderFactory([
            OverlapBlockEncoderFactory(
                standard_row_pair_factories(),
                use_complementarity_drop=use_comp_drop,
            )
        ])
    with DualOut(tex_path) as out:
        # ── Derive setup display strings from ctx ────────────────────────
        ctx_title = ",  ".join(
            f"{t.name}=({','.join(t.labels)})" for t in ctx.types
        )
        ptypes_text = ctx_title
        ptypes_tex  = r",\quad ".join(
            f"{tex_escape(t.name)} $= \\{{{', '.join(t.labels)}\\}}$"
            for t in ctx.types
        )
        g_parts_name = " × ".join(f"S_{t.name}" for t in ctx.types)
        g_parts_n    = " × ".join(f"S_{len(t.labels)}" for t in ctx.types)
        g_text = (f"G = {g_parts_name} = {g_parts_n}"
                  f"  (order {ctx.the_group.order()})")
        g_parts_name_tex = r" \times ".join(
            rf"S_{{\text{{{tex_escape(t.name)}}}}}" for t in ctx.types
        )
        g_parts_n_tex = r" \times ".join(
            rf"S_{{{len(t.labels)}}}" for t in ctx.types
        )
        g_tex = (rf"$G = {g_parts_name_tex} = {g_parts_n_tex}$,"
                 rf" order $= {ctx.the_group.order()}$")

        out.write_tex_header(title_main=title_main, subtitle=subtitle,
                             operations=plan.operations)

        # Caller-supplied preamble (e.g. "About this document").
        if intro_callback is not None:
            intro_callback(out)

        # ── Setup ────────────────────────────────────────────────────────
        out.section("Setup")

        out.kv("Particle types", ptypes_text, ptypes_tex)
        out.kv("Group G",        g_text,      g_tex)

        out.blank()
        out.line("Operations:\n")
        _pl = list("xyzw")
        for op in plan.operations:
            pl   = _pl[:op.rank]
            args = ",".join(pl)
            sym  = op.argument_symmetry.name.lower()
            out.line(
                f"  {op.name}({args})  (rank {op.rank}, {sym})",
                tex=rf"${op_tex_sample(op, pl)}$ \quad (rank {op.rank}, {sym})\\"
            )

        out.blank()
        out.line("Event vectors (3-D):", tex=r"\textbf{Event vectors (3-D):}\\")
        for lbl, vec in event.items():
            vec_text = f"[{', '.join(f'{v:.2f}' for v in vec)}]"
            vec_tex  = r",\ ".join(f"{v:.2f}" for v in vec)
            out.kv(f"  {lbl}", vec_text, tex=rf"$\vc{{{lbl}}} = ({vec_tex})$\\")

        # ── Build encoders for each active phase ─────────────────────────
        orbit_enc   = orbit_factory.build(plan)  if orbit_factory  is not None else None
        phase2_enc  = phase2_factory.build(plan) if phase2_factory is not None else None
        phase3_enc  = phase3_factory.build(plan) if phase3_factory is not None else None

        phase1_vals = orbit_enc.encode(event).values  if orbit_enc  is not None else None
        phase2_vals = phase2_enc.encode(event).values if phase2_enc is not None else None
        phase3_vals = phase3_enc.encode(event).values if phase3_enc is not None else None

        fo_list    = _repS_fn(ctx, plan.operations)
        repS_atoms = [fo.canonical_representative() for fo in fo_list]
        phase1_results: dict = {}

        # ── Phase 1 encoding ─────────────────────────────────────────────
        if orbit_enc is not None:
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

                if isinstance(enc, HalfSortEncoder):
                    orbit_atoms = enc._representatives
                    method = "HalfSortEncoder  (stores |eval| of sign=+1 reps)"
                    method_tex = r"\texttt{HalfSortEncoder} — stores $|\text{eval}|$ of sign$=+1$ representatives"
                elif isinstance(enc, SortEncoder):
                    orbit_atoms = enc._orbit
                    method = "SortEncoder  (stores eval of full orbit)"
                    method_tex = r"\texttt{SortEncoder} — stores eval of full orbit"
                else:
                    orbit_atoms = enc._orbit
                    method = f"[{enc.method_name}] Encoder  (working method not known)"
                    method_tex = r""

                _sl = list("xyzw")[:fo.operation.rank]
                out.subsection(
                    f"FlavouredOperator: {fo.operation.name}  flavour={tuple(fo.flavour.counts)}",
                    tex=f"FlavouredOperator: ${op_tex_sample(fo.operation, _sl)}$  flavour $= {tuple(fo.flavour.counts)}$"
                )
                out.kv("  Canonical representative", atom_text(fo.canonical_representative()),
                       tex=f"$\\displaystyle{atom_tex(fo.canonical_representative())}$")
                out.kv("  Encoder", method, tex=method_tex + r"\\")
                out.kv("  Orbit atoms", "  ".join(atom_text(a) for a in enc._orbit),
                       tex=r",\ ".join(f"${atom_tex(a)}$" for a in enc._orbit))
                if isinstance(enc, HalfSortEncoder):
                    out.kv("  represented by ", "  ".join(atom_text(a) for a in orbit_atoms),
                       tex=r",\ ".join(f"${atom_tex(a)}$" for a in orbit_atoms))

                eval_vals = [evaluate(a, event) for a in orbit_atoms]
                out.kv("  Eval values",
                       "  ".join(fmtf(v) for v in eval_vals),
                       tex=r",\ ".join(tex_num(v) for v in eval_vals))

                out.kv(f"  Encoded [{chunk.size} reals]",
                       "  ".join(fmtf(v) for v in chunk),
                       tex=r",\ ".join(tex_num(v) for v in chunk))

                cursor += enc.output_dim

        # ── Phase 2 encoding ─────────────────────────────────────────────
        if phase2_enc is not None:
            _drop_label = "enabled" if use_comp_drop else "disabled"
            out.section(f"Phase 2 Encoding  (OverlapBlocks, complementarity-drop {_drop_label})")

            p2cursor = 0
            for spec, block_enc in phase2_enc._block_encoders:
                pf0 = block_enc._selections[0].pf
                block_title = (f"{pf0.op_u.name}[{tuple(pf0.flavour_u.counts)}] "
                               f"× {pf0.op_v.name}[{tuple(pf0.flavour_v.counts)}]")
                _qs = list("xyzw")
                block_title_tex = (f"${op_tex_sample(pf0.op_u, _qs)}"
                                   f" \\times "
                                   f"{op_tex_sample(pf0.op_v, _qs)}$"
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
                    stored_dim = 0 if sel.is_comp_drop else enc_dim

                    out.line(f"    PairFlavour  overlap={tuple(pf.overlap)}"
                             f"  type={method}  kind={kind}  dim={enc_dim}")
                    out._t(f"\\textbf{{PairFlavour}} overlap$={tuple(pf.overlap)}$,"
                           f" type=\\texttt{{{tex_escape(method)}}}, kind=\\texttt{{{tex_escape(kind)}}},"
                           f" encoded dim$={enc_dim}$\\\\")

                    sel_vals = block_slice[sel_cursor : sel_cursor + stored_dim]
                    if stored_dim > 0:
                        out.kv("      Encoded values",
                               "  ".join(fmtf(v) for v in sel_vals),
                               tex=r",\ ".join(tex_num(v) for v in sel_vals))
                    elif sel.is_comp_drop:
                        out.line("      [not stored — reconstructed from complement]",
                                 tex=r"\textit{not stored --- reconstructed from complement}\\")

                    _emit_decoded_pair_multiset(out, decoded)

                    sel_cursor += stored_dim

                p2cursor += block_enc.output_dim

        # ── Phase 3 encoding (one-shot whole-table) ──────────────────────
        if phase3_enc is not None:
            _emit_phase3_section(out, plan, ctx, event, phase3_enc, phase3_vals)

        # ── Alignment decoder (needs Phase 1 + Phase 2 outputs together) ──
        if orbit_enc is not None and phase2_enc is not None:
            # Rebuild all_pair_decoded list for alignment decoder
            all_pair_decoded = []
            p2cursor = 0
            for _, block_enc in phase2_enc._block_encoders:
                block_slice  = phase2_vals[p2cursor : p2cursor + block_enc.output_dim]
                decoded_list = block_enc.decode(block_slice, phase1_results)
                all_pair_decoded.extend(decoded_list)
                p2cursor += block_enc.output_dim

            out.section("Alignment Decoder  (recovers the G-orbit of repS evaluations)")

            decoded_align = decode_alignment(plan, phase1_results, all_pair_decoded, atol=1e-10)

            n_cols = ctx.the_group.order()
            n_reps = len(repS_atoms)
            out.kv("  |repS|", str(n_reps), tex=rf"$|repS| = {n_reps}$\\")
            out.kv("  |G|",    str(n_cols), tex=rf"$|G| = {n_cols}$\\")
            out.kv("  Decoded vectors", str(len(decoded_align.vectors)),
                   tex=rf"decoded column vectors: ${len(decoded_align.vectors)}$\\")

            the_group = ctx.the_group
            gt_vectors = [
                tuple(evaluate(g.apply(atom), event) for atom in repS_atoms)
                for g in the_group.all_group_elements()
            ]

            out.blank()
            out.line("  Decoded alignment table  (rows=repS atoms, cols=group elements):")
            out._t(r"\textbf{Decoded alignment table} (rows$=\text{repS}$ atoms, columns$=g\in G$):\\")

            col_spec = "l" + "r" * n_cols
            col_hdrs_text = ["atom"] + [f"col {c}" for c in range(n_cols)]
            col_hdrs_tex  = ["Atom"]  + [f"$g_{{{c}}}$" for c in range(n_cols)]
            out.begin_table(col_spec, caption="Decoded alignment table")
            out.table_row(col_hdrs_text, col_hdrs_tex)
            out.table_mid()

            if decoded_align.vectors:
                for r, atom in enumerate(repS_atoms):
                    row_vals = [vec[r] for vec in decoded_align.vectors]
                    cells_text = [atom_text(atom)] + [fmtf(v, 7) for v in row_vals]
                    cells_tex  = [f"${atom_tex(atom)}$"] + [tex_num(v) for v in row_vals]
                    out.table_row(cells_text, cells_tex)

            out.end_table()

            out.blank()
            out.line("  Ground truth  (direct evaluation of g·repS on event):")
            out._t(r"\textbf{Ground truth} (direct evaluation of $g \cdot \text{repS}$ on event):\\")
            out.begin_table(col_spec, caption="Ground truth alignment table")
            out.table_row(col_hdrs_text, col_hdrs_tex)
            out.table_mid()
            for r, atom in enumerate(repS_atoms):
                gt_row = [vec[r] for vec in gt_vectors]
                cells_text = [atom_text(atom)] + [fmtf(v, 7) for v in gt_row]
                cells_tex  = [f"${atom_tex(atom)}$"] + [tex_num(v) for v in gt_row]
                out.table_row(cells_text, cells_tex)
            out.end_table()

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
            status_tex  = (r"\textcolor{green!60!black}{\textbf{MATCH}}"
                           r" (multisets equal, atol$=10^{-8}$)") if ok \
                     else (r"\textcolor{red}{\textbf{MISMATCH}}"
                           r" — decoded multiset differs from ground truth")
            out.line(f"  {status_text}", tex=status_tex)

            out.blank()
            out.line("  Max |error| per repS row (decoded vs ground truth):")
            out._t(r"\textbf{Max $|$error$|$ per repS row:}\\")
            for r, atom in enumerate(repS_atoms):
                dec_row = sorted(vec[r] for vec in decoded_align.vectors)
                gt_row  = sorted(vec[r] for vec in gt_vectors)
                max_err = max(abs(d - g) for d, g in zip(dec_row, gt_row))
                out.kv(f"    {atom_text(atom)}",
                       f"max_err = {max_err:.2e}",
                       tex=f"${atom_tex(atom)}$: max err $= {max_err:.2e}$\\\\")

        # ── Summary ──────────────────────────────────────────────────────
        out.section("Summary")
        n_reps   = len(repS_atoms)
        n_cols   = ctx.the_group.order()
        running_total = 0
        if orbit_enc is not None:
            d1 = orbit_enc.output_dim
            running_total += d1
            out.kv("  Phase 1 encoded length", str(d1),
                   tex=rf"Phase 1 encoded length: ${d1}$ reals\\")
        if phase2_enc is not None:
            d2 = phase2_enc.output_dim
            running_total += d2
            out.kv("  Phase 2 encoded length", str(d2),
                   tex=rf"Phase 2 encoded length: ${d2}$ reals\\")
        if phase3_enc is not None:
            d3 = phase3_enc.output_dim
            running_total += d3
            out.kv("  Phase 3 encoded length", str(d3),
                   tex=rf"Phase 3 encoded length: ${d3}$ reals\\")
        out.kv("  Total encoded length",   str(running_total),
               tex=rf"Total encoded length: ${running_total}$ reals\\")
        out.kv("  repS size",              str(n_reps),
               tex=rf"$|\text{{repS}}| = {n_reps}$ atoms\\")
        out.kv("  |G|",                    str(n_cols),
               tex=rf"$|G| = {n_cols}$ group elements\\")
        out.line(f"\n  TeX output written to: {tex_path}")
        out.line("  Compile with:  pdflatex " + os.path.basename(tex_path))


# ─────────────────────────────────────────────────────────────────────────────
# Phase-3 narrative
# ─────────────────────────────────────────────────────────────────────────────


def _emit_phase3_section(out, plan, ctx, event: dict, phase3_enc, phase3_vals) -> None:
    """Render the Phase 3 narrative: rep atoms, alignment table T, encoder
    descriptor, and encoded-output summary.

    Phase 3 encoders consume the full rep alignment table T(event) as a
    multiset of |G| rows in R^{|rep|}.  Unlike Phase 1 and Phase 2 they
    don't have per-orbit / per-pair structure to traverse — there's a
    single composite output vector — so the narrative is shorter.

    Layout:
      * Section header
      * The rep atom list (one row per Phase 3 column)
      * The alignment table T (transposed: rows=atoms, cols=group elements)
        — same shape as the alignment-decoder ground truth elsewhere in the
        document, but emitted here for the Phase-3-only case where there is
        no alignment-decoder section.
      * The Phase 3 segment descriptors (kind, method, length per sub-encoder)
      * A small preview of the encoded output (first/last few values, range).
    """
    # Lazy import to avoid cycles at module load time.
    from symcoder.encoders.phase3_common import (
        rep_atoms_for, build_alignment_table,
    )

    rep_atoms = rep_atoms_for(plan)
    the_group = ctx.the_group
    n_cols    = the_group.order()
    n_rep     = len(rep_atoms)

    out.section(
        f"Phase 3 Encoding  (one-shot whole-table; rep size {n_rep}, |G|={n_cols})"
    )

    # ── rep atoms ────────────────────────────────────────────────────────
    out.line(
        f"  rep atoms ({n_rep}, the unpruned set used by Phase 3 — every "
        f"distinct atom, NOT just one canonical representative per orbit):",
        tex=(rf"\textbf{{rep atoms ({n_rep}):}} the unpruned set used by "
             rf"Phase 3 — every distinct atom across all G-orbits, not just "
             rf"one canonical representative per orbit.\\"),
    )
    for j, atom in enumerate(rep_atoms):
        out.kv(f"    [{j:>2d}]", atom_text(atom),
               tex=f"$[{j}]\\ \\ {atom_tex(atom)}$\\\\")

    # ── alignment table T ────────────────────────────────────────────────
    out.blank()
    out.line(
        "  Alignment table T  (rows = rep atoms, columns = group elements):",
        tex=(r"\textbf{Alignment table $T$} "
             r"(rows $=$ rep atoms, columns $= g \in G$):\\"),
    )

    T = build_alignment_table(plan, event, the_group, rep_atoms)  # (n_cols, n_rep)
    # Transposed display: rows = atom, columns = group element.
    col_spec     = "l" + "r" * n_cols
    col_hdrs_t   = ["atom"] + [f"col {c}" for c in range(n_cols)]
    col_hdrs_x   = ["Atom"] + [f"$g_{{{c}}}$" for c in range(n_cols)]
    out.begin_table(col_spec, caption="Phase 3 input alignment table T")
    out.table_row(col_hdrs_t, col_hdrs_x)
    out.table_mid()
    for j, atom in enumerate(rep_atoms):
        row_vals = T[:, j].tolist()
        cells_t  = [atom_text(atom)] + [fmtf(v, 7) for v in row_vals]
        cells_x  = [f"${atom_tex(atom)}$"] + [tex_num(v) for v in row_vals]
        out.table_row(cells_t, cells_x)
    out.end_table()

    # ── Segment descriptors (one per Phase 3 sub-encoder) ────────────────
    out.blank()
    tree = phase3_enc.describe(start_offset=0)
    out.line(
        f"  Phase 3 segments ({len(tree.segments)}):",
        tex=rf"\textbf{{Phase 3 segments ({len(tree.segments)}):}}\\",
    )
    for seg in tree.segments:
        out.kv(
            f"    [{seg.start:>5d}:{seg.stop:<5d}]  kind={seg.kind}  "
            f"method={seg.method_name}  len={seg.length}",
            "",
            tex=(rf"$[{seg.start}{{:}}{seg.stop}]$  "
                 rf"kind \texttt{{{tex_escape(seg.kind)}}}  "
                 rf"method \texttt{{{tex_escape(seg.method_name or '-')}}}  "
                 rf"length ${seg.length}$"
                 rf"\\"),
        )

    # ── Encoded-output preview ───────────────────────────────────────────
    import numpy as _np
    vals = _np.asarray(phase3_vals)
    out.blank()
    out.line(
        f"  Encoded output: {vals.size} reals, "
        f"range [{vals.min():+.4e}, {vals.max():+.4e}], "
        f"max|·| = {float(_np.max(_np.abs(vals))):+.4e}.",
        tex=(rf"\textbf{{Encoded output:}} ${vals.size}$ reals, "
             rf"range $[{vals.min():+.4e},\ {vals.max():+.4e}]$, "
             rf"max$|\cdot| = {float(_np.max(_np.abs(vals))):+.4e}$.\\"),
    )

    # Show up to first 8 and last 8 entries for context.
    head_n = min(8, vals.size)
    tail_n = min(8, vals.size - head_n)
    out.line("  First values:",
             tex=r"\textit{First values:}\\")
    out.kv("    ",
           "  ".join(fmtf(v) for v in vals[:head_n]),
           tex=r",\ ".join(tex_num(v) for v in vals[:head_n]))
    if tail_n > 0:
        out.line(f"  Last {tail_n} values:",
                 tex=rf"\textit{{Last {tail_n} values:}}\\")
        out.kv("    ",
               "  ".join(fmtf(v) for v in vals[-tail_n:]),
               tex=r",\ ".join(tex_num(v) for v in vals[-tail_n:]))


# ─────────────────────────────────────────────────────────────────────────────
# Phase-2 decoded output rendering
# ─────────────────────────────────────────────────────────────────────────────


def _emit_decoded_pair_multiset(out: DualOut, decoded) -> None:
    """Render a Phase-2-decoded PairFlavour result as TWO independent multisets.

    A previous version of this code listed each ``decoded.atom_pairs`` entry
    next to the corresponding ``decoded.pairs`` entry with an arrow, which
    suggested a 1-to-1 correspondence — there is none.  ``atom_pairs`` is
    an *orbit* of atom-pairs (algebraic / symbolic); ``pairs`` is an
    *independent orbit* of evaluated real pairs.  Both have the same length,
    but the ordering within each is arbitrary and they should not be zipped.
    """
    n_atoms = len(decoded.atom_pairs)
    n_pairs = len(decoded.pairs)

    if n_atoms == 0 and n_pairs == 0:
        return

    if n_atoms > 0:
        out.line("      Atom-pair orbit (multiset; ordering arbitrary):",
                 tex=r"\quad Atom-pair orbit (multiset; ordering arbitrary):\\")
        shown = decoded.atom_pairs[:8]
        for ap in shown:
            u_atom, v_atom = ap
            out.line(f"        ({atom_text(u_atom)}, {atom_text(v_atom)})",
                     tex=(f"$({atom_tex(u_atom)},\\ {atom_tex(v_atom)})$\\\\"))
        if n_atoms > 8:
            out.line(f"        ... ({n_atoms - 8} more atom-pairs)",
                     tex=rf"\quad \ldots ({n_atoms - 8} more atom-pairs)\\")

    if n_pairs > 0:
        out.line("      Eval-pair multiset (u_val, v_val):",
                 tex=r"\quad Eval-pair multiset $(u, v)$:\\")
        shown = decoded.pairs[:8]
        for pair in shown:
            u_val, v_val = pair
            out.line(f"        ({fmtf(u_val)}, {fmtf(v_val)})",
                     tex=(f"$({tex_num(u_val)},\\ {tex_num(v_val)})$\\\\"))
        if n_pairs > 8:
            out.line(f"        ... ({n_pairs - 8} more eval-pairs)",
                     tex=rf"\quad \ldots ({n_pairs - 8} more eval-pairs)\\")
