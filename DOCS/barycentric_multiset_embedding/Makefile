PDFLATEX := pdflatex
PYTHON   := ../../venv/bin/python

.PHONY: regenerate_pdfs clean tables

# All locally-buildable PDFs in this directory.  Add new entries here when
# new .tex sources are introduced.  (Third-party PDFs that are checked in
# but not built from .tex sources — e.g. Davis-Boletin-... — are NOT
# listed and are left untouched by this Makefile.)
regenerate_pdfs: main.pdf

# ─── Generic .tex → .pdf rule ─────────────────────────────────────────────
#
# Two pdflatex passes per file: the first writes the .aux / .toc files,
# the second resolves \tableofcontents, \ref/\label and any other forward
# references using them.  Both .tex files in this directory currently use
# section labels and at least one has a table of contents, so the second
# pass is needed for correctness even though it is sometimes redundant.
#
# SOURCE_DATE_EPOCH=0 pins all embedded PDF timestamps to the Unix epoch
# so the output is byte-identical across repeated runs.  This matches the
# convention used in symcoder/Makefile.

%.pdf: %.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode $<
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode $<

clean:
	rm -f *.log *.aux *.out *.toc

# ─── Change-marked ("track changes") diff PDF ─────────────────────────────
#
# main-CONEDIFF.pdf shows the current main.tex against a baseline revision,
# with removed text struck through in red and added text underlined in blue.
#
#   make diff                  build it
#   make main-CONEDIFF.pdf     the same thing, named explicitly
#   make diff DIFFBASE=main    compare against some other revision instead
#   make diffclean             remove the diff and its scratch files
#
# The baseline lives in git rather than in a file, so make cannot stat it and
# cannot tell on its own whether it has moved.  We therefore extract it into
# main-CONEDIFF-base.tex on every invocation --- a `git show`, which costs
# milliseconds --- but only *replace* that file when the extracted text
# actually differs.  Its timestamp thus advances only on a genuine change, so
# the expensive steps (latexdiff, then two pdflatex passes) are skipped by the
# ordinary prerequisite check whenever neither the baseline nor main.tex has
# moved since the last build.

DIFFBASE  := HEAD
LATEXDIFF := latexdiff
DIFFOPTS  := --append-safecmd="bvec,bary,R,Sn" \
             --config="PICTUREENV=(?:picture|DIFnomarkup|tikzpicture)"

.PHONY: diff diffclean FORCE

diff: main-CONEDIFF.pdf

FORCE:

# The baseline text, extracted from git.  FORCE makes the recipe run on every
# invocation; the compare-and-swap below leaves the timestamp alone unless the
# content genuinely changed.
main-CONEDIFF-base.tex: FORCE
	@git show $(DIFFBASE):./main.tex > $@.tmp || { rm -f $@.tmp; exit 1; }
	@if cmp -s $@.tmp $@; then \
	  rm -f $@.tmp; \
	else \
	  mv $@.tmp $@; \
	  echo "baseline main.tex re-extracted from $(DIFFBASE)"; \
	fi

# latexdiff marks up baseline against working copy.  --append-safecmd protects
# this document's own single-argument macros from being broken apart by the
# markup, and PICTUREENV keeps latexdiff out of the tikzpicture, which it would
# otherwise rewrite into something that does not compile.
main-CONEDIFF.tex: main-CONEDIFF-base.tex main.tex
	$(LATEXDIFF) $(DIFFOPTS) main-CONEDIFF-base.tex main.tex > $@

# main-CONEDIFF.pdf itself needs no rule of its own: the generic %.pdf: %.tex
# rule above builds it from main-CONEDIFF.tex, two passes and all.

diffclean:
	rm -f main-CONEDIFF.pdf main-CONEDIFF.tex main-CONEDIFF-base.tex \
	      main-CONEDIFF-base.tex.tmp main-CONEDIFF.log main-CONEDIFF.aux \
	      main-CONEDIFF.out main-CONEDIFF.toc
