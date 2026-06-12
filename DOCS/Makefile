PDFLATEX := pdflatex
PYTHON   := ../venv/bin/python

.PHONY: regenerate_pdfs clean tables

# All locally-buildable PDFs in this directory.  Add new entries here when
# new .tex sources are introduced.  (Third-party PDFs that are checked in
# but not built from .tex sources — e.g. Davis-Boletin-... — are NOT
# listed and are left untouched by this Makefile.)
regenerate_pdfs: encoder.pdf parity_blindness_and_repS_concern.pdf phase3_vandermonde.pdf simplicial_complex_embedders_v1.pdf simplicial_complex_embedders_v2.pdf simplicial_complex_embedders_v3.pdf

# ─── Generated worked-example tables (v3) ─────────────────────────────────
#
# The Algorithm-2 worked-example tables in v3 are generated from EncDec.py
# so that the document stays in lock-step with the code.  To change the
# example multiset, edit M in generate_simplex_example_tables.py and re-run.
# sx2_stepD_main.tex stands in for all the sx2_*.tex fragments the script
# writes in one run.

generated/sx2_stepD_main.tex: generate_simplex_example_tables.py ../EncDec.py
	$(PYTHON) generate_simplex_example_tables.py

tables: generated/sx2_stepD_main.tex

simplicial_complex_embedders_v3.pdf: generated/sx2_stepD_main.tex

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
