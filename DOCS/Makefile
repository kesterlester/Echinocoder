PDFLATEX := pdflatex

.PHONY: regenerate_pdfs clean

# All locally-buildable PDFs in this directory.  Add new entries here when
# new .tex sources are introduced.  (Third-party PDFs that are checked in
# but not built from .tex sources — e.g. Davis-Boletin-... — are NOT
# listed and are left untouched by this Makefile.)
regenerate_pdfs: encoder.pdf parity_blindness_and_repS_concern.pdf phase3_vandermonde.pdf simplicial_complex_embedders.pdf

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
