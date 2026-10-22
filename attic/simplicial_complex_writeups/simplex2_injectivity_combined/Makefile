# Makefile for the combined writeup (master.pdf).
#
# master.tex \input's the two original bodies (part1_embedder, part2_injectivity),
# the new bridge, and the embedder appendix.  The worked-example tables come in
# through the 'generated' symlink -> ../generated, so 'make tables' in DOCS/
# keeps them in lock-step with EncDec.py.  This Makefile does NOT regenerate
# those tables (run 'make tables' in DOCS/ for that); it only builds the PDF.
#
# Three pdflatex passes: .aux/.toc on the first, \tableofcontents and the
# forward/backward \ref between the embedder, bridge and proof on the second,
# and hyperref/page-ref settling on the third.  SOURCE_DATE_EPOCH=0 pins
# embedded timestamps for byte-reproducible output (matches DOCS/Makefile).

PDFLATEX := pdflatex
PARTS := master.tex part1_embedder.tex bridge.tex part2_injectivity.tex appendix_code.tex

.PHONY: all clean

all: master.pdf

master.pdf: $(PARTS) generated/sx2_stepD_main.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode master.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode master.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode master.tex

clean:
	rm -f *.log *.aux *.out *.toc
