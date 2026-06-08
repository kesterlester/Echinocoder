PDFLATEX := pdflatex

.PHONY: all clean

all: algorithm_description.pdf

# SOURCE_DATE_EPOCH=0 pins all embedded timestamps to the Unix epoch so the
# PDF is identical on repeated runs (no stored date/time metadata).
# Two pdflatex passes are needed: the first writes the .toc and .aux files;
# the second reads them back to resolve the table of contents and \ref links.
algorithm_description.pdf: algorithm_description.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode algorithm_description
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode algorithm_description

clean:
	rm -f algorithm_description.aux algorithm_description.log algorithm_description.out algorithm_description.toc algorithm_description.pdf
