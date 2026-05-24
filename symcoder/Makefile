PYTHON   := ../venv/bin/python
PDFLATEX := pdflatex

.PHONY: regenerate_pdfs

regenerate_pdfs: demo_roundtrip.pdf

demo_roundtrip.tex: demo_roundtrip.py
	$(PYTHON) $<

# SOURCE_DATE_EPOCH=0 pins all embedded timestamps to the Unix epoch so the
# PDF is identical on repeated runs (no stored date/time metadata).
demo_roundtrip.pdf: demo_roundtrip.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode $<
