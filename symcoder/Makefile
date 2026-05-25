PYTHON   := ../venv/bin/python
PDFLATEX := pdflatex

.PHONY: regenerate_pdfs

regenerate_pdfs: demo_roundtrip.pdf demo_shadow_event_parity_+1.pdf demo_shadow_event_parity_-1.pdf

# ─── Round-trip demo ──────────────────────────────────────────────────────

demo_roundtrip.tex: demo_roundtrip.py $(wildcard *.py) $(wildcard ../symatom/*.py) $(wildcard operations/*.py) $(wildcard encoders/*.py)
	$(PYTHON) $<

# SOURCE_DATE_EPOCH=0 pins all embedded timestamps to the Unix epoch so the
# PDF is identical on repeated runs (no stored date/time metadata).
demo_roundtrip.pdf: demo_roundtrip.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode $<

# ─── Parity-blindness demo  (one Python invocation = both .tex files) ─────
#
# demo_shadow_event.py runs the full encode/decode pipeline twice — once on
# a chosen "shadow event" E, once on its parity flip P(E) = {l: -v} — and
# writes the two narratives to two .tex files.  Put both files on the LHS
# of one rule so make knows the single recipe produces both.

demo_shadow_event_parity_+1.tex demo_shadow_event_parity_-1.tex: demo_shadow_event.py $(wildcard *.py) $(wildcard tests/_parity_shadow_events.py) $(wildcard ../symatom/*.py) $(wildcard operations/*.py) $(wildcard encoders/*.py)
	$(PYTHON) $<

demo_shadow_event_parity_+1.pdf: demo_shadow_event_parity_+1.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode $<

demo_shadow_event_parity_-1.pdf: demo_shadow_event_parity_-1.tex
	SOURCE_DATE_EPOCH=0 $(PDFLATEX) -interaction=nonstopmode $<
