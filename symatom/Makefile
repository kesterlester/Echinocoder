ROOT    := $(abspath $(CURDIR)/..)
PYTHON  := $(ROOT)/venv/bin/python
PYTEST  := $(ROOT)/venv/bin/pytest

.PHONY: demo test update-demo-output help

help:
	@echo "Targets:"
	@echo "  demo                Run the demo and print to terminal"
	@echo "  test                Run the symatom test suite"
	@echo "  update-demo-output  Regenerate demo_output.txt"

demo:
	cd $(ROOT) && $(PYTHON) -m symatom.demo

test:
	cd $(ROOT) && $(PYTEST) symatom/ -v

update-demo-output:
	cd $(ROOT) && $(PYTHON) -m symatom.demo > symatom/demo_output.txt
	@echo "demo_output.txt updated"
