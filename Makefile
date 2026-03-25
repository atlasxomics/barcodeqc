.PHONY: test test-e2e

PYTHON ?= $(if $(wildcard .venv/bin/python),./.venv/bin/python,python3)

test:
	$(PYTHON) -m pytest -m 'not e2e'

test-e2e:
	BARCODEQC_RUN_E2E=1 $(PYTHON) -m pytest tests/test_qc_e2e.py -q
