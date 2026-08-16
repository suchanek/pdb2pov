#
# Makefile for the Python port of pdb2pov.
#
# Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.
# Subject to the GNU License.
#
# There is nothing to build -- the package is pure Python with no
# dependencies.  These targets are for testing and for the two conversions
# the C tree's `make check` runs, so the two can be compared by eye as well
# as by the differential tests.
#
# Targets:
#   make test       run the test suite
#   make check      run the suite, then convert the bundled crambin
#   make lint       compile every module, and check the public API imports
#   make clean      remove caches and check artefacts
#

PYTHON  ?= python3
PYTEST  ?= $(PYTHON) -m pytest
SRC      = src
RUN      = PYTHONPATH=$(SRC) $(PYTHON) -m pdb2pov

.PHONY: all test check lint clean

all: test

test:
	PYTHONPATH=$(SRC) $(PYTEST) -q

#
# The same conversions the C tree's `make check` performs, so the two can be
# diffed by hand.  The differential tests do this automatically; this target
# is for when you want to look at the files.
#
check: test
	@echo "--- converting the bundled crambin ---"
	$(RUN) ../1CRN check_py_vdw -v -p --quiet
	$(RUN) ../1CRN check_py_bs  -b -d 1.9 -p --quiet
	$(RUN) ../1CRN check_py_obj -b -d 1.9 -o --quiet
	$(RUN) ../elements check_py_elements -v -p --quiet
	@ls -l check_py_*.pov check_py_obj.inc
	@echo "--- every element must resolve to its own texture ---"
	@if grep -qE 'Atom_X\b' check_py_elements.pov; then \
	  echo "FAIL: an element in the table has no dedicated texture"; exit 1; \
	else \
	  echo "OK: every element resolved to its own texture"; \
	fi
	@#
	@# The differential tests already diff against the C program when it is
	@# present.  Say so either way, so a skip is not mistaken for a pass.
	@#
	@if [ -f ../pdb2pov.c ]; then \
	  echo "OK: the C sources are present, so the differential tests ran"; \
	else \
	  echo "NOTE: no C sources alongside; the differential tests were skipped"; \
	fi

lint:
	$(PYTHON) -m compileall -q $(SRC) tests
	PYTHONPATH=$(SRC) $(PYTHON) -c "import pdb2pov; \
	  [getattr(pdb2pov, n) for n in pdb2pov.__all__]; \
	  print('OK:', len(pdb2pov.__all__), 'public names')"

clean:
	rm -f check_py_*.pov check_py_*.inc
	rm -rf .pytest_cache
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
