#
# Makefile for pdb2pov -- Brookhaven PDB to POV-Ray scene conversion.
#
# Copyright (c) 1993-2026 Eric G. Suchanek, Ph.D.
# Subject to the GNU License.
#
# The 1993 Makefile carried a PORTFLAGS set that force-included a prototype
# header, disabled _FORTIFY_SOURCE and the stack protector, and pinned the
# compiler to -std=gnu89, all so the K&R sources would build on a 64-bit host.
# The sources are prototyped C17 now, so none of that is needed: the compiler
# can see every declaration, and the buffer overrun that FORTIFY tripped on is
# fixed rather than suppressed.
#
# Targets:
#   make            build pdb2pov
#   make check      build and run a conversion of the bundled 1CRN.pdb
#   make clean      remove objects
#   make distclean  remove objects, the binary and check artefacts
#

CC      ?= cc
CSTD    ?= -std=c17
OPT     ?= -O2
WARN    ?= -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes \
           -Wmissing-prototypes -Wpointer-arith -Wcast-qual -Wwrite-strings

CFLAGS  += $(CSTD) $(OPT) $(WARN)

# libm is not linked implicitly by every toolchain.  macOS gets it via
# libSystem, which is why the omission went unnoticed for thirty years; on
# GNU/Linux the link fails on sincos and hypot without it.
LDLIBS  += -lm

PROG    = pdb2pov
OBJ     = pdb2pov.o util.o
HEADERS = pdb2pov.h pdb2pov_errors.h pdb2pov_usage.h util.h

INCLUDES = atoms2.inc atoms_cpk.inc atoms_vdw.inc atoms_covalent.inc \
           atoms_glass2.inc

.PHONY: all check clean distclean

all: $(PROG)

$(PROG): $(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

pdb2pov.o: pdb2pov.c $(HEADERS)
util.o: util.c util.h

#
# A smoke test: trim the bundled crambin structure to coordinate records --
# the parser predates several PDB conventions and REMARK records will upset
# it -- then convert it both ways.
#
check: $(PROG)
	@grep -E "^(ATOM|HETATM|END)" 1CRN.pdb > check_crambin.pdb
	./$(PROG) check_crambin check_crambin_vdw -v -p
	./$(PROG) check_crambin check_crambin_bs -b -d 1.9 -p
	./$(PROG) check_crambin check_crambin_obj -b -d 1.9 -o
	./$(PROG) check_crambin check_crambin_leg -b -d 1.9 -p --legacy-elements
	./$(PROG) check_crambin check_crambin_chn -b -d 1.9 -p --chain A
	@echo "--- generated ---"
	@ls -l check_crambin_*.pov check_crambin_obj.inc
	@echo "--- crambin must be unaffected by the 2.1 parser changes ---"
	@grep -v 'Prepared by' check_crambin_bs.pov \
	   | sed 's/check_crambin_bs/N/g' > check_a.tmp
	@grep -v 'Prepared by' check_crambin_leg.pov \
	   | sed 's/check_crambin_leg/N/g' > check_b.tmp
	@if diff check_a.tmp check_b.tmp > /dev/null; then \
	  echo "OK: default and --legacy-elements agree on crambin"; \
	  rm -f check_a.tmp check_b.tmp; \
	else \
	  echo "FAIL: crambin differs between default and --legacy-elements"; \
	  rm -f check_a.tmp check_b.tmp; exit 1; \
	fi

clean:
	rm -f $(OBJ)

distclean: clean
	rm -f $(PROG) check_crambin.pdb check_crambin_*.pov check_crambin_*.inc check_*.tmp
