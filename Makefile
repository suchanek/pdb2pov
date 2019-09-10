#
#
# $Header: EGS:pdb2povf/RCS/Makefile,v 1.3 1993/11/26 10:24:17 eric Exp eric $
#
# Makefile for thee CPK rendering package
#
# See the docs for instructions on building this package.
#
# Copyright 1992 Eric G. Suchanek, Ph.D.
#
#
#
# $Log: Makefile,v $
# Revision 1.3  1993/11/26  10:24:17  eric
# added some new depends
#
# Revision 1.2  1993/11/26  00:26:03  eric
# additional depends for multiple cpu types
#
# Revision 1.1  1993/11/13  11:02:08  eric
# Initial revision
#
#

.precious: pdb2pov.c

DEBUG = -g
DELETE=rm

#
# The cflags for compiling pdb2pov under UNIX hosts
#


CC = cc
OPT = -O
CFLAGS = $(OPT) -c 


#
# The flags for linking
#

LDFLAGS	 = 
DBGFLAGS = 
LLIBS	= 
STARTUP	=

VERS=119

#
# The default make rule set
#

.SUFFIXES:
.SUFFIXES: .c .o .o30 .q .d .i .h .pr .o40 .guide .cp .cps .fn .fns .vr .vrs .dvi .texinfo

.c.o: 
	$(CC) $(CFLAGS) -o $*.o $*.c


.texinfo.dvi:
	virtex20 &plain $*.texinfo

.texinfo.guide:
	makeinfo --amiga $*

.cp.cps:
	texindex $*.cp

.vr.vrs:
	texindex $*.vr

.fn.fns:
	texindex $*.fn

P = pdb2povf/
SRC = pdb2pov.c m68040.c asyncio.c util.c pdb2light.h pdb2pov_errors.h pdb2pov_usage.h
HEADERS = pdb2light.h pdb2pov_errors.h pdb2pov_usage.h
SAMPLES = crambin.pdb crambin.pov
INCLUDE_T = atoms2.inc atoms_cpk.inc atoms_vdw.inc atoms_covalent.inc atoms_glass2.inc

PROG_DOC = pdb2pov.dvi pdb2pov.dvi.info pdb2pov.guide pdb2pov.guide.info 

FILES = $(SRC) makefile Makefile.unix $(INCLUDE_T) $(PROG_DOC)

OBJ = pdb2pov.o util.o

PROGS = pdb2pov
DOC = pdb2pov.dvi Pdb2POV.dvi.info pdb2pov.guide pdb2pov.guide.info About.Sids About.Sids.info 

BABEL = 
P = 

SRC_ARCHIVE = pdb2pov_src_$(VERS).tar
ARCHIVE = pdb2pov_$(VERS).tar

pdb2pov: $(OBJ)
	$(CC) -o pdb2pov $(OBJ) $(LLIBS)


pdb2pov.o:   $(HEADERS)

all: $(PROGS) $(PROG_DOC)


$(SRC_ARCHIVE): $(FILES)
	tar cvof $(SRC_ARCHIVE) $(FILES)

src_arc: $(SRC_ARCHIVE)

clean:
	@$(DELETE) $(OBJ)

realclean:
	@$(DELETE) $(PROGS)

IND2 = pdb2pov.cps pdb2pov.vrs pdb2pov.fns
IND = 

#pdb2pov.dvi:  pdb2pov.texinfo $(IND2)
#pdb2pov.guide:  pdb2pov.texinfo $(IND2)

doc: $(DOC)

everything: cpk arc src

