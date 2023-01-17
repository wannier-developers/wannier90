include ../../make.inc

include ../../src/Makefile.header

ifeq ($(COMMS),mpi08)
DEFS = -DMPI -DMPI08
COMPILER = $(MPIF90)
else ifeq ($(COMMS),mpih)
DEFS = -DMPI -DMPIH
COMPILER = $(MPIF90)
else ifeq ($(COMMS),mpi90)
DEFS = -DMPI -DMPI90
COMPILER = $(MPIF90)
else ifeq ($(COMMS),mpi)
# default to f90 style "use mpi"
DEFS = -DMPI -DMPI90
COMPILER = $(MPIF90)
else
DEFS =
COMPILER = $(F90)
endif

a.out: demo.f90 $(LIBRARYV2)
	$(COMPILER) $(FCOPTS) $(DEFS) -I/home/jerome/wannier90/src/obj/ $< $(LDOPTS) $(LIBRARYV2) $(LIBS) -o $@
