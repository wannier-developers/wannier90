
include ../make.inc

include ../Makefile.header

ifdef COMMS
COMMS := $(strip $(COMMS))
PYMOD = wan90$(COMMS)
F90WRAPEXTRA = --f77exec=mpif77 --f90exec=mpif90
else
COMMS = serial
PYMOD = wan90
endif

# Define this if not already defined
F90WRAP ?= f90wrap
F2PY ?= f2py-f90wrap

WRAPSRC = ../src/constants.F90 ../src/types.F90 ../src/wannier90_readwrite.F90 \
        ../src/wannier90_types.F90 ../src/comms.F90 ../src/library_interface.F90 \
        ../src/library_extra.F90

# location of modules from main code
BUILD = $(realpath ../src/obj)

all: wrap python

wrap: $(WRAPSRC)
	$(F90WRAP) -m $(PYMOD) -k kind_map $(WRAPSRC)

python: cabort.o
	$(F2PY) -c -m _$(PYMOD) $(F90WRAPEXTRA) -I$(BUILD) f90wrap_*.f90 cabort.o -L$(realpath ..) -l$(DYNLIBBASE)

# I expected this to only be needed for f2py-f90wrap, but seems to be wrong
cabort.o: cabort.c
	gcc -c -fPIC cabort.c

.PHONY: all wrap python
