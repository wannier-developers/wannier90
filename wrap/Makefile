
include ../make.inc

# Do this before Makefile.header so we can annotate dynlib with COMMS
ifdef COMMS
COMMS := $(strip $(COMMS))
PYMOD = wan90$(COMMS)
F90WRAPEXTRA = --f77exec=mpif77 --f90exec=mpif90
else
COMMS = serial
PYMOD = wan90
endif

# Contains definition of OBJS, OBJS_POST, LIBRARY, DYNLIBRARY
include ../src/Makefile.header

# location of modules from main code
BUILD = ../src/obj

all: wrap python

wrap: $(WRAPSRC)
	$(F90WRAP) -m $(PYMOD) -k kind_map $(WRAPSRC)

python: cabort.o
	f2py -c -m _$(PYMOD) $(F90WRAPEXTRA) -I$(BUILD) f90wrap_*.f90 cabort.o -L.. -l$(DYNLIBBASE)

# I expected this to only be needed for f2py-f90wrap, but seems to be wrong
cabort.o: cabort.c
	gcc -c -fPIC cabort.c
