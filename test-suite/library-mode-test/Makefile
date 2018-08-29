include ../../make.inc

include ../../src/Makefile.header

ifdef COMMS
COMMS := $(strip $(COMMS))
else
COMMS = serial
endif
# Only works in serial for now
# Otherwise, in the future, put back the logic here from 
# Makefile.2
#COMMS = serial

ifeq ($(COMMS),mpi)
TEMP1 = -DMPI
TEMP2 = $(MPIF90)
else
TEMP1 =
TEMP2 = $(F90)
endif

librarytest: POSTOPTS = $(TEMP1) 
librarytest: COMPILER = $(TEMP2)
librarytest: test_library_serial.x

test_library_serial.x: test_library.F90 $(LIBRARY)
	$(COMPILER) $< $(POSTOPTS) $(LDOPTS) $(LIBS) $(LIBRARY) -o $@
#	$(F90) -fstrict-aliasing  -fno-omit-frame-pointer -fno-realloc-lhs -fcheck=bounds,do,recursion,pointer -ffree-form -Wall -Waliasing -Wsurprising -Wline-truncation -Wno-tabs -Wno-uninitialized -Wno-unused-dummy-argument -Wno-unused -Wno-character-truncation -O1 -g -fbacktrace -fPIC $< ../../libwannier.a -lblas -llapack -o $@

clean:
	rm -f test_library_serial.x *~

.PHONY: clean librarytest
