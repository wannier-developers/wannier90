ROOT = ../..
include $(ROOT)/make.inc
# COMPILER (F90 or MPIF90), POSTOPTS and DYNLIBBASE (serial/MPI library name) per COMMS
include $(ROOT)/Makefile.header

# the fortran compiler drives the link (the library carries the fortran and
# MPI runtimes) and the C++ runtime is added explicitly; under MPI the demo
# itself calls MPI_Init and sets the communicator, so the C++ compile needs
# mpi.h and -DMPI (from POSTOPTS)
ifdef USEMPI
CXX = mpicxx
endif
CXXOPTS ?= -O2 -g -Wall

wannier_c.x: wannier.o c_interface_check
	$(COMPILER) $(LDOPTS) $< -o $@ -L$(ROOT) -l$(DYNLIBBASE) $(LIBS) -lstdc++ -Wl,-rpath,$(abspath $(ROOT))
wannier.o: wannier.cc $(ROOT)/src/wannier90.h
	$(CXX) $(CXXOPTS) $(POSTOPTS) -I$(ROOT)/src -c $<

# the C interface enters libwannier90 only when the library is built with
# WANNIER90_WITH_C set; test the library itself rather than the variable,
# since make.inc may have changed after the build
c_interface_check:
	@nm -g --defined-only $(wildcard $(ROOT)/lib$(DYNLIBBASE).a $(ROOT)/lib$(DYNLIBBASE).so*) 2>/dev/null \
	  | grep -qw w90_create || \
	  { echo "lib$(DYNLIBBASE) absent or built without the C interface: make libs WANNIER90_WITH_C=1"; false; }

clean:
	@rm -f wannier.o wannier_c.x
.PHONY: clean c_interface_check
