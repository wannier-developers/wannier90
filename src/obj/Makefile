# -*- MakeFile -*-
#  Should be no need to change below this line
#

include ../../make.inc

# Contains definition of OBJS, OBJS_POST, LIBRARY, DYNLIBRARY, ...
include ../../Makefile.header

# Define the path where to actually look for files (the source .f90 files)
VPATH=..

OBJS  = constants.o io.o utility.o types.o hamiltonian.o overlap.o kmesh.o disentangle.o \
		ws_distance.o wannierise.o plot.o transport.o sitesym.o comms.o wannier90_types.o \
		wannier90_readwrite.o readwrite.o error.o error_base.o library_interface.o \
		library_extra.o
ifdef WANNIER90_WITH_C
OBJS  += c_interface.o wannier90.o
endif

OBJS_POST  = ws_distance.o types.o kmesh.o io.o comms.o utility.o get_oper.o constants.o \
 		postw90_common.o wan_ham.o spin.o dos.o berry.o gyrotropic.o kpath.o kslice.o \
		boltzwann.o geninterp.o postw90_types.o postw90_readwrite.o readwrite.o \
		error.o error_base.o tetrahedron.o

wannier: ../../wannier90.x
w90chk2chk: ../../w90chk2chk.x
w90spn2spn: ../../w90spn2spn.x
post: mpi_test ../../postw90.x

../../wannier90.x: wannier_prog.F90 ../../$(STATICLIBRARY) $(STATICLIBRARYEXTRA)
	$(COMPILER) $(POSTOPTS) $(FCOPTS) $< $(LDOPTS) ../../$(STATICLIBRARY) $(STATICLIBRARYEXTRA) $(LIBS) -o $@

../../w90chk2chk.x: w90chk2chk.F90 $(OBJS)
	$(COMPILER) $(POSTOPTS) $(FCOPTS) $< $(LDOPTS) $(OBJS) $(LIBS) -o $@

../../w90spn2spn.x: w90spn2spn.F90 $(OBJS)
	$(COMPILER) $(POSTOPTS) $(FCOPTS) $< $(LDOPTS) $(OBJS) $(LIBS) -o $@

../../postw90.x: $(POSTDIR)postw90.F90 $(OBJS_POST)
	$(COMPILER) $(POSTOPTS) $(FCOPTS) $< $(LDOPTS) $(OBJS_POST) $(LIBS) -o $@

../../$(STATICLIBRARY): $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

mpi_test :
ifeq ($(COMMS),mpi)

ifndef MPIF90

	$(info  ) \
	$(info ********************* STOP **********************) \
	$(info *  COMMS=mpi but MPIF90 is not set in make.inc  *) \
	$(info *  Set MPIF90 to your parallel fortran compiler *) \
	$(info *  On many systems MPIF90=mpif90 will work      *) \
	$(info ********************* STOP **********************) \
	$(info  ) \
	$(error MPIF90 is not set)

endif

endif


staticlib: ../../$(STATICLIBRARY)

dynlib: ../../$(DYNLIBRARY) ../../$(DYNLIBRARYEXTRA)

../../$(DYNLIBRARY): $(OBJS)
	$(COMPILER) $(SHAREDLIBFLAGS) $(POSTOPTS) $(LDOPTS) $(OBJS) -o $@ $(LIBS)

../../$(DYNLIBRARYEXTRA): $(OBJSEXTRA)
	$(COMPILER) $(SHAREDLIBFLAGS) $(POSTOPTS) $(LDOPTS) $(OBJSEXTRA) -o $@ -L.. -l$(DYNLIBBASE)



clean:
	rm -f *.o *.mod *.MOD *.obj

library_interface.o : library_interface.F90 \
	comms.o constants.o error_base.o io.o kmesh.o overlap.o plot.o transport.o \
	types.o wannier90_readwrite.o wannier90_types.o wannierise.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

library_extra.o : library_extra.F90 io.o library_interface.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

c_interface.o : c_interface.F90 \
	library_interface.o comms.o constants.o error_base.o io.o kmesh.o overlap.o plot.o transport.o \
	types.o wannier90_readwrite.o wannier90_types.o wannierise.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $< --std=f2018

wannier90.o: wannier90.c
	 $(CC) $(CFLAGS) -c $<

constants.o: constants.F90
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

io.o: io.F90 types.o constants.o error.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

utility.o: utility.F90 constants.o io.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

wannier90_types.o: wannier90_types.F90 constants.o io.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

wannier90_readwrite.o: wannier90_readwrite.F90 constants.o io.o types.o wannier90_types.o readwrite.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

readwrite.o: readwrite.F90 types.o constants.o io.o utility.o comms.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

types.o: types.F90 constants.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

hamiltonian.o: hamiltonian.F90 ws_distance.o constants.o io.o utility.o types.o wannier90_types.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

overlap.o: overlap.F90 constants.o io.o utility.o types.o sitesym.o wannier90_types.o comms.o error.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

kmesh.o: kmesh.F90 constants.o io.o utility.o types.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

error_base.o: error_base.F90
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

error.o: error.F90 error_base.o comms.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

disentangle.o: disentangle.F90 constants.o io.o types.o sitesym.o comms.o wannier90_types.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

wannierise.o: wannierise.F90 hamiltonian.o constants.o io.o utility.o types.o sitesym.o comms.o wannier90_types.o wannier90_readwrite.o error.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

plot.o: plot.F90 constants.o io.o utility.o types.o hamiltonian.o ws_distance.o wannier90_types.o comms.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

ws_distance.o: ws_distance.F90 constants.o io.o types.o utility.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

transport.o: transport.F90 constants.o io.o types.o hamiltonian.o wannier90_types.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

sitesym.o : sitesym.F90 utility.o types.o constants.o wannier90_types.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

comms.o: comms.F90 constants.o error_base.o ../../make.inc
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

boltzwann.o: $(POSTDIR)boltzwann.F90 constants.o types.o io.o utility.o postw90_common.o get_oper.o wan_ham.o comms.o spin.o dos.o postw90_types.o ws_distance.o postw90_readwrite.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

geninterp.o: $(POSTDIR)geninterp.F90 constants.o types.o io.o get_oper.o postw90_common.o comms.o utility.o wan_ham.o postw90_types.o ws_distance.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

get_oper.o: $(POSTDIR)get_oper.F90 types.o constants.o comms.o postw90_common.o io.o utility.o postw90_types.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

spin.o: $(POSTDIR)spin.F90 comms.o types.o constants.o utility.o postw90_common.o get_oper.o io.o wan_ham.o postw90_types.o ws_distance.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

dos.o: $(POSTDIR)dos.F90 comms.o types.o postw90_types.o constants.o utility.o io.o wan_ham.o postw90_common.o get_oper.o spin.o ws_distance.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

kpath.o: $(POSTDIR)kpath.F90 comms.o types.o constants.o io.o spin.o berry.o postw90_types.o ws_distance.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

kslice.o: $(POSTDIR)kslice.F90 comms.o types.o constants.o io.o spin.o berry.o postw90_types.o ws_distance.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

berry.o: $(POSTDIR)berry.F90 comms.o types.o constants.o utility.o postw90_common.o get_oper.o io.o spin.o wan_ham.o postw90_types.o ws_distance.o tetrahedron.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

gyrotropic.o: $(POSTDIR)gyrotropic.F90 comms.o types.o constants.o utility.o postw90_common.o get_oper.o io.o spin.o wan_ham.o berry.o postw90_types.o ws_distance.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

wan_ham.o: $(POSTDIR)wan_ham.F90 types.o postw90_types.o constants.o utility.o postw90_common.o get_oper.o ws_distance.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

postw90_common.o: $(POSTDIR)postw90_common.F90 postw90_types.o ws_distance.o comms.o types.o utility.o constants.o io.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

postw90_types.o: $(POSTDIR)postw90_types.F90 comms.o types.o utility.o constants.o io.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

postw90_readwrite.o: $(POSTDIR)postw90_readwrite.F90 comms.o postw90_types.o types.o utility.o constants.o io.o error.o wannier90_readwrite.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

tetrahedron.o: $(POSTDIR)tetrahedron.F90 constants.o utility.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $<

.PHONY: wannier staticlib dynlib post clean mpi_test w90spn2spn
