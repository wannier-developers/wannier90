
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=mpif90
FCOPTS=-O3 -assume byterecl -g -traceback -I/mnt/buildbot/binary/openmpi313_intel17/include/
LDOPTS=-O3

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR =  /mnt/buildbot/binary/intel17/compilers_and_libraries/linux/mkl/lib/intel64
LIBS   = -L/mnt/buildbot/binary/openmpi313_intel17/lib/ -lmpi_mpifh -lmpi  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

