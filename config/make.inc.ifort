
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=mpiifort
FCOPTS=-O2
LDOPTS=-O2

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = /opt/intel/mkl/lib/intel64
LIBS   =  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

#=======================
# ATLAS Blas and LAPACK
#=======================
#LIBDIR = /usr/local/lib
#LIBS = -L$(LIBDIR)  -llapack -lf77blas -lcblas -latlas
