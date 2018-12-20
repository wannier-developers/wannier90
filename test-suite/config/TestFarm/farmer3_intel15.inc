#=====================================================
# For Linux with intel
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=mpif90
FCOPTS=-O0 -assume byterecl -g -traceback -I/home/farmer3/local/mvapich222a_intel15/include/
LDOPTS=-O0

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = /home/farmer3/local/intel15/composer_xe_2015.3.187/mkl/lib/intel64
LIBS   = -L/home/farmer3/local/mvapich222a_intel15/lib64/ -lmpi -lmpifort  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

