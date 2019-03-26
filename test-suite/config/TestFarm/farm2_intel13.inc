#=====================================================
# For Linux with intel
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=ifort
FCOPTS=-O2 -assume byterecl -g -traceback -I/home/farm2/local/mpich314_intel13/include/
LDOPTS=-O2

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = /home/farm2/local/intel13/composer_xe_2013.1.117/mkl/lib/intel64
LIBS   = -L/home/farm2/local/mpich314_intel13/lib/ -lmpifort  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

