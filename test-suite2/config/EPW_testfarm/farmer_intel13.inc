#=====================================================
# For Linux with intel 
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=mpif90
FCOPTS=-O2 -assume byterecl -g -traceback -I/home/farmer/local/openmpi188_intel13/include/
LDOPTS=-O2

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = /home/farmer/local/intel13/composer_xe_2013.5.192/mkl/lib/intel64
LIBS   = -L/home/farmer/local/openmpi188_intel13/lib/ -lmpi_mpifh -lmpi  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

