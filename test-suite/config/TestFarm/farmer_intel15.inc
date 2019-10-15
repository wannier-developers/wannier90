
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=mpif90
FCOPTS=-O3 -assume byterecl -g -traceback -I/home/farmer/local/openmpi188_intel15/include/
LDOPTS=-O3

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR =  /home/farmer/local/intel15/composer_xe_2015.3.187/mkl/lib/intel64/
LIBS   = -L/home/farmer/local/openmpi188_intel15/lib/ -lmpi_mpifh -lmpi  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

