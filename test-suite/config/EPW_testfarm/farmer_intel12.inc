
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=ifort
FCOPTS=-O0 -assume byterecl -g -traceback -I/home/farmer/local/openmpi188_intel12/include/
LDOPTS=-O0

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR =  /home/farmer/local/intel11/composer_xe_2011_sp1.8.273/mkl/lib/intel64
LIBS   = -L/home/farmer/local/openmpi188_intel12/lib/ -lmpi_mpifh -lmpi  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread
