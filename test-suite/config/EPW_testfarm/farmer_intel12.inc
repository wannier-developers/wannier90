
#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
F90 = ifort
#COMMS=mpi
#MPIF90=ifort
FCOPTS=-O2 -assume byterecl -g -traceback 
LDOPTS=-O2

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR =  /home/farmer/local/intel11/composer_xe_2011_sp1.8.273/mkl/lib/intel64
LIBS   = -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread
