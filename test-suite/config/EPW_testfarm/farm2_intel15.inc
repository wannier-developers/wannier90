#=====================================================
# For Linux with intel
#=====================================================
F90 = ifort
COMMS=mpi
MPIF90=ifort
FCOPTS=-O2 -assume byterecl -g -traceback -I/home/farm2/local/impi51_intel15/compilers_and_libraries_2016.0.079/linux/mpi/intel64/include/
LDOPTS=-O2

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================

LIBDIR = /home/farm2/local/intel15/composer_xe_2015.3.187/mkl/lib/intel64/
LIBS   = -L/home/farm2/local/impi51_intel15/compilers_and_libraries_2016.0.079/linux/mpi/intel64/lib/ -lmpi -lmpifort  -L$(LIBDIR) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

