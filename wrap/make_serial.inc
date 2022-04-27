#===================
# gfortran
#===================
F90 = gfortran

#Add next two lines for parallel postw90. Choose the correct name of the mpi f90 wrapper
#COMMS  = mpi
#MPIF90 = mpgfortran #mpif90

FCOPTS = -g -O0 -fPIC -Wall -Wextra -fcheck=all -fbacktrace -finit-real=nan -finit-integer=0x7fffffff -DW90DEV
LDOPTS = $(FCOPTS)

MODFLAG = -J
#LIBS = -llapack -lblas

LIBS = -Wl,--start-group /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_gf_lp64.a /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_core.a /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -ldl
