#===================
# gfortran
#===================
F90 = gfortran

#Add next two lines for parallel postw90. Choose the correct name of the mpi f90 wrapper
#COMMS  = mpi
MPIF90 = mpif90

FCOPTS = -O3 -g -x f95-cpp-input -I/home/farmer4/local/openmpi110_gcc485/include
LDOPTS = -g -pthread

#Next two lines are good for debugging
#FCOPTS = -Wall -Waliasing -Wsurprising -Wline-truncation -Wno-tabs -Wno-uninitialized -Wno-unused-dummy-argument -Wno-unused -Wno-character-truncation -O1 -g -fbacktrace
#LDOPTS = -Wall -Waliasing -Wsurprising -Wline-truncation -Wno-tabs -Wno-uninitialized -Wno-unused-dummy-argument -Wno-unused -Wno-character-truncation -O1 -g -fbacktrace


#=======================
# ATLAS Blas and LAPACK
#=======================
LIBDIR = /home/farmer4/local/OpenBLAS-0.2.18_gcc485/lib/
LIBS = -L$(LIBDIR)  -llapack -lblas


#=======================
# NETLIB LAPACK and BLAS
#=======================
#LIBDIR = /rscratch/jry20/LIB/gfortran
#LIBS = -L$(LIBDIR)  -llapack -lblas

