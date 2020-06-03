#=====
#pgf90 
#=====

F90 = pgf90

#Add next two lines for parallel postw90. Choose the correct name of the mpi f90 wrapper
COMMS = mpi
MPIF90 =  mpif90

FCOPTS = -O2 -g f95-cpp-input -I/mnt/buildbot/binary/mvapich23_pgi185/include
LDOPTS = -O2

# Debug options
#FCOPTS  = -g -Mbounds -Mchkptr -Minform=warn -Mstandard 
#LDOPTS = -g -Mbounds -Mchkptr -Minform=warn -Mstandard
LIBDIR =  /mnt/buildbot/binary/mvapich23_pgi185/lib
LIBS   = -llapack -lblas  -L$(LIBDIR) -lmpi

