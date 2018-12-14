#=====
#pgf90 
#=====

F90 = pgf90

#Add next two lines for parallel postw90. Choose the correct name of the mpi f90 wrapper
#COMMS = mpi
MPIF90 =  mpif90

FCOPTS = -O2 -fast -I/home/farm2/local/pgi_1507/linux86-64/15.7/mpi/mpich/include
LDOPTS = -O2

# Debug options
#FCOPTS  = -g -Mbounds -Mchkptr -Minform=warn -Mstandard 
#LDOPTS = -g -Mbounds -Mchkptr -Minform=warn -Mstandard
LIBDIR =  /home/farm2/local/pgi_1507/linux86-64/15.7/mpi/mpich/lib
LIBS   = -llapack -lblas  -L$(LIBDIR) -lmpi

