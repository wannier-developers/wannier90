#=======================
# NAG 
#======================= 

F90 = nagfor
MPIF90 = nagfor
FCOPTS = -O0 -DNAG

#=======================
# ATLAS Blas and LAPACK
#=======================
LIBDIR = /home/farm2/local/nag24/lib/NAG_Fortran/
LIBS = -L$(LIBDIR)  -llapack -lf77blas -lcblas -latlas

#=======================
# NETLIB LAPACK and BLAS
#=======================

#LIBDIR = /usr/local/lib
#LIBS = -L$(LIBDIR)  -llapack -lblas
LIBS = -llapack -lblas

# ACML for blas and lapack
#LIBS = -lacml

