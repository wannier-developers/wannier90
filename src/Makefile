#  Should be no need to change below this line
#

include ../../make.sys

OBJS  =  constants.o io.o utility.o parameters.o hamiltonian.o overlap.o \
	 kmesh.o disentangle.o wannierise.o plot.o transport.o

OBJS2  =  wannier_lib.o

OBJS_POST  = parameters.o kmesh.o io.o comms.o utility.o get_oper.o constants.o wanint_common.o wan_ham.o spin_wanint.o dos_wanint.o berry_wanint.o kpath_plot.o slice_plot.o

LIBRARY = ../../libwannier.a

POSTDIR = ../postw90/

ifdef COMMS
COMMS := $(strip $(COMMS))
else
COMMS = serial
#COMPILER = $(F90)
endif
#ifeq ($(COMMS),mpi)
#POSTOPTS = -DMPI
#COMPILER = $(MPIF90)
#else
#POSTOPTS =
#COMPILER = $(F90)
#endif

#comms_test_prog:
#POSTOPTS = 
#COMPILER = $(F90)

#comms_test_post:
#	$(eval 
ifeq ($(COMMS),mpi)
TEMP1 = -DMPI
TEMP2 = $(MPIF90)
else
TEMP1 =
TEMP2 = $(F90)
endif


#	@( if [ $(COMMS) == "mpi" ] ; \ 
#		then TEMP1 = $(MPIF90) ; \
#		TEMP2 = -DMPI ; \
#	else 
#		TEMP1 =  ; \
#		TEMP2 =  ; \
#	fi ) ;

prog libs: POSTOPTS = 
prog libs: COMPILER = $(F90)
prog: $(OBJS) 
	$(COMPILER) ../wannier_prog.F90 $(LDOPTS) $(OBJS) $(LIBS) -o ../../wannier90.x

post: POSTOPTS = $(TEMP1)
post: COMPILER = $(TEMP2)
post: mpi_test $(OBJS_POST)
	$(COMPILER) $(POSTDIR)postw90.F90 $(POSTOPTS) $(LDOPTS) $(OBJS_POST) $(LIBS) -o ../../postw90.x

mpi_test :
ifeq ($(COMMS),mpi)

ifndef MPIF90

	$(info  ) \
	$(info ********************* STOP **********************) \
	$(info *  COMMS=mpi but MPIF90 is not set in make.sys  *) \
	$(info *  Set MPIF90 to your parallel fortran compiler *) \
	$(info *  On many systems MPIF90=mpif90 will work      *) \
	$(info ********************* STOP **********************) \
	$(info  ) \
	$(error MPIF90 is not set) 

endif

endif 



libs:     $(LIBRARY)

$(LIBRARY): $(OBJS) $(OBJS2) 
	$(AR) $(ARFLAGS) $(LIBRARY) $(OBJS2) $(OBJS)

clean:
	rm -f *.o *.mod *.MOD *.obj

constants.o: ../constants.F90
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../constants.F90

io.o: ../io.F90 constants.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../io.F90

utility.o: ../utility.F90  constants.o io.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../utility.F90

parameters.o: ../parameters.F90 constants.o io.o utility.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../parameters.F90

hamiltonian.o: ../hamiltonian.F90 constants.o io.o utility.o parameters.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../hamiltonian.F90

overlap.o: ../overlap.F90 constants.o io.o utility.o parameters.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../overlap.F90

kmesh.o: ../kmesh.F90 constants.o io.o utility.o parameters.o 
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../kmesh.F90

disentangle.o: ../disentangle.F90 constants.o io.o parameters.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../disentangle.F90

wannierise.o: ../wannierise.F90 constants.o io.o utility.o parameters.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../wannierise.F90

plot.o: ../plot.F90 constants.o io.o utility.o parameters.o hamiltonian.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../plot.F90

transport.o: ../transport.F90 constants.o io.o parameters.o hamiltonian.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c ../transport.F90

comms.o: $(POSTDIR)comms.F90 constants.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)comms.F90

get_oper.o: $(POSTDIR)get_oper.F90 parameters.o constants.o comms.o wanint_common.o io.o utility.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)get_oper.F90

spin_wanint.o: $(POSTDIR)spin_wanint.F90 comms.o parameters.o constants.o utility.o wanint_common.o get_oper.o io.o wan_ham.o
	 $(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)spin_wanint.F90

dos_wanint.o: $(POSTDIR)dos_wanint.F90 comms.o parameters.o constants.o utility.o io.o wan_ham.o wanint_common.o get_oper.o spin_wanint.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)dos_wanint.F90

kpath_plot.o: $(POSTDIR)kpath_plot.F90 comms.o parameters.o constants.o io.o spin_wanint.o berry_wanint.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)kpath_plot.F90

slice_plot.o: $(POSTDIR)slice_plot.F90 comms.o parameters.o constants.o io.o spin_wanint.o berry_wanint.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)slice_plot.F90

berry_wanint.o: $(POSTDIR)berry_wanint.F90 comms.o parameters.o constants.o utility.o wanint_common.o get_oper.o io.o spin_wanint.o wan_ham.o
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)berry_wanint.F90

wan_ham.o: $(POSTDIR)wan_ham.F90 parameters.o constants.o utility.o wanint_common.o 
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)wan_ham.F90

wanint_common.o: $(POSTDIR)wanint_common.F90 comms.o parameters.o utility.o constants.o io.o 
	$(COMPILER) $(POSTOPTS) $(FCOPTS) -c $(POSTDIR)wanint_common.F90

wannier_lib.o: ./constants.o ./io.o ./utility.o ./parameters.o \
	./hamiltonian.o ./kmesh.o ./overlap.o ./disentangle.o \
	./wannierise.o ./plot.o ./transport.o ../wannier_lib.F90
	 $(COMPILER) $(FCOPTS) -c ../wannier_lib.F90


.PHONY : test_mpi 




