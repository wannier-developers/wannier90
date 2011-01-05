#  Should be no need to change below this line
#

include ../make.sys

OBJS  =  constants.o io.o utility.o parameters.o hamiltonian.o overlap.o \
	 kmesh.o disentangle.o wannierise.o plot.o transport.o

OBJS2  =  wannier_lib.o

OBJS_PAR  = constants.o io.o utility.o parameters.o hamiltonian.o overlap.o kmesh.o comms.o wanint_common.o 


LIBRARY = ../libwannier.a

prog: $(OBJS) 
	$(F90) wannier_prog.F90 $(LDOPTS) $(OBJS) $(LIBS) -o ../wannier90.x

wanpar: $(OBJS_PAR)
	$(F90) wanpar.F90 $(LDOPTS) $(OBJS_PAR) $(LIBS) -o ../wanpar.mpi


libs:     $(LIBRARY)

$(LIBRARY): $(OBJS) $(OBJS2) 
	$(AR) $(ARFLAGS) $(LIBRARY) $(OBJS2) $(OBJS)



clean:
	rm -f *.o *.mod *.MOD *.obj

constants.o: constants.F90
	 $(F90) $(FCOPTS) -c constants.F90

io.o: io.F90 constants.o
	 $(F90) $(FCOPTS) -c io.F90

utility.o: utility.F90  constants.o io.o
	 $(F90) $(FCOPTS) -c utility.F90

parameters.o: parameters.F90 constants.o io.o utility.o
	 $(F90) $(FCOPTS) -c parameters.F90

hamiltonian.o: hamiltonian.F90 constants.o io.o utility.o parameters.o
	 $(F90) $(FCOPTS) -c hamiltonian.F90

overlap.o: overlap.F90 constants.o io.o utility.o parameters.o
	 $(F90) $(FCOPTS) -c overlap.F90

kmesh.o: kmesh.F90 constants.o io.o utility.o parameters.o 
	 $(F90) $(FCOPTS) -c kmesh.F90

disentangle.o: disentangle.F90 constants.o io.o parameters.o
	 $(F90) $(FCOPTS) -c disentangle.F90

wannierise.o: wannierise.F90 constants.o io.o utility.o parameters.o
	 $(F90) $(FCOPTS) -c wannierise.F90

plot.o: plot.F90 constants.o io.o utility.o parameters.o hamiltonian.o
	 $(F90) $(FCOPTS) -c plot.F90

transport.o: transport.F90 constants.o io.o parameters.o hamiltonian.o
	 $(F90) $(FCOPTS) -c transport.F90

comms.o: comms.F90 constants.o
	 $(F90) $(FCOPTS) -c comms.F90

wanint_common.o: wanint_common.F90 comms.o parameters.o utility.o constants.o io.o 
	$(F90) $(FCOPTS) -c wanint_common.F90

wannier_lib.o: ./constants.o ./io.o ./utility.o ./parameters.o \
	./hamiltonian.o ./kmesh.o ./overlap.o ./disentangle.o \
	./wannierise.o ./plot.o ./transport.o wannier_lib.F90
	 $(F90) $(FCOPTS) -c wannier_lib.F90


