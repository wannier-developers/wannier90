#  Should be no need to change below this line
#

include ../make.sys

OBJS  =  parameters.o io.o utility.o overlap.o kmesh.o  \
         disentangle.o wannierise.o plot.o

OBJS2  =  wannier_lib.o

LIBRARY = ../libwannier.a

prog: $(OBJS) 
	$(F90) wannier_prog.F90 $(LDOPTS) $(OBJS) $(LIBS) -o ../wannier90.x

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

overlap.o: overlap.F90 parameters.o constants.o io.o
	 $(F90) $(FCOPTS) -c overlap.F90

kmesh.o: kmesh.F90 utility.o parameters.o constants.o io.o
	 $(F90) $(FCOPTS) -c kmesh.F90

disentangle.o: disentangle.F90 parameters.o constants.o io.o
	 $(F90) $(FCOPTS) -c disentangle.F90

wannierise.o: wannierise.F90 parameters.o io.o utility.o constants.o io.o
	 $(F90) $(FCOPTS) -c wannierise.F90

plot.o: plot.F90
	 $(F90) $(FCOPTS) -c plot.F90

wannier_lib.o: ./plot.o ./wannierise.o ./overlap.o \
	./disentangle.o ./kmesh.o ./io.o ./parameters.o \
	./constants.o wannier_lib.F90
	 $(F90) $(FCOPTS) -c wannier_lib.F90


