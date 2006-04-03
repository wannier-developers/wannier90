#  Should be no need to change below this line
#

include ../make.sys

OBJS  =  parameters.o io.o utility.o overlap.o kmesh.o  \
         disentangle.o wannierise.o plot.o


all: $(OBJS) 
	$(F90) wannier_prog.F90 $(LDOPTS) $(OBJS) $(LIBS) -o ../wannier90.x


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

kmesh.o: kmesh.F90 parameters.o constants.o io.o
	 $(F90) $(FCOPTS) -c kmesh.F90

disentangle.o: disentangle.F90 parameters.o constants.o io.o
	 $(F90) $(FCOPTS) -c disentangle.F90

wannierise.o: wannierise.F90 parameters.o io.o utility.o constants.o io.o
	 $(F90) $(FCOPTS) -c wannierise.F90

plot.o: plot.F90
	 $(F90) $(FCOPTS) -c plot.F90


