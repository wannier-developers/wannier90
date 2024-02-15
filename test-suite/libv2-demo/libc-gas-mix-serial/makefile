ROOT=../../../
CXXFLAGS=-Wall -I$(ROOT)/src -pedantic -g
include $(ROOT)/make.inc
XLIBS=-lstdc++ $(ROOT)/libwan2.a
a.out: wannier.o demo.o 
	$(MPIF90) wannier.o demo.o $(LIBS) $(XLIBS)

