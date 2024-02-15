ROOT=../../../
CXX=mpicxx
CXXFLAGS=-Wall -I$(ROOT)/src -g
include $(ROOT)/make.inc
XLIBS=-lstdc++ $(ROOT)/libwan2.a -lmpi -lmpi_cxx
a.out: wannier.o demo.o 
	$(MPIF90) wannier.o demo.o $(LIBS) $(XLIBS)

