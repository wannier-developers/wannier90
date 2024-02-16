ROOT=../../../
include $(ROOT)/make.inc
XINCS=-J$(ROOT)/src/obj
XLIBS=$(ROOT)/libwan2.a
a.out: demo.o 
	$(MPIF90) demo.o $(LIBS) $(XLIBS)
demo.o: demo.f90
	$(MPIF90) $(FCOPTS) $(XINCS) -c $<

