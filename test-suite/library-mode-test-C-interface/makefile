include ../../make.inc
.cc.o:
	$(CXX) $(FCOPTS) -g -I ../../src -c $<  -fpermissive
wannier_c.x : wannier.o ../../libwannier90.a
	$(CXX) $(LDOPTS) -g wannier.o ../../libwannier90.a $(LIBS) -lgfortran -o $@ -fpermissive
clean:
	@rm -f wannier.o wannier_c.x
