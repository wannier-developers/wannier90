include ../../make.inc
.cc.o:
	$(CXX) -I ../../src -c $< 
wannier_c.x : wannier.o ../../libwan2.a
	$(CXX)  wannier.o ../../libwan2.a $(LIBS) -lgfortran -o $@
clean:
	@rm wannier.o wannier_c.x
