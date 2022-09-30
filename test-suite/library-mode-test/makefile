ROOT=../../
include $(ROOT)/make.inc
# module include path specifier "-J" is gcc specific
# compare_results compares the printed final centres and spreads with those obtained with the standalone wannier90.x executable (see directory "ref")
results.dat: demo.x 
	./demo.x > $@
	./compare_results.py
demo.x: demo.o 
	$(MPIF90) $< -o $@ $(LIBS) -L$(ROOT) -lwannier90
demo.o: demo.f90
	$(MPIF90) $(FCOPTS) -J$(ROOT)/src/obj -c $<
clean:
	@rm -f demo.{o,x} results.dat
