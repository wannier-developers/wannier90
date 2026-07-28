ROOT=../../
include ../../make.inc
# -DMPI -DMPI90/-DMPIH/-DMPI08 preprocessor flags (POSTOPTS)
include ../../Makefile.header

# module include path specifier "-J" is gcc specific
# compare_results compares the printed final centres and spreads with those obtained with the standalone wannier90.x executable (see directory "ref")
results.dat: demo.x
	./demo.x > $@
	./compare_results.py
demo.x: demo.o
	$(COMPILER) $< -o $@ $(POSTOPTS) $(LIBS) -L$(ROOT) -l$(DYNLIBBASE)
demo.o: demo.F90
	$(COMPILER) $(POSTOPTS) -J$(ROOT)/src/obj -c $<
clean:
	@rm -f demo.{o,x} results.dat
