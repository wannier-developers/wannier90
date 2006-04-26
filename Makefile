

default: wannier

all: wannier

doc: thedoc

wannier: 
	(cd src ; make all)

debug: $(OBJS) 
	$(F90) wannier_prog.f90 $(LDOPTS) $(OBJS) $(LIBS) -o wannier90.x


clean:
	touch make.sys ;\
	cd src ; make clean ;\
	cd ../tests ; make clean ; \
	cd ../doc/user_guide ; make clean ; \
	cd ../tutorial ; make clean 	

thedoc:
	cd doc/user_guide ; make guide ; \
	cd ../tutorial ; make tutorial

dist:
	@(tar cf - \
		./src/*.?90 \
		./tests/run_test.pl \
		./tests/test*/wannier.win \
		./tests/test*/des.dat \
		./tests/test*/wannier.eig \
		./tests/test*/wannier.?mn \
		./tests/test*/stnd* \
		./examples/example1/UNK* \
		./examples/*/*.win \
		./examples/example[2-4]/*.eig \
                ./examples/example[1-4]/*.*mn \
		./examples/example5/*.scf \
		./examples/example5/*.nscf \
		./examples/example5/*.pw2wan \
		./examples/example5/*.UPF \
		./config/* \
		./utility/*.pl \
                ./doc/*/*.tex \
                ./doc/*/*.eps \
                ./doc/*/*.fig \
		./doc/*.pdf \
		./*/Makefile \
		./*/*/Makefile \
		./Makefile \
		./LICENCE \
		./README* \
        | gzip -c > \
                ./wannier90.tar.gz)

test:   wannier
	(cd tests ; make test )

dist-lite:
	@(tar cf - \
		./src/*.?90 \
		./tests/run_test.pl \
		./tests/test*/wannier.win \
		./tests/test*/des.dat \
		./tests/test*/wannier.eig \
		./tests/test*/wannier.?mn \
		./tests/test*/stnd* \
		./config/* \
		./*/Makefile \
		./utility/*.pl \
		./Makefile \
		./LICENCE \
		./README.* \
	| gzip -c > \
		./wannier90.tar.gz)


