

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

dist-full:
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
		./examples/*/*.eig \
                ./examples/*/*.*mn \
		./examples/*/*.scf \
		./examples/*/*.nscf \
		./examples/*/*.pw2wan \
		./examples/*/*.UPF \
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
		./README.* \
        | gzip -c > \
                ./wannier90.tar.gz)

test:   wannier
	(cd tests ; make test )

dist:
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


