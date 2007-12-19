

default: wannier

all: wannier lib

doc: thedoc

wannier: 
	(cd src ; make prog)

debug: $(OBJS) 
	$(F90) wannier_prog.f90 $(LDOPTS) $(OBJS) $(LIBS) -o wannier90.x

lib:
	(cd src ; make libs)

libs: lib

clean:
	cd src ; rm -f *.o *.mod *.MOD *.obj ;\
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
		./examples/example[5-9]/*.scf \
		./examples/example1[0-3]/*.scf \
		./examples/example[5-9]/*.nscf \
		./examples/example1[0-3]/*.nscf \
		./examples/example[5-9]/*.pw2wan \
		./examples/example1[0-3]/*.pw2wan \
		./examples/example[5-9]/*.UPF \
		./examples/example1[0-3]/*.UPF \
		./config/make.sys* \
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
		./CHANGE.log \
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
		./CHANGE.log \
	| gzip -c > \
		./wannier90.tar.gz)


