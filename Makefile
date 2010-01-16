

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

veryclean: clean
	rm -f wannier90.x libwannier.a ; \
	cd doc ; rm -f user_guide.pdf tutorial.pdf ; \
	cd user_guide ; rm -f user_guide.ps ; \
	cd ../tutorial ; rm -f tutorial.ps

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
		./examples/README \
		./examples/example01/UNK* \
		./examples/*/*.win \
		./examples/example0[2-4]/*.eig \
                ./examples/example0[1-4]/*.*mn \
		./examples/example0[5-9]/*.scf \
		./examples/example1[0-3]/*.scf \
		./examples/example0[5-9]/*.nscf \
		./examples/example1[0-3]/*.nscf \
		./examples/example0[5-9]/*.pw2wan \
		./examples/example1[0-3]/*.pw2wan \
		./examples/example1[4-5]/defected/*.scf \
		./examples/example1[4-5]/defected/*.nscf \
		./examples/example1[4-5]/defected/*.win \
		./examples/example1[4-5]/defected/*.pw2wan \
                ./examples/example1[4-5]/periodic/*.scf \
                ./examples/example1[4-5]/periodic/*.nscf \
                ./examples/example1[4-5]/periodic/*.win \
                ./examples/example1[4-5]/periodic/*.pw2wan \
		./examples/example08/README \
		./pseudo/*.UPF \
		./pwscf/README \
		./pwscf/v3.2.3/*.f90 \
		./pwscf/v4.0/*.f90 \
		./pwscf/v4.1/*.f90 \
		./config/make.sys* \
		./utility/*.pl \
		./utility/PL_assessment/*.f90 \
		./utility/PL_assessment/README \
                ./doc/*/*.tex \
                ./doc/*/*.eps \
                ./doc/*/*.fig \
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


