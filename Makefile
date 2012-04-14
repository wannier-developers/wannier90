REALMAKEFILE=../../src/Makefile.2

default: wannier post

all: wannier lib post

doc: thedoc

serialobjs: objdir
	(cd src/obj ; make -f $(REALMAKEFILE) serialobjs)

wannier: objdir serialobjs
	(cd src/obj ; make -f $(REALMAKEFILE) wannier)

lib: objdir serialobjs
	(cd src/obj ; make -f $(REALMAKEFILE) libs)

libs: lib

post: objdirp
	(cd src/objp ; make -f $(REALMAKEFILE) post)

clean:
	rm -f *~
	rm -f src/*~
	@( if [ -d src/obj ] ; \
		then cd src/obj ; \
		make -f $(REALMAKEFILE) clean ; \
		cd ../ ; rm -rf obj ; cd ../ ; \
	fi ; \
	if [ -d src/objp ] ; \
		then cd src/objp ; \
		make -f $(REALMAKEFILE) clean ; \
		cd ../ ; rm -rf objp ; cd ../ ; \
	fi )
	make -C tests clean
	make -C doc/user_guide clean
	make -C doc/tutorial clean

veryclean: clean
	rm -f wannier90.x postw90.x libwannier.a ; \
	cd doc ; rm -f user_guide.pdf tutorial.pdf ; \
	cd user_guide ; rm -f user_guide.ps ; \
	cd ../tutorial ; rm -f tutorial.ps 

thedoc:
	cd doc/user_guide ; make guide ; \
	cd ../tutorial ; make tutorial

dist:
	@(tar cf - \
		./src/*.?90 \
		./src/postw90/*.?90 \
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
		./*/Makefile.2 \
		./*/*/Makefile \
		./Makefile \
		./LICENCE \
		./README* \
		./CHANGE.log \
        | gzip -c > \
                ./wannier90.tar.gz)

test:   default
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

objdir: 
	@( if [ ! -d src/obj ] ; \
		then mkdir src/obj ; \
	fi ) ;

objdirp: 
	@( if [ ! -d src/objp ] ; \
		then mkdir src/objp ; \
	fi ) ;

.PHONY: wannier default all doc lib libs post clean veryclean thedoc dist test dist-lite objdir objdirp serialobjs
