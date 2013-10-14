ifndef ROOTDIR
ROOTDIR=.
endif

REALMAKEFILE=../Makefile.2

default: wannier post

all: wannier lib post w90chk2chk

doc: thedoc

serialobjs: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) serialobjs)

w90chk2chk: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) w90chk2chk)

wannier: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) wannier)

lib: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) libs)

w90pov:
	(cd $(ROOTDIR)/utility/w90pov && $(MAKE) )

w90vdw:
	(cd $(ROOTDIR)/utility/w90vdw && $(MAKE) )

libs: lib

post: objdirp
	(cd $(ROOTDIR)/src/objp && $(MAKE) -f $(REALMAKEFILE) post)

clean:
	cd $(ROOTDIR) && rm -f *~
	cd $(ROOTDIR) && rm -f src/*~
	@( cd $(ROOTDIR) && if [ -d src/obj ] ; \
		then cd src/obj && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf obj ; \
	fi )
	@( cd $(ROOTDIR) && if [ -d src/objp ] ; \
		then cd src/objp && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf objp ; \
	fi )
	$(MAKE) -C $(ROOTDIR)/tests clean
	$(MAKE) -C $(ROOTDIR)/doc/user_guide clean
	$(MAKE) -C $(ROOTDIR)/doc/tutorial clean
	$(MAKE) -C $(ROOTDIR)/utility/w90pov clean

veryclean: clean
	cd $(ROOTDIR) && rm -f wannier90.x postw90.x libwannier.a
	cd $(ROOTDIR)/doc && rm -f user_guide.pdf tutorial.pdf
	cd $(ROOTDIR)/doc/user_guide && rm -f user_guide.ps
	cd $(ROOTDIR)/doc/tutorial && rm -f tutorial.ps 

thedoc:
	$(MAKE) -C $(ROOTDIR)/doc/user_guide 
	$(MAKE) -C $(ROOTDIR)/doc/tutorial 

dist:
	@(cd $(ROOTDIR) && tar -cz --transform='s,^\./,wannier90-2.0.0/,' -f wannier90-2.0.0.tar.gz \
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
		./examples/example1[6-9]/*.scf \
		./examples/example1[6-9]/*.nscf \
		./examples/example1[6-9]/*.pw2wan \
		./pseudo/*.UPF \
		./pwscf/README \
		./pwscf/v3.2.3/*.f90 \
		./pwscf/v4.0/*.f90 \
		./pwscf/v4.1/*.f90 \
		./pwscf/v5.0/*.f90 \
		./config/make.sys* \
		./utility/*.pl \
		./utility/PL_assessment/*.f90 \
		./utility/PL_assessment/README \
		./utility/w90vdw/w90vdw.f90 \
		./utility/w90vdw/README \
		./utility/w90vdw/doc/Makefile \
		./utility/w90vdw/doc/w90vdw.tex \
		./utility/w90vdw/examples/benzene_s_val/benzene_s_val.* \
                ./utility/w90vdw/examples/benzene_s_val/ref/benzene_s_val.* \
                ./utility/w90vdw/examples/benzene_s_cond/benzene_s_cond.* \
                ./utility/w90vdw/examples/benzene_s_cond/ref/benzene_s_cond.* \
		./utility/w90pov/doc/*.tex \
		./utility/w90pov/doc/*.pdf \
		./utility/w90pov/doc/figs/*.png \
		./utility/w90pov/src/*.f90 \
		./utility/w90pov/src/*.c \
		./utility/w90pov/examples/*/*.gz \
		./utility/w90pov/examples/*/*.inp \
		./utility/w90pov/examples/*/*.inc \
		./utility/w90pov/examples/*/*.pov \
		./utility/w90pov/examples/*/ref/*.png \
		./utility/w90pov/README \
		./utility/w90chk2chk/README \
                ./doc/*/*.tex \
                ./doc/*/*.eps \
                ./doc/*/*.fig \
                ./doc/user_guide.pdf \
                ./doc/tutorial.pdf \
		./doc/wannier90.bib \
		./*/Makefile \
		./*/Makefile.2 \
		./*/*/Makefile \
		./Makefile \
		./LICENCE \
		./README* \
		./CHANGE.log \
	)

test:   default
	(cd $(ROOTDIR)/tests && $(MAKE) test )

dist-lite:
	@(cd $(ROOTDIR) && tar -cz --transform='s,^\./,wannier90/,' -f wannier90.tar.gz \
		./src/*.?90 \
		./src/postw90/*.?90 \
		./tests/run_test.pl \
		./tests/test*/wannier.win \
		./tests/test*/des.dat \
		./tests/test*/wannier.eig \
		./tests/test*/wannier.?mn \
		./tests/test*/stnd* \
		./config/* \
		./*/Makefile \
		./utility/*.pl \
		./*/Makefile \
		./*/Makefile.2 \
		./*/*/Makefile \
		./Makefile \
		./LICENCE \
		./README.* \
		./CHANGE.log \
	)

objdir: 
	@( cd $(ROOTDIR) && if [ ! -d src/obj ] ; \
		then mkdir src/obj ; \
	fi ) ;

objdirp: 
	@( cd $(ROOTDIR) && if [ ! -d src/objp ] ; \
		then mkdir src/objp ; \
	fi ) ;

.PHONY: wannier default all doc lib libs post clean veryclean thedoc dist test dist-lite objdir objdirp serialobjs
