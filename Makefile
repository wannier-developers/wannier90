ifndef ROOTDIR
ROOTDIR=.
endif

REALMAKEFILE=../Makefile.2

TAR := $(shell if which gnutar 1>/dev/null 2> /dev/null; then echo gnutar; else echo tar; fi )

default: wannier post

all: wannier lib post w90chk2chk w90pov w90vdw w90spn2spn

doc: thedoc

serialobjs: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) serialobjs)

w90chk2chk: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) w90chk2chk)

w90spn2spn: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) w90spn2spn)

wannier: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) wannier)

lib: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) libs)

dynlib: objdir serialobjs
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) dynlibs)

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
	$(MAKE) -C $(ROOTDIR)/doc/user_guide clean
	$(MAKE) -C $(ROOTDIR)/doc/tutorial clean
	$(MAKE) -C $(ROOTDIR)/utility/w90pov clean
	$(MAKE) -C $(ROOTDIR)/utility/w90vdw clean
	cd $(ROOTDIR)/test-suite && ./clean_tests

veryclean: clean
	cd $(ROOTDIR) && rm -f wannier90.x postw90.x libwannier.a w90chk2chk.x w90spn2spn.x
	cd $(ROOTDIR)/doc && rm -f user_guide.pdf tutorial.pdf
	cd $(ROOTDIR)/doc/user_guide && rm -f user_guide.ps
	cd $(ROOTDIR)/doc/tutorial && rm -f tutorial.ps 
	cd $(ROOTDIR)/test-suite && ./clean_tests -i

thedoc:
	$(MAKE) -C $(ROOTDIR)/doc/user_guide 
	$(MAKE) -C $(ROOTDIR)/doc/tutorial 

# For now hardcoded to 3.0.0, and using HEAD
# Better to get the version from the io.F90 file and use
# the tag (e.g. v3.0.0) instead of HEAD
dist: 
	cd $(ROOTDIR) && git archive HEAD --prefix=wannier90-3.0.0/ -o wannier90-3.0.0.tar.gz

dist-legacy:
	@(cd $(ROOTDIR) && $(TAR) -cz --transform='s,^\./,wannier90-3.0/,' -f wannier90-3.0.tar.gz \
		./src/*.?90 \
		./src/postw90/*.?90 \
		./autodoc/README.txt \
		./autodoc/*.md \
		./autodoc/media/favicon*png \
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
                ./examples/example16-noqe/Si.amn \
                ./examples/example16-noqe/Si.mmn \
                ./examples/example16-noqe/Si.eig \
                ./examples/example16-withqe/Si.scf \
                ./examples/example16-withqe/Si.nscf \
                ./examples/example16-withqe/Si.pw2wan \
		./examples/example1[7-9]/*.scf \
		./examples/example1[7-9]/*.nscf \
		./examples/example1[7-9]/*.pw2wan \
		./examples/example20/*.scf \
		./examples/example20/*.nscf \
		./examples/example20/*.pw2wan \
		./examples/example20/SrMnO3/SrMnO3-d.pw2wan \
		./examples/example20/SrMnO3/SrMnO3-d.win \
		./examples/example20/SrMnO3/SrMnO3-eg.pw2wan \
		./examples/example20/SrMnO3/SrMnO3-eg.win \
		./examples/example20/SrMnO3/SrMnO3.nscf \
		./examples/example20/SrMnO3/SrMnO3.scf \
		./examples/example20/SrMnO3/SrMnO3-t2g.pw2wan \
		./examples/example20/SrMnO3/SrMnO3-t2g.win \
		./examples/example2[1-2]/README \
		./examples/example2[1-2]/*/*.scf \
		./examples/example2[1-2]/*/*.nscf \
		./examples/example2[1-2]/*/*.win \
		./examples/example2[1-2]/*/*.sym \
		./examples/example2[1-2]/*/*.pw2wan \
		./pseudo/*.UPF \
		./pwscf/README \
		./pwscf/v*/*.f90 \
		./pwscf/v*/README \
		./config/make.inc* \
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
		./doc/wannier90.bib \
		./*/Makefile \
		./*/Makefile.2 \
		./*/*/Makefile \
		./Makefile \
		./LICENSE \
		./README* \
		./CHANGE.log \
	)

test-serial: w90chk2chk wannier post  
	(cd $(ROOTDIR)/test-suite && ./run_tests --category=default )

test-parallel: w90chk2chk wannier post 
	(cd $(ROOTDIR)/test-suite && ./run_tests --category=default --numprocs=4 )

# Alias
tests: test-serial test-parallel

dist-lite:
	@(cd $(ROOTDIR) && $(TAR) -cz --transform='s,^\./,wannier90/,' -f wannier90.tar.gz \
		./src/*.?90 \
		./src/postw90/*.?90 \
		./config/* \
		./*/Makefile \
		./utility/*.pl \
		./*/Makefile \
		./*/Makefile.2 \
		./*/*/Makefile \
		./Makefile \
		./LICENSE \
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

.PHONY: wannier default all doc lib libs post clean veryclean thedoc dist test-serial test-parallel dist-lite objdir objdirp serialobjs tests w90spn2spn
