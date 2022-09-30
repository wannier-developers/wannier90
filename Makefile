ifndef ROOTDIR
ROOTDIR=.
endif

# include make.inc to determine if (last) build was serial or parallel via def/undef COMMS
include make.inc

REALMAKEFILE=../Makefile.2

TAR := $(shell if which gnutar 1>/dev/null 2> /dev/null; then echo gnutar; else echo tar; fi )

.NOTPARALLEL:
default: wannier post

PREFIX ?= /usr

VERSION_MAJOR = 3
VERSION_MINOR = 1
VERSION_PATCH = 0

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)
VERSION_SHORT = $(VERSION_MAJOR).$(VERSION_MINOR)

install: default
	install -d $(DESTDIR)$(PREFIX)/bin/
	for x in wannier90.x postw90.x w90chk2chk.x w90spn2spn.x ; do \
		if [ -f "$$x" ]; then install -m755 "$$x" "$(DESTDIR)$(PREFIX)/bin/$$x"; fi; \
	done
	if [ -f "utility/w90pov/w90pov" ]; then install -m755 "utility/w90pov/w90pov" "$(DESTDIR)$(PREFIX)/bin/w90pov"; fi;
	if [ -f "utility/w90vdw/w90vdw.x" ]; then install -m755 "utility/w90vdw/w90vdw.x" "$(DESTDIR)$(PREFIX)/bin/w90vdw.x"; fi;
	install -d $(DESTDIR)$(PREFIX)/lib/
	if [ -f "$(LIBRARYV2)" ]; then install -m644 "$(LIBRARYV2)" "$(DESTDIR)$(PREFIX)/lib/$(LIBRARYV2)"; fi;
	if [ -f "$(LIBRARYV2)" ]; then $(MAKE) pkgconfig; fi;

all: wannier lib post w90chk2chk w90pov w90vdw w90spn2spn

doc: thedoc

w90chk2chk: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) w90chk2chk)

w90spn2spn: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) w90spn2spn)

wannier: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) wannier)

lib: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) libs)

dynlib: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) dynlibs)

w90pov:
	(cd $(ROOTDIR)/utility/w90pov && $(MAKE) )

w90vdw:
	(cd $(ROOTDIR)/utility/w90vdw && $(MAKE) )

w90py: dynlib
	(cd $(ROOTDIR)/wrap && $(MAKE) )

libs: lib

PKGCONFIG_FILENAME = wannier.pc
pkgconfig:
	$(file > $(PKGCONFIG_FILENAME),prefix=$(DESTDIR)$(PREFIX))
	$(file >> $(PKGCONFIG_FILENAME),exec_prefix=$(DESTDIR)$(PREFIX)/bin)
	$(file >> $(PKGCONFIG_FILENAME),libdir=$(DESTDIR)$(PREFIX)/lib)
	$(file >> $(PKGCONFIG_FILENAME),includedir=$(DESTDIR)$(PREFIX)/include)
	$(file >> $(PKGCONFIG_FILENAME),)
	$(file >> $(PKGCONFIG_FILENAME),Name: wannier)
	$(file >> $(PKGCONFIG_FILENAME),Description: Compute maximally-localised Wannier functions.)
	$(file >> $(PKGCONFIG_FILENAME),Requires: )
	$(file >> $(PKGCONFIG_FILENAME),Version: $(VERSION))
	$(file >> $(PKGCONFIG_FILENAME),Libs: -L$${libdir} -lwannier)
	$(file >> $(PKGCONFIG_FILENAME),Cflags: -I$${includedir})
	install -D -m644 "$(PKGCONFIG_FILENAME)" "$(DESTDIR)$(PREFIX)/lib/pkgconfig/$(PKGCONFIG_FILENAME)"
	cd $(ROOTDIR) && rm -f $(PKGCONFIG_FILENAME)

post: objdir
	(cd $(ROOTDIR)/src/obj && $(MAKE) -f $(REALMAKEFILE) post)

clean:
	cd $(ROOTDIR) && rm -f *~
	cd $(ROOTDIR) && rm -f src/*~
	cd $(ROOTDIR) && rm -f $(PKGCONFIG_FILENAME)
	@( cd $(ROOTDIR) && if [ -d src/obj ] ; \
		then cd src/obj && \
		$(MAKE) -f $(REALMAKEFILE) clean && \
		cd ../ && rm -rf obj ; \
	fi )
	$(MAKE) -C $(ROOTDIR)/utility/w90pov clean
	$(MAKE) -C $(ROOTDIR)/utility/w90vdw clean
	cd $(ROOTDIR)/test-suite && ./clean_tests

veryclean: clean
	cd $(ROOTDIR) && rm -f wannier90.x postw90.x w90chk2chk.x w90spn2spn.x libwannier90.{a,so.4} libwannier90_mpi.{a,so.4} *.{gcda,gcno}
	cd $(ROOTDIR)/test-suite && ./clean_tests -i

thedoc:
	@(echo "The latex user_guide and tutorials have been migrated to markdown \
	format, for more details see 'docs/README.md' file.")

# For now hardcoded to 3.1.0, and using HEAD
# Better to get the version from the io.F90 file and use
# the tag (e.g. v3.1.0) instead of HEAD
dist:
	cd $(ROOTDIR) && git archive HEAD --prefix=wannier90-3.1.0/ -o wannier90-3.1.0.tar.gz

dist-legacy:
	@(cd $(ROOTDIR) && $(TAR) -cz --transform='s,^\./,wannier90-3.1/,' -f wannier90-3.1.tar.gz \
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
	(cd $(ROOTDIR)/test-suite && ./run_tests --category=par --numprocs=4 )

# Alias
ifdef COMMS
tests: test-serial test-parallel
else
tests: test-serial
endif

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

.PHONY: wannier default all doc lib libs post clean veryclean thedoc dist test-serial test-parallel dist-lite objdir objdirp tests w90spn2spn install pkgconfig
