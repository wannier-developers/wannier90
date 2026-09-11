ifndef ROOTDIR
ROOTDIR=.
endif

# include make.inc to determine if (last) build was serial or parallel via def/undef COMMS
include make.inc

# Contains definition of OBJS, OBJS_POST, LIBRARY, DYNLIBRARY, ...
include Makefile.header

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

	install -d $(DESTDIR)$(PREFIX)/include/
	for m in $(addprefix src/obj/,$(LIBMODS)); do \
		if [ -f "$$m" ]; then install -m644 "$$m" "$(DESTDIR)$(PREFIX)/include/"; fi; \
	done
	install -d $(DESTDIR)$(PREFIX)/lib/
	if [ -f "$(STATICLIBRARY)" ]; then install -m644 "$(STATICLIBRARY)" "$(DESTDIR)$(PREFIX)/lib/$(STATICLIBRARY)"; fi;
	if [ -f "$(DYNLIBRARY)" ]; then install -m644 "$(DYNLIBRARY)" "$(DESTDIR)$(PREFIX)/lib/$(DYNLIBRARY)"; fi;
	if [ -f "$(STATICLIBRARY)" ]; then $(MAKE) pkgconfig; fi;

all: wannier libs post w90chk2chk w90pov w90vdw w90spn2spn

doc: thedoc

w90chk2chk:
	$(MAKE) -C src/obj w90chk2chk

w90spn2spn:
	$(MAKE) -C src/obj w90spn2spn

wannier:
	$(MAKE) -C src/obj wannier

# General rule to make the wannier90.x, postw90.x, w90chk2chk.x and w90spn2spn.x executables
# Internally it uses ../$@ because in the src/ directory, the executable is created one level up
# (i.e. in the root directory)
%.x:
	$(MAKE) -C src/obj ../$@

staticlib:
	$(MAKE) -C src/obj staticlib

dynlib:
	$(MAKE) -C src/obj dynlib

w90pov:
	$(MAKE) -C $(ROOTDIR)/utility/w90pov

w90vdw:
	$(MAKE) -C $(ROOTDIR)/utility/w90vdw

w90py: libs
	$(MAKE) -C $(ROOTDIR)/test-suite/library/py-f90wrap

libs: staticlib dynlib

PKGCONFIG_FILENAME = $(DYNLIBBASE).pc
pkgconfig:
	{ \
	  echo "prefix=$(DESTDIR)$(PREFIX)"; \
	  echo "exec_prefix=$(DESTDIR)$(PREFIX)/bin"; \
	  echo "libdir=$(DESTDIR)$(PREFIX)/lib"; \
	  echo "includedir=$(DESTDIR)$(PREFIX)/include"; \
	  echo ""; \
	  echo "Name: $(DYNLIBBASE)"; \
	  echo "Description: $(LIBDESCRIPTION)."; \
	  echo "Requires: "; \
	  echo "Version: $(VERSION)"; \
	  echo 'Libs: -L$${libdir} -l'"$(DYNLIBBASE)"; \
	  echo 'Cflags: -I$${includedir}'; \
	} > "$(PKGCONFIG_FILENAME)"
	install -d $(DESTDIR)$(PREFIX)/lib/pkgconfig/
	install -D -m644 "$(PKGCONFIG_FILENAME)" "$(DESTDIR)$(PREFIX)/lib/pkgconfig/$(PKGCONFIG_FILENAME)"
	cd $(ROOTDIR) && rm -f $(PKGCONFIG_FILENAME)

post:
	$(MAKE) -C src/obj post

clean:
	cd $(ROOTDIR) && rm -f *~
	cd $(ROOTDIR) && rm -f src/*~
	cd $(ROOTDIR) && rm -f $(PKGCONFIG_FILENAME)
	cd $(ROOTDIR) && $(MAKE) -C src/obj clean
	$(MAKE) -C $(ROOTDIR)/utility/w90pov clean
	$(MAKE) -C $(ROOTDIR)/utility/w90vdw clean

# Note: .x.dSYM are directories (hence the -r option to rm) and are only created on macOS (when compiling with certain flags, e.g. debug), so they are not always present
veryclean: clean
	cd $(ROOTDIR) && rm -rf wannier90.x postw90.x w90chk2chk.x w90spn2spn.x libwannier90.{a,so.4} libwannier90.{a,so.4} *.{gcda,gcno} *.x.dSYM

thedoc:
	@(echo "The latex user_guide and tutorials have been migrated to markdown \
	format, for more details see 'docs/README.md' file.")

# Create a tarball of the current repo without git files (useful to distribute the code outside of git)
dist:
	cd $(ROOTDIR) && git archive HEAD --prefix=wannier90-current/ -o wannier90-current.tar.gz

# The test suite is driven by pytest; see test-suite/README.md.
# Requires Python >= 3.10 with pytest and PyYAML:
#   pip install -r test-suite/requirements.txt
PYTHON ?= python3

test-serial: w90chk2chk wannier post
	(cd $(ROOTDIR)/test-suite && $(PYTHON) -m pytest tests )

test-parallel: w90chk2chk wannier post
	(cd $(ROOTDIR)/test-suite && $(PYTHON) -m pytest tests --nprocs=4 \
		-m "wannier90 or postw90 or checkpoint or parallel" )

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

.PHONY: wannier default all doc libs staticlib dynlib post clean veryclean thedoc dist test-serial test-parallel dist-lite tests w90spn2spn install pkgconfig
