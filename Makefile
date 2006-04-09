

default: wannier

all: wannier

spec: thespec

wannier: 
	(cd src ; make all)

debug: $(OBJS) 
	$(F90) wannier_prog.f90 $(LDOPTS) $(OBJS) $(LIBS) -o wannier90.x


clean:
	touch make.sys ;\
	cd src ; make clean ;\
	cd ../tests ; make clean ; \
	cd ../spec ; make clean	

thespec:
	(cd spec ; make spec)

dist:
	@(tar cf - \
		./src/*.?90 \
		./tests/run_test.pl \
		./tests/test*/wannier.win \
		./tests/test*/des.dat \
		./tests/test*/wannier.eig \
		./tests/test*/wannier.*mn \
		./examples/*/wannier.eig \
		./tests/test*/stnd* \
		./examples/*/wannier.win \
		./examples/*/old* \
		./examples/*/new* \
                ./examples/*/*.*mn \
                ./examples/*/*.nnkp \
                ./examples/*/README \
		./config/* \
                ./spec/*.tex \
                ./spec/*.eps \
                ./spec/*.fig \
		./*/Makefile \
		./Makefile \
        | gzip -c > \
                ./wannier90.tar.gz)

test:   wannier
	(cd tests ; make test )

dist-lite:
	@(tar cf - \
                ./src/?90 \
		./src/Makefile \
                ./Makefile \
        | gzip -c > \
                ./wannier90.tar.gz)


