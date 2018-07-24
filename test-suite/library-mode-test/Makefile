test_library.x: test_library.F90 ../../libwannier.a
	mpif90.openmpi -fstrict-aliasing  -fno-omit-frame-pointer -fno-realloc-lhs -fcheck=bounds,do,recursion,pointer -ffree-form -Wall -Waliasing -Wsurprising -Wline-truncation -Wno-tabs -Wno-uninitialized -Wno-unused-dummy-argument -Wno-unused -Wno-character-truncation -O1 -g -fbacktrace -fPIC $< ../../libwannier.a -lblas -llapack -o $@

clean:
	rm -f test_library.x *~

.PHONY: clean
