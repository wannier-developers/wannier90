
Demonstrator for basic functionality of the Wannier90 v4 library using the C-interface.  The Wannier90 library must be built with WANNIER90_WITH_C=1 to provide this interface (see the installation instructions).

The program parses .win files for unit cell definitions, kpoints and number of Wannier functions before reading M and A matrices and eigenvalues from usual Wannier90 input files and perfoms a wannierisation.

Other variables in the .win file are also read and used (via the call to cinput_reader() function), such that the results should be identical to the standalone executable.

This executable can be invoked in various of the test directories, taking the name of the .win as the only argument (not all features and functions of the real wannier90.x have been duplicated in this demo, of course), eg:

(export EXE=`pwd`/wannier_c.x; cd ../tests/testw90_example01; eval $EXE gaas.win )
(export EXE=`pwd`/wannier_c.x; cd ../tests/testw90_example02; eval $EXE lead.win )
(export EXE=`pwd`/wannier_c.x; cd ../tests/testw90_example03; eval $EXE silicon.win )
(export EXE=`pwd`/wannier_c.x; cd ../tests/testw90_example04; eval $EXE copper.win )
(export EXE=`pwd`/wannier_c.x; cd ../tests/testw90_example05; eval $EXE diamond.win )
(export EXE=`pwd`/wannier_c.x; cd ../tests/testw90_example07; eval $EXE silane.win )
