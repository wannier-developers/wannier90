
This crude file demonstrates basic functionality of the Wannier90 library (version 2).

It parses .win files for unit cell definitions, kpoints and number of Wannier functions before reading M and A matrices and eigenvalues from usual Wannier90 input files and perfoms a wannierisation.

Other variables in the .win file are also read and used (via the call to cinput_reader() function), such that the results should be identical to the standalone executable.

This executable can be invoked in the test directories, taking the name of the .win as the only argument, eg:

( cd ../tests/testw90_example01; ../../libv2-demo2/wannier_c.x gaas.win )
( cd ../tests/testw90_example02; ../../libv2-demo2/wannier_c.x lead.win )
( cd ../tests/testw90_example03; ../../libv2-demo2/wannier_c.x silicon.win )
( cd ../tests/testw90_example04; ../../libv2-demo2/wannier_c.x copper.win )
( cd ../tests/testw90_example05; ../../libv2-demo2/wannier_c.x diamond.win )
( cd ../tests/testw90_example07; ../../libv2-demo2/wannier_c.x silane.win )
