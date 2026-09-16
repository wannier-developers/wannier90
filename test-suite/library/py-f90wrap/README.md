This directory contains a Python wrapping of the Fortran library.

For f90-wrap, which this depends on see:

- [github](https://github.com/jameskermode/f90wrap)
- [DOI 10.1088/1361-648X/ab82d2](https://iopscience.iop.org/article/10.1088/1361-648X/ab82d2)

To experiment, please study the python files: serial-example.py , mpi-example.py , example-dos.py

- Make sure f90wrap is installed
- Specifiy `F90WRAP` in the build configuration (make.inc)
- Edit the makefile as appropriate.
