#ifndef __wannier_hh__
#define __wannier_hh__

#include <mpi.h>
#include <complex>

void wannier_setup(void*&, int*, double*, int*, int, int, int, double*, int*, int&, int*, MPI_Comm);
void wannier_run(void*, std::complex<double>*, double*, std::complex<double>*, std::complex<double>*, double*, double*);

#endif
