#ifndef __wannier_hh__
#define __wannier_hh__

#include <complex>

void wannier_setup(void*&, double*, int*, int, int, int, double*, int&, int*);
void wannier_run(void*, std::complex<double>*, double*, std::complex<double>*, std::complex<double>*, double*, double*);

#endif
