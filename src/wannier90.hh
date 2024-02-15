#ifndef __wannier90_hh__
#define __wannier90_hh__

#include <ISO_Fortran_binding.h>
#include <complex>
#include <cstring>

extern "C" {
void cchkpt(void*, CFI_cdesc_t*);

void cinput_reader_special(void*, CFI_cdesc_t*, int&);
void cinput_reader(void*, int&);
void cinput_setopt(void*, CFI_cdesc_t*, int&, int);
void cset_option_float33(void*, CFI_cdesc_t*, void*);
void cset_option_float(void*, CFI_cdesc_t*, double);
void cset_option_floatxy(void*, CFI_cdesc_t*, void*, int, int);
void cset_option_int3(void*, CFI_cdesc_t*, void*);
void cset_option_int(void*, CFI_cdesc_t*, int);
void cset_option_intx(void*, CFI_cdesc_t*, int*, int);

void coverlaps(void*, int&);
void cdisentangle(void*, int&);
void cwannierise(void*, int&);

void* w90_create();
void w90_delete(void*);

void cset_kpoint_distribution(void*, int*);
void cset_parallel_comms(void*, int);

void cset_a_matrix(void*, std::complex<double>*);
void cset_eigval(void*, double*);
void cset_m_local(void*, std::complex<double>*);
void cset_u_matrix(void*, std::complex<double>*);
void cset_u_opt(void*, std::complex<double>*);

void cget_nn(void*, int&);
void cget_nnkp(void*, int*);
void cget_gkpb(void*, int*);
void cget_centres(void*, void*);
void cget_spreads(void*, void*);
}

#ifdef MPI_VERSION
#include <mpi.h>
void cinput_setopt(void* blob, std::string seed, int& ierr, MPI_Comm comm) {
        int fcomm = MPI_Comm_c2f(comm);
        CFI_cdesc_t stringdesc;
        char* seedc = (char*)seed.c_str(); // discarding constness
        CFI_establish(&stringdesc, seedc, CFI_attribute_other, CFI_type_char, strlen(seedc), 0, NULL);
        ierr = 0;
        cinput_setopt(blob, &stringdesc, ierr, fcomm);
}
#else
void cinput_setopt(void* blob, std::string seed, int& ierr) {
        CFI_cdesc_t stringdesc;
        char* seedc = (char*)seed.c_str(); // discarding constness
        CFI_establish(&stringdesc, seedc, CFI_attribute_other, CFI_type_char, strlen(seedc), 0, NULL);
        ierr = 0;
        cinput_setopt(blob, &stringdesc, ierr, 0);
}
#endif

// see https://community.intel.com/t5/Intel-Fortran-Compiler/C-interoperablilty-and-character-strings/td-p/1084167
void cset_option(void* blob, std::string key, int x) {
        CFI_cdesc_t stringdesc;
        char* keyc = (char*)key.c_str(); // discarding constness
        CFI_establish(&stringdesc, keyc, CFI_attribute_other, CFI_type_char, strlen(keyc), 0, NULL);
        cset_option_int(blob, &stringdesc, x);
}
void cset_option(void* blob, std::string key, int x[3]) {
        CFI_cdesc_t stringdesc;
        char* keyc = (char*)key.c_str(); // discarding constness
        CFI_establish(&stringdesc, keyc, CFI_attribute_other, CFI_type_char, strlen(keyc), 0, NULL);
        cset_option_int3(blob, &stringdesc, &x);
}
void cset_option(void* blob, std::string key, int *x, int y) {
        CFI_cdesc_t stringdesc;
        char* keyc = (char*)key.c_str(); // discarding constness
        CFI_establish(&stringdesc, keyc, CFI_attribute_other, CFI_type_char, strlen(keyc), 0, NULL);
        cset_option_intx(blob, &stringdesc, x, y);
}
void cset_option(void* blob, std::string key, double x[][3]) {
        CFI_cdesc_t stringdesc;
        char* keyc = (char*)key.c_str(); // discarding constness
        CFI_establish(&stringdesc, keyc, CFI_attribute_other, CFI_type_char, strlen(keyc), 0, NULL);
        cset_option_float33(blob, &stringdesc, &x);
}
void cset_option(void* blob, std::string key, double x) {
        CFI_cdesc_t stringdesc;
        char* keyc = (char*)key.c_str(); // discarding constness
        CFI_establish(&stringdesc, keyc, CFI_attribute_other, CFI_type_char, strlen(keyc), 0, NULL);
        cset_option_float(blob, &stringdesc, x);
}
void cset_option(void* blob, std::string key, double* x, int i1, int i2) {
        CFI_cdesc_t stringdesc;
        char* keyc = (char*)key.c_str(); // discarding constness
        CFI_establish(&stringdesc, keyc, CFI_attribute_other, CFI_type_char, strlen(keyc), 0, NULL);
        cset_option_floatxy(blob, &stringdesc, x, i1, i2);
}
#endif
