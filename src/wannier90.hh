#ifndef __wannier90_hh__
#define __wannier90_hh__

#include <ISO_Fortran_binding.h>
#include <complex>
#include <cstring>

extern "C" {
void cchkpt(void*, CFI_cdesc_t*);

void cinput_reader_special(void*, CFI_cdesc_t*);
void cinput_reader(void*);
void cinput_setopt(void*, CFI_cdesc_t*, int);
void cset_option_float33(void*, CFI_cdesc_t*, void*);
void cset_option_float(void*, CFI_cdesc_t*, double);
void cset_option_floatxy(void*, CFI_cdesc_t*, void*, int, int);
void cset_option_int3(void*, CFI_cdesc_t*, void*);
void cset_option_int(void*, CFI_cdesc_t*, int);
void cset_option_intx(void*, CFI_cdesc_t*, int*, int);

void coverlaps(void*);
void cdisentangle(void*);
void cwannierise(void*);

void* w90_create();
void w90_delete(void*);

void cset_kpoint_distribution(void*, int*);
void cset_parallel_comms(void*, int);

void cset_a_matrix(void*, std::complex<double>*);
void cset_eigval(void*, double*);
void cset_m_matrix_local(void*, std::complex<double>*);
void cset_m_orig(void*, std::complex<double>*);
void cset_u_matrix(void*, std::complex<double>*);
void cset_u_opt(void*, std::complex<double>*);

void cget_nn(void*, int&);
void cget_nnkp(void*, int*);
void cget_centres(void*, void*);
void cget_spreads(void*, void*);
}

#ifdef MPI_VERSION
#include <mpi.h>
void cset_parallel_comms(void* blob, MPI_Comm comm) {
        int fcomm = MPI_Comm_c2f(comm);
        cset_parallel_comms(blob, fcomm); // translate to fortran integer and set communicator
}
void cinput_setopt(void* blob, std::string seed, MPI_Comm comm) {
        int fcomm = MPI_Comm_c2f(comm);
        CFI_cdesc_t stringdesc;
        char* seedc = (char*)seed.c_str(); // discarding constness
        CFI_establish(&stringdesc, seedc, CFI_attribute_other, CFI_type_char, strlen(seedc), 0, NULL);
        cinput_setopt(blob, &stringdesc, fcomm);
}
#else
void cinput_setopt(void* blob, std::string seed, int comm) {
        CFI_cdesc_t stringdesc;
        char* seedc = (char*)seed.c_str(); // discarding constness
        CFI_establish(&stringdesc, seedc, CFI_attribute_other, CFI_type_char, strlen(seedc), 0, NULL);
        cinput_setopt(blob, &stringdesc, 0);
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

void cinput_reader_special(void* blob, std::string seed) {
        CFI_cdesc_t stringdesc;
        char* seedc = (char*)seed.c_str(); // discarding constness
        CFI_establish(&stringdesc, seedc, CFI_attribute_other, CFI_type_char, strlen(seedc), 0, NULL);
        cinput_reader_special(blob, &stringdesc);
}
/*void* getglob(const std::string seed) {
        CFI_cdesc_t stringdesc;
        char* seedc = (char*)seed.c_str(); // discarding constness
        CFI_establish(&stringdesc, seedc, CFI_attribute_other, CFI_type_char, strlen(seedc), 0, NULL);
        return getglob(&stringdesc);
}*/
void cchkpt(void* blob, std::string text2) {
        CFI_cdesc_t desc2;
        char* text2c = (char*)text2.c_str(); // discarding constness
        CFI_establish(&desc2, text2c, CFI_attribute_other, CFI_type_char, strlen(text2c), 0, NULL);
        cchkpt(blob, &desc2);
}
#endif
