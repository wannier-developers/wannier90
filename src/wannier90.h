#ifndef WANNIER90_H
#define WANNIER90_H

#include <stdbool.h>
#include <ISO_Fortran_binding.h>

/* CFI_cdesc_t provides a mechanism for passing variable length char* arrays */
// see https://community.intel.com/t5/Intel-Fortran-Compiler/C-interoperablilty-and-character-strings/td-p/1084167

/* the internals of the fortran data structure are not visible from C
 * provide a type wrapping a void* pointer to memory allocated by fortran
 */
typedef struct {
    void* cptr;
} w90_data;

/* Convenience constructor for creating and initializing w90_data */

void w90_create(w90_data*);
void w90_delete(w90_data*);

void w90_set_option_double2d_f(w90_data, CFI_cdesc_t*, void*, int, int);
void w90_set_option_double2d(w90_data blob, char* key, void* data, int x, int y);

void w90_set_option_double1d_f(w90_data, CFI_cdesc_t*, void*, int);
void w90_set_option_double1d(w90_data blob, char* key, void* data, int x);

void w90_set_option_double_f(w90_data, CFI_cdesc_t*, double);
void w90_set_option_double(w90_data blob, char* key, double data);

void w90_set_option_int2d_f(w90_data, CFI_cdesc_t*, void*, int, int);
void w90_set_option_int2d(w90_data blob, char* key, void* data, int x, int y);

void w90_set_option_int1d_f(w90_data, CFI_cdesc_t*, void*, int);
void w90_set_option_int1d(w90_data blob, char* key, void* data, int x);

void w90_set_option_int_f(w90_data, CFI_cdesc_t*, int);
void w90_set_option_int(w90_data blob, char* key, int data);

void w90_set_option_logical_f(w90_data, CFI_cdesc_t*, bool);
void w90_set_option_logical(w90_data blob, char* key, bool data);

void w90_set_option_text_f(w90_data blob, CFI_cdesc_t*, CFI_cdesc_t* data);
void w90_set_option_text(w90_data blob, char* key, char* data);

void w90_set_option_text2d_f(w90_data blob, CFI_cdesc_t*, CFI_cdesc_t* data);
void w90_set_option_text2d(w90_data blob, char* key, char** data, int n);

void w90_input_setopt_f(w90_data, CFI_cdesc_t*, int*);
void w90_input_setopt(w90_data blob, char* seed, int* ierr);

void w90_input_reader(w90_data, int*);

void w90_get_proj(w90_data, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, int*);

void w90_get_nk(w90_data, void*);
void w90_get_nw(w90_data, void*);
void w90_get_nn(w90_data, void*);
void w90_get_nnkp(w90_data, void*);
void w90_get_gkpb(w90_data, void*);

void w90_set_eigval(w90_data, void*);
void w90_set_m_local(w90_data, void*);
void w90_set_u_matrix(w90_data, void*);
void w90_set_u_opt(w90_data, void*);

void w90_project_overlap(w90_data, int*);
void w90_print_info(w90_data, int*);
void w90_disentangle(w90_data, int*);
void w90_wannierise(w90_data, int*);

void w90_get_centres(w90_data, void*);
void w90_get_spreads(w90_data, void*);

#ifdef MPI_VERSION
void w90_set_comm(w90_data blob, MPI_Comm comm);
#endif
#endif
