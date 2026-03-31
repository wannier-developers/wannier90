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

void w90_create(w90_data* w90_obj);
void w90_delete(w90_data* w90_obj);

void w90_set_option_double2d_f(w90_data w90_obj, CFI_cdesc_t* keyword, double* arg_cptr, int x, int y);
void w90_set_option_double2d(w90_data blob, const char* key, double* data, int x, int y);

void w90_set_option_double1d_f(w90_data w90_obj, CFI_cdesc_t* keyword, double* arg_cptr, int x);
void w90_set_option_double1d(w90_data blob, const char* key, double* data, int x);

void w90_set_option_double_f(w90_data w90_obj, CFI_cdesc_t* keyword, double value);
void w90_set_option_double(w90_data blob, const char* key, double data);

void w90_set_option_int2d_f(w90_data w90_obj, CFI_cdesc_t* keyword, int* arg_cptr, int x, int y);
void w90_set_option_int2d(w90_data blob, const char* key, int* data, int x, int y);

void w90_set_option_int1d_f(w90_data w90_obj, CFI_cdesc_t* keyword, int* arg_cptr, int x);
void w90_set_option_int1d(w90_data blob, const char* key, int* data, int x);

void w90_set_option_int_f(w90_data w90_obj, CFI_cdesc_t* keyword, int value);
void w90_set_option_int(w90_data blob, const char* key, int data);

void w90_set_option_logical_f(w90_data w90_obj, CFI_cdesc_t* keyword, bool value);
void w90_set_option_logical(w90_data blob, const char* key, bool data);

void w90_set_option_text_f(w90_data blob, CFI_cdesc_t* keyword, CFI_cdesc_t* text);
void w90_set_option_text(w90_data blob, const char* key, const char* data);

void w90_set_option_text2d_f(w90_data blob, CFI_cdesc_t* keyword, CFI_cdesc_t* text);
void w90_set_option_text2d(w90_data blob, const char* key, const char* const* data, int n);

void w90_input_setopt_f(w90_data w90_obj, CFI_cdesc_t* seedname, int* ierr);
void w90_input_setopt(w90_data blob, const char* seed, int* ierr);

void w90_input_reader(w90_data w90_obj, int* ierr);

void w90_get_proj(w90_data w90_obj, int* n, double* site, int* l, int* m, int* s, int* rad, double* x, double* z, double* sqa, double* zona, int* ierr);
void w90_get_nproj(w90_data w90_obj, int* n);
void w90_get_num_excl_bands(w90_data w90_obj, int* num_excl_bands);
void w90_get_excl_bands(w90_data w90_obj, int* excl_bands, int* ierr);

void w90_get_nk(w90_data w90_obj, int* n);
void w90_get_nw(w90_data w90_obj, int* n);
void w90_get_nn(w90_data w90_obj, int* n, int* ierr);
void w90_get_nnkp(w90_data w90_obj, int* nnkp, int* ierr);
void w90_get_gkpb(w90_data w90_obj, int* gkpb, int* ierr);

void w90_set_eigval(w90_data w90_obj, double* eigval);
void w90_set_m_local(w90_data w90_obj, _Complex double* m_matrix_local);
void w90_set_u_matrix(w90_data w90_obj, _Complex double* u_matrix);
void w90_set_u_opt(w90_data w90_obj, _Complex double* u_matrix_opt);

void w90_project_overlap(w90_data w90_obj, int* ierr);
void w90_print_info(w90_data w90_obj, int* ierr);
void w90_disentangle(w90_data w90_obj, int* ierr);
void w90_wannierise(w90_data w90_obj, int* ierr);

void w90_get_centres(w90_data w90_obj, double* centres);
void w90_get_spreads(w90_data w90_obj, double* spreads);

#ifdef MPI
void w90_set_comm(w90_data blob, MPI_Comm comm);
#endif
#endif
