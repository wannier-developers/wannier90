#include "wannier90.h"
#include <string.h>

/* Function implementations */

void w90_set_option_double2d(w90_data blob, char* key, void* data, int x, int y) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_double2d_f(blob, &desc, data, x, y);
}

void w90_set_option_double1d(w90_data blob, char* key, void* data, int x) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_double1d_f(blob, &desc, data, x);
}

void w90_set_option_double(w90_data blob, char* key, double data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_double_f(blob, &desc, data);
}

void w90_set_option_int2d(w90_data blob, char* key, void* data, int x, int y) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_int2d_f(blob, &desc, data, x, y);
}

void w90_set_option_int1d(w90_data blob, char* key, void* data, int x) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_int1d_f(blob, &desc, data, x);
}

void w90_set_option_int(w90_data blob, char* key, int data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_int_f(blob, &desc, data);
}

void w90_set_option_logical(w90_data blob, char* key, bool data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_logical_f(blob, &desc, data);
}

void w90_set_option_text(w90_data blob, char* key, char* data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    CFI_cdesc_t desc1;
    CFI_establish(&desc1, data, CFI_attribute_other, CFI_type_char, strlen(data), 0, NULL);
    w90_set_option_text_f(blob, &desc, &desc1);
}

void w90_input_setopt(w90_data blob, char* seed, int* ierr) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, seed, CFI_attribute_other, CFI_type_char, strlen(seed), 0, NULL);
    *ierr = 0;
    w90_input_setopt_f(blob, &desc, ierr);
}

#ifdef MPI_VERSION
#include <mpi.h>
void w90_set_comm(w90_data blob, MPI_Comm comm) {
    int fcomm = MPI_Comm_c2f(comm);
    w90_set_comm_f(blob.cptr, fcomm);
}
#endif
