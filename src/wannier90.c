#include "wannier90.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Function implementations */

void w90_set_option_double2d(w90_data blob, const char* key, double* data, int x, int y) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_double2d_f(blob, &desc, data, x, y);
}

void w90_set_option_double1d(w90_data blob, const char* key, double* data, int x) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_double1d_f(blob, &desc, data, x);
}

void w90_set_option_double(w90_data blob, const char* key, double data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_double_f(blob, &desc, data);
}

void w90_set_option_int2d(w90_data blob, const char* key, int* data, int x, int y) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_int2d_f(blob, &desc, data, x, y);
}

void w90_set_option_int1d(w90_data blob, const char* key, int* data, int x) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_int1d_f(blob, &desc, data, x);
}

void w90_set_option_int(w90_data blob, const char* key, int data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_int_f(blob, &desc, data);
}

void w90_set_option_logical(w90_data blob, const char* key, bool data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    w90_set_option_logical_f(blob, &desc, data);
}

void w90_set_option_text(w90_data blob, const char* key, const char* data) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);
    CFI_cdesc_t desc1;
    CFI_establish(&desc1, (char*)data, CFI_attribute_other, CFI_type_char, strlen(data), 0, NULL);
    w90_set_option_text_f(blob, &desc, &desc1);
}

void w90_set_option_text2d(w90_data blob, const char* key, const char* const* data, int n) {
    /* The data should be an array of strings, where each string is null-terminated.
    int n = 2; // 2 strings
    const char *data[2] = {"Ga", "N"}; // Array of strings
    */
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)key, CFI_attribute_other, CFI_type_char, strlen(key), 0, NULL);

    // Fortran expects each string to be fixed length, so we need to determine
    // the maximum string length and pad the strings accordingly.
    size_t max_len = 0;
    // 1. Compute the lengths and find the maximum
    for (int i = 0; i < n; i++) {
        size_t current_len = strlen(data[i]);
        if (current_len > max_len) {
            max_len = current_len;
        }
    }
    // 2. Create a flat buffer for Fortran
    // Fortran expects contiguous memory for an array of strings.
    char *flat_array = malloc(n * max_len);
    if (!flat_array) {
        fprintf(stderr, "Error: memory allocation failed in w90_set_option_text2d\n");
        return;
    }
    
    for (int i = 0; i < n; i++) {
        // Copy string and pad the rest with spaces (Fortran style)
        strncpy(&flat_array[i * max_len], data[i], max_len);
        // Fill remaining space with blanks so Fortran 'trim' works
        for (size_t j = strlen(data[i]); j < max_len; j++) {
            flat_array[i * max_len + j] = ' ';
        }
    }
    // 3. Establish the Descriptor
    CFI_CDESC_T(1) desc_data;
    CFI_index_t extent[] = {n}; // number of elements
    CFI_cdesc_t *desc_data_ptr = (CFI_cdesc_t *)&desc_data;
    // We use max_len as the element size
    CFI_establish(desc_data_ptr, flat_array, CFI_attribute_other, 
                  CFI_type_char, max_len, 1, extent);

    w90_set_option_text2d_f(blob, &desc, desc_data_ptr);

    free(flat_array);
}

void w90_input_setopt(w90_data blob, const char* seed, int* ierr) {
    CFI_cdesc_t desc;
    CFI_establish(&desc, (char*)seed, CFI_attribute_other, CFI_type_char, strlen(seed), 0, NULL);
    *ierr = 0;
    w90_input_setopt_f(blob, &desc, ierr);
}

#ifdef MPI
#include <mpi.h>
void w90_set_comm(w90_data blob, MPI_Comm comm) {
    int fcomm = MPI_Comm_c2f(comm);
    w90_set_comm_f(blob.cptr, fcomm);
}
#endif
