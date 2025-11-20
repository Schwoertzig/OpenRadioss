#ifndef VECTOR_DOUBLE_H
#define VECTOR_DOUBLE_H

#include <stdint.h>
#include "my_real.h"
#include "vector_int.h"

typedef struct {
    my_real* data;
    uint64_t size, capacity;
} Vector_double;

void push_back_vec_double(Vector_double* v, const my_real* point);
void double_capacity_vec_double(Vector_double* v);
Vector_double* alloc_empty_vec_double();
Vector_double* alloc_with_init_vec_double(const my_real* points, uint64_t size);
Vector_double* alloc_with_capacity_vec_double(uint64_t size);
void dealloc_vec_double(Vector_double* v);
my_real* get_ith_elem_vec_double(const Vector_double* v, uint64_t i);
void set_ith_elem_vec_double(Vector_double* v, uint64_t i, my_real* d);
void copy_vec_double(const Vector_double* src, Vector_double* dest);
void print_vec_double(const Vector_double* v);
void sort_vec_double(const Vector_double* src, Vector_double* sorted_array, Vector_uint* permutation_inds);



#endif