#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

static __inline__ int mat_get(const int*mat, unsigned c_len, unsigned i, unsigned j)
{
    return mat[i * c_len + j];
}

static __inline__ void mat_set(int*mat, unsigned c_len, unsigned i, unsigned j, int val)
{
    mat[i * c_len + j] = val;
}

static __inline__ int max(int a, int b) {
    return a > b ? a : b;
}

static __inline__ void mat_debug(const int*mat, unsigned r_len, unsigned c_len)
{
    unsigned i, j;
    for (i = 0; i < r_len; i++) {
        printf("|");
        for (j = 0; j < c_len - 1; j++) {
            printf("%d ", mat_get(mat, c_len, i, j));
        }
        printf("%d|\n", mat_get(mat, c_len, i, c_len - 1));
    }
}

#endif  /* UTILS_H */
