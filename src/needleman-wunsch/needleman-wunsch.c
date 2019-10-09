#include "needleman-wunsch.h"

#include <stdlib.h>

#define DONE 0
#define LEFT 1
#define UP   2
#define DIAG 3

static inline int mat_get(const int*mat, unsigned c_len, unsigned i, unsigned j)
{
    return mat[i * c_len + j];
}

static inline void mat_set(int*mat, unsigned c_len, unsigned i, unsigned j, int val)
{
    mat[i * c_len + j] = val;
}

static inline int max(int a, int b) {
    return a > b ? a : b;
}

#define NEEDLEMAN_WUNSCH_DEBUG

#ifdef NEEDLEMAN_WUNSCH_DEBUG

#include <stdio.h>

void mat_debug(const int*mat, unsigned r_len, unsigned c_len)
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
#endif  /* NEEDLEMAN_WUNSCH_DEBUG */

enum NEEDLEMAN_WUNSCH_CODES needleman_wunsch_run(const char*a, unsigned a_len, const char*b, unsigned b_len,
                                                 char**a_aligned, char**b_aligned, unsigned*aligned_len,
                                                 int*score, int (*scoring_function)(char a, char b), int G)
{
    enum NEEDLEMAN_WUNSCH_CODES r;

    unsigned i, j;
    int*D,*PTR;
    
    D = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));
    if (D == NULL) {
        r = NEEDLEMAN_WUNSCH_BAD_ALLOC;
        goto err0;
    }

    PTR = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));
    if (PTR == NULL) {
        r = NEEDLEMAN_WUNSCH_BAD_ALLOC;
        goto err1;
    }

    /* initialization of matrix. */
    mat_set(D, b_len + 1, 0, 0, 0);
    mat_set(PTR, b_len + 1, 0, 0, DONE);
    for (i = 1; i < (a_len + 1); i++) {
        mat_set(D, b_len + 1, i, 0, i * G);
        mat_set(PTR, b_len + 1, i, 0, UP);
    }
    for (i = 1; i < (b_len + 1); i++) {
        mat_set(D, b_len + 1, 0, i, i * G);
        mat_set(PTR, b_len + 1, 0, i, LEFT);
    }

    /* dynamic programming step. */
    for (i = 1; i < a_len + 1; i++) {
        for (j = 1; j < b_len + 1; j++) {
            int score_up   = mat_get(D, b_len + 1, i - 1, j) + G;
            int score_left = mat_get(D, b_len + 1, i, j - 1) + G;
            int score_diag = mat_get(D, b_len + 1, i - 1, j - 1) + scoring_function(a[i - 1], b[j - 1]);
            int score_max = max(score_up, max(score_left, score_diag));
            mat_set(D, b_len + 1, i, j, score_max);
            if (score_max == score_up) {
                mat_set(PTR, b_len + 1, i, j, UP);
            } else if (score_max == score_left) {
                mat_set(PTR, b_len + 1, i, j, LEFT);
            } else {
                mat_set(PTR, b_len + 1, i, j, DIAG);
            }
        }
    }

#ifdef NEEDLEMAN_WUNSCH_DEBUG
    printf("D:\n");
    mat_debug(D, a_len + 1, b_len + 1);
    printf("PTR:\n");
    mat_debug(PTR, a_len + 1, b_len + 1);
#endif

    /* reverse. */
    (*aligned_len) = 0;
    (*a_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
    (*b_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));    
    
    i = a_len;
    j = b_len;

    while ((i != 0) && (j != 0)) {
        printf("%d\n", mat_get(PTR, b_len + 1, i, j));
        if (mat_get(PTR, b_len + 1, i, j) == DIAG) {
            (*a_aligned)[(*aligned_len)] = a[i - 1];
            (*b_aligned)[(*aligned_len)] = b[j - 1];
            (*aligned_len)++;
            i--;
            j--;
        } else if (mat_get(PTR, b_len + 1, i, j) == LEFT) {
            (*a_aligned)[(*aligned_len)] = '-';
            (*b_aligned)[(*aligned_len)] = b[j - 1];
            (*aligned_len)++;
            j--;
        } else if (mat_get(PTR, b_len + 1, i, j) == UP) {
            (*a_aligned)[(*aligned_len)] = a[i - 1];
            (*b_aligned)[(*aligned_len)] = '-';
            (*aligned_len)++;
            i--;
        }
    }

    for (i = 0; i < (*aligned_len) / 2; i++) {
        char tmp = (*a_aligned)[i];
        (*a_aligned)[i] = (*a_aligned)[(*aligned_len) - 1 - i];
        (*a_aligned)[(*aligned_len) - 1 - i] = tmp;
        tmp = (*b_aligned)[i];
        (*b_aligned)[i] = (*b_aligned)[(*aligned_len) - 1 - i];
        (*b_aligned)[(*aligned_len) - 1 - i] = tmp;
    }

    (*score) = mat_get(D, b_len + 1, a_len, b_len);

#ifdef NEEDLEMAN_WUNSCH_DEBUG
    printf("a_aligned:\n%.*s\n", (*aligned_len), (*a_aligned));
    printf("b_aligned:\n%.*s\n", (*aligned_len), (*b_aligned));
    printf("score: %d\n", (*score));
#endif    

    free(PTR);
    free(D);

    return NEEDLEMAN_WUNSCH_Ok;

 err1:
    free(D);
 err0:
    return r;
}
