#include "affine.h"

#include "utils.h"

#include <stdlib.h>
#include <limits.h>

#define DONE 0
#define LEFT 1
#define UP   2
#define DIAG 3

/* #define AFFINE_DEBUG */

#define MIN -1000000

void affine_run(const char*a, unsigned a_len, const char*b, unsigned b_len,
                char**a_aligned, char**b_aligned, unsigned*aligned_len,
                int*score, int (*scoring_function)(char a, char b), int gap, int gap_serial)
{
    unsigned i, j;

    int*G,*H,*V,*PTR;

    G = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));
    H = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));
    V = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));
    PTR = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));

    /* initialization of matrix. */
    for (i = 0; i < a_len + 1; i++) {
        mat_set(G, b_len + 1, i, 0, gap + ((i > 0) ? (i - 1) : 0) * gap_serial);
        mat_set(V, b_len + 1, i, 0, mat_get(G, b_len + 1, i, 0));
        mat_set(H, b_len + 1, i, 0, MIN);
        mat_set(PTR, b_len + 1, i, 0, UP);
    }
    for (i = 0; i < b_len + 1; i++) {
        mat_set(H, b_len + 1, 0, i, gap + ((i > 0) ? (i - 1) : 0) * gap_serial);
        mat_set(V, b_len + 1, 0, i, mat_get(H, b_len + 1, 0, i));
        mat_set(G, b_len + 1, 0, i, MIN);
        mat_set(PTR, b_len + 1, 0, i, LEFT);
    }

    mat_set(V, b_len + 1, 0, 0, 0);
    mat_set(PTR, b_len + 1, 0, 0, DONE);

    /* dynamic programming step. */
    for (i = 1; i < a_len + 1; i++) {
        for (j = 1; j < b_len + 1; j++) {
            int tmp;
            int max_val = INT_MIN;
            int diag;
            int from = DONE;

            tmp = max(mat_get(V, b_len + 1, i - 1, j) + gap, mat_get(G, b_len + 1, i - 1, j) + gap_serial);
            mat_set(G, b_len + 1, i, j, tmp);

            tmp = max(mat_get(V, b_len + 1, i, j - 1) + gap, mat_get(H, b_len + 1, i, j - 1) + gap_serial);
            mat_set(H, b_len + 1, i, j, tmp);

            diag = mat_get(V, b_len + 1, i - 1, j - 1) + scoring_function(a[i - 1], b[j - 1]);
            if (diag > max_val) {
                from = DIAG;
                max_val = diag;
            }

            tmp = mat_get(H, b_len + 1, i, j);
            if (tmp > max_val) {
                from = LEFT;
                max_val = tmp;
            }

            tmp = mat_get(G, b_len + 1, i, j);
            if (tmp > max_val) {
                from = UP;
                max_val = tmp;
            }

            mat_set(PTR, b_len + 1, i, j, from);
            mat_set(V, b_len + 1, i, j, max_val);
        }
    }

#ifdef AFFINE_DEBUG
    printf("G:\n");
    mat_debug(G, a_len + 1, b_len + 1);
    printf("H:\n");
    mat_debug(H, a_len + 1, b_len + 1);
    printf("V:\n");
    mat_debug(V, a_len + 1, b_len + 1);
    printf("PTR:\n");
    mat_debug(PTR, a_len + 1, b_len + 1);
#endif

    /* reverse. */
    (*aligned_len) = 0;
    (*a_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
    (*b_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));    
    
    i = a_len;
    j = b_len;

    while ((i != 0) || (j != 0)) {
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

    (*score) = mat_get(V, b_len + 1, a_len, b_len);

#ifdef AFFINE_DEBUG
    printf("%u\n", (*aligned_len));
    printf("a_aligned:\n%.*s\n", (*aligned_len), (*a_aligned));
    printf("b_aligned:\n%.*s\n", (*aligned_len), (*b_aligned));
    printf("score: %d\n", (*score));
#endif

    free(PTR);
    free(V);
    free(H);
    free(G);
}
