#include "smith-waterman.h"

#include "utils.h"

#include <stdlib.h>
#include <limits.h>

#define DONE 0
#define LEFT 1
#define UP   2
#define DIAG 3

#define SMITH_WATERMAN_DEBUG

void smith_waterman_run(const char*a, unsigned a_len, const char*b, unsigned b_len,
                        char**a_aligned, char**b_aligned, unsigned*aligned_len,
                        int*score, int (*scoring_function)(char a, char b), int G)
{
    unsigned i, j;
    int*D,*PTR;

    unsigned i_1 = 0, j_1 = 0;
    int max_val = INT_MIN;
    
    D = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));
    PTR = (int*) calloc((a_len + 1) * (b_len + 1), sizeof(int));

    /* initialization of matrix. */
    mat_set(D, b_len + 1, 0, 0, 0);
    mat_set(PTR, b_len + 1, 0, 0, DONE);
    for (i = 1; i < (a_len + 1); i++) {
        mat_set(D, b_len + 1, i, 0, 0);
        mat_set(PTR, b_len + 1, i, 0, UP);
    }
    for (i = 1; i < (b_len + 1); i++) {
        mat_set(D, b_len + 1, 0, i, 0);
        mat_set(PTR, b_len + 1, 0, i, LEFT);
    }

    /* dynamic programming step. */
    for (i = 1; i < a_len + 1; i++) {
        for (j = 1; j < b_len + 1; j++) {
            int score_up   = mat_get(D, b_len + 1, i - 1, j) + G;
            int score_left = mat_get(D, b_len + 1, i, j - 1) + G;
            int score_diag = mat_get(D, b_len + 1, i - 1, j - 1) + scoring_function(a[i - 1], b[j - 1]);
            int score_max = max(0, max(score_up, max(score_left, score_diag)));
            mat_set(D, b_len + 1, i, j, score_max);
            if (score_max == score_up) {
                mat_set(PTR, b_len + 1, i, j, UP);
            } else if (score_max == score_left) {
                mat_set(PTR, b_len + 1, i, j, LEFT);
            } else {
                mat_set(PTR, b_len + 1, i, j, DIAG);
            }

            if (score_max > max_val) {
                printf("\n%d\n", score_max);
                max_val = score_max;
                i_1 = i;
                j_1 = j;
            }
        }
    }

#ifdef SMITH_WATERMAN_DEBUG
    printf("D:\n");
    mat_debug(D, a_len + 1, b_len + 1);
    printf("PTR:\n");
    mat_debug(PTR, a_len + 1, b_len + 1);
#endif

    /* reverse. */
    (*aligned_len) = 0;
    (*a_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
    (*b_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));    
    
    i = i_1;
    j = j_1;

    while ((i != 0) && (j != 0) && (mat_get(D, b_len + 1, i, j) != 0)) {
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

    (*score) = mat_get(D, b_len + 1, i_1, j_1);

#ifdef SMITH_WATERMAN_DEBUG
    printf("a_aligned:\n%.*s\n", (*aligned_len), (*a_aligned));
    printf("b_aligned:\n%.*s\n", (*aligned_len), (*b_aligned));
    printf("score: %d\n", (*score));
#endif

    free(PTR);
    free(D);
}
