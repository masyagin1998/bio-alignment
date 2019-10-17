#include "hirschberg.h"
#include "needleman-wunsch.h"

#include "utils.h"

#include <stdlib.h>
#include <string.h>

/* #define HIRSCHBERG_DEBUG */

#ifdef HIRSCHBERG_DEBUG

#include <stdio.h>

void print_arr_char(const char*arr, unsigned arr_len)
{
    printf("%.*s\n", arr_len, arr);
}

void print_arr_int(const int*arr, unsigned arr_len)
{
    unsigned i;
    for (i = 0; i < arr_len; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

#endif  /* HIRSCHBERG_DEBUG */

static void calculate_score(const char*a, unsigned a_len, const char*b, unsigned b_len,
							int**score, unsigned*score_len,
                            int (*scoring_function)(char a, char b), int G)
{
    unsigned i, j, k;
    unsigned tmp_len = 0;
    int tmp[2][b_len + 1];

    tmp[0][tmp_len] = 0;
    tmp[1][tmp_len] = 0;
    tmp_len++;

    for (j = 1; j < b_len + 1; j++) {
        tmp[0][tmp_len] = j * G;
        tmp[1][tmp_len] = 0;
        tmp_len++;
    }

    for (i = 1; i < a_len + 1; i++) {
        tmp[1][0] = tmp[0][0] + G;
        for (j = 1; j < b_len + 1; j++) {
            tmp[1][j] = max(max(tmp[0][j - 1] + scoring_function(a[i - 1], b[j - 1]), tmp[1][j - 1] + G), tmp[0][j] + G);
        }
        for (k = 0; k < tmp_len; k++) {
            tmp[0][k] = tmp[1][k];
        }
    }

    (*score_len) = tmp_len;
    (*score) = (int*) malloc((*score_len) * sizeof(int));
    for (j = 0; j < b_len + 1; j++) {
        (*score)[j] = tmp[1][j];
    }

#ifdef HIRSCHBERG_DEBUG
        printf("l:\n");
        print_arr_int((*score), (*score_len));
#endif  /* HIRSCHBERG_DEBUG */            
}

static void reverse_arr_int(int*arr, unsigned len)
{
    unsigned i;
    for (i = 0; i < (len / 2); i++) {
        int tmp = arr[i];
        arr[i] = arr[len - 1 - i];
        arr[len - 1 - i] = tmp;
    }
}

static void reverse_arr_char(char*arr, unsigned len)
{
    unsigned i;
    for (i = 0; i < (len / 2); i++) {
        char tmp = arr[i];
        arr[i] = arr[len - 1 - i];
        arr[len - 1 - i] = tmp;
    }
}

void hirschberg_run_inner(char*a, unsigned a_len, char*b, unsigned b_len,
                          char**a_aligned, char**b_aligned, unsigned*aligned_len,
                          int*score, int (*scoring_function)(char a, char b), int G)
{
    unsigned i;

#ifdef HIRSCHBERG_DEBUG
    printf("a:\n");
    print_arr_char(a, a_len);
    printf("b:\n");
    print_arr_char(b, b_len);    
#endif  /* HIRSCHBERG_DEBUG */

    if (a_len == 0) {
        (*aligned_len) = 0;
        (*a_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
        (*b_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
        
        for (i = 1; i < b_len + 1; i++) {
            (*a_aligned)[i - 1] = '-';
            (*b_aligned)[i - 1] = b[i - 1];
            (*aligned_len)++;
        }
    } else if (b_len == 0) {
        (*aligned_len) = 0;
        (*a_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
        (*b_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
        
        for (i = 1; i < a_len + 1; i++) {
            (*a_aligned)[i - 1] = a[i - 1];
            (*b_aligned)[i - 1] = '-';
            (*aligned_len)++;
        }
    } else if ((a_len == 1) || (b_len == 1)) {
        needleman_wunsch_run(a, a_len, b, b_len, a_aligned, b_aligned, aligned_len, score, scoring_function, G);
    } else {
		int*score1;
        unsigned score1_len;
        int*score2;
        unsigned score2_len;

		unsigned scores_sum_len;

		unsigned mid;
		int m;

        char*a_aligned_part_first,*b_aligned_part_first;
        unsigned a_aligned_part_first_len;
        char*a_aligned_part_second,*b_aligned_part_second;
        unsigned a_aligned_part_second_len;		
		
        (*aligned_len) = 0;
        (*a_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
        (*b_aligned) = (char*) malloc((a_len + b_len) * sizeof(char));
        
        calculate_score(a, a_len / 2, b, b_len, &score1, &score1_len, scoring_function, G);

        reverse_arr_char(a + a_len / 2, a_len - a_len / 2);
        reverse_arr_char(b, b_len);
        calculate_score(a + (a_len / 2), a_len - a_len / 2, b, b_len, &score2, &score2_len, scoring_function, G);
        reverse_arr_char(a + a_len / 2, a_len - a_len / 2);
        reverse_arr_char(b, b_len);

        reverse_arr_int(score2, score2_len);

        int scores_sum[score1_len];
		scores_sum_len = score1_len;

        for (i = 0; i < score1_len; i++) {
            scores_sum[i] = score1[i] + score2[i];
        }

        free(score1);
        free(score2);        

#ifdef HIRSCHBERG_DEBUG
        printf("scores_sum:\n");
        print_arr_int(scores_sum, scores_sum_len);
#endif  /* HIRSCHBERG_DEBUG */        

		mid = 0;
		m = scores_sum[0];

        for (i = 1; i < scores_sum_len; i++) {
            if (scores_sum[i] > m) {
                m = scores_sum[i];
                mid = i;
            }
        }

#ifdef HIRSCHBERG_DEBUG
        printf("mid:\n");
        printf("%u\n", mid);
#endif  /* HIRSCHBERG_DEBUG */

        hirschberg_run_inner(a, a_len / 2, b, mid, &a_aligned_part_first, &b_aligned_part_first, &a_aligned_part_first_len, score, scoring_function, G);
        hirschberg_run_inner(a + a_len / 2, a_len - a_len / 2, b + mid, b_len - mid, &a_aligned_part_second, &b_aligned_part_second, &a_aligned_part_second_len, score, scoring_function, G);

        (*aligned_len) = a_aligned_part_first_len + a_aligned_part_second_len;
        for (i = 0; i < (*aligned_len); i++) {
            if (i < a_aligned_part_first_len) {
                (*a_aligned)[i] = a_aligned_part_first[i];
                (*b_aligned)[i] = b_aligned_part_first[i];
            } else {
                (*a_aligned)[i] = a_aligned_part_second[i - a_aligned_part_first_len];
                (*b_aligned)[i] = b_aligned_part_second[i - a_aligned_part_first_len];
            }
        }

        free(a_aligned_part_first);
        free(b_aligned_part_first);
        free(a_aligned_part_second);
        free(b_aligned_part_second);
    }
}

void hirschberg_run(const char*a, unsigned a_len, const char*b, unsigned b_len,
                    char**a_aligned, char**b_aligned, unsigned*aligned_len,
                    int*score, int (*scoring_function)(char a, char b), int G)
{
	unsigned i;
	
	/* because of const pointer qualifiers. */
    char a1[a_len];
    char b1[b_len];
    memcpy(a1, a, a_len * sizeof(char));
    memcpy(b1, b, b_len * sizeof(char));

    hirschberg_run_inner(a1, a_len, b1, b_len, a_aligned, b_aligned, aligned_len, score, scoring_function, G);

	(*score) = 0;
	for (i = 0; i < (*aligned_len); i++) {
		if (((*a_aligned)[i] == '-') || ((*b_aligned)[i] == '-')) {
			(*score) += G;
		} else {
			(*score) += scoring_function((*a_aligned)[i], (*b_aligned)[i]);
		}
	}
}
