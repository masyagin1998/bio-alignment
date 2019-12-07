#ifndef AFFINE_H
#define AFFINE_H

void affine_run(const char*a, unsigned a_len, const char*b, unsigned b_len,
                char**a_aligned, char**b_aligned, unsigned*aligned_len,
                int*score, int (*scoring_function)(char a, char b), int G, int G_serial);

#endif  /* AFFINE_H */
