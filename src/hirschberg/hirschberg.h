#ifndef HIRSCHBERG_H
#define HIRSCHBERG_H

void hirschberg_run(const char*a, unsigned a_len, const char*b, unsigned b_len,
                    char**a_aligned, char**b_aligned, unsigned*aligned_len,
                    int*score, int (*scoring_function)(char a, char b), int G);

#endif  /* HIRSCHBERG_H */
