#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

enum NEEDLEMAN_WUNSCH_CODES
{
    NEEDLEMAN_WUNSCH_Ok = 0,
    NEEDLEMAN_WUNSCH_BAD_ALLOC = -1,
};

enum NEEDLEMAN_WUNSCH_CODES needleman_wunsch_run(const char*a, unsigned a_len, const char*b, unsigned b_len,
                                                 char**a_aligned, char**b_aligned, unsigned*aligned_len,
                                                 int*score, int (*scoring_function)(char a, char b), int G);

int scoring_function_amino_acids_blosum62(char a, char b);
int scoring_function_nucleotides_dna_full(char a, char b);
int scoring_function_default(char a, char b);

#endif  /* NEEDLEMAN_WUNSCH_H */
