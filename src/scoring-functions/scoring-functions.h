#ifndef SCORING_FUNCTIONS_WUNSCH_H
#define SCORING_FUNCTIONS_WUNSCH_H

int scoring_function_amino_acids_blosum62(char a, char b);
int scoring_function_nucleotides_dna_full(char a, char b);
int scoring_function_default(char a, char b);

#endif  /* SCORING_FUNCTIONS_WUNSCH_H */
