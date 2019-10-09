#include "needleman-wunsch.h"

#include <stdio.h>
#include <stdlib.h>

#define BLOSUM62_SIZE 24

const int BLOSUM62[BLOSUM62_SIZE][BLOSUM62_SIZE] =
{
        /* A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
 /* A */ { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
 /* R */ {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
 /* N */ {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
 /* D */ {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
 /* C */ { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
 /* Q */ {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
 /* E */ {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
 /* G */ { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
 /* H */ {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
 /* I */ {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
 /* L */ {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
 /* K */ {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
 /* M */ {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
 /* F */ {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
 /* P */ {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
 /* S */ { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
 /* T */ { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
 /* W */ {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
 /* Y */ {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
 /* V */ { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
 /* B */ {-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
 /* Z */ {-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
 /* X */ { 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
 /* * */ {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};

static int idx_blosum62(char a)
{
    switch (a) {
    case 'A':
        return 0;
    case 'R':
        return 1;
    case 'N':
        return 2;
    case 'D':
        return 3;
    case 'C':
        return 4;
    case 'Q':
        return 5;
    case 'E':
        return 6;
    case 'G':
        return 7;
    case 'H':
        return 8;
    case 'I':
        return 9;
    case 'L':
        return 10;
    case 'K':
        return 11;
    case 'M':
        return 12;
    case 'F':
        return 13;
    case 'P':
        return 14;
    case 'S':
        return 15;
    case 'T':
        return 16;
    case 'W':
        return 17;
    case 'Y':
        return 18;
    case 'V':
        return 19;
    case 'B':
        return 20;
    case 'Z':
        return 21;
    case 'X':
        return 22;
    default:
        return 23;
    }
}

int scoring_function_amino_acids_blosum62(char a, char b)
{
    return BLOSUM62[idx_blosum62(a)][idx_blosum62(b)];
}

#define DNA_FULL_SIZE 4

const int DNA_FULL[DNA_FULL_SIZE][DNA_FULL_SIZE] =
{
 /*        A   T   G   C */
 /* A */ { 5, -4, -4, -4},
 /* T */ {-4,  5, -4, -4},
 /* G */ {-4, -4,  5, -4},
 /* C */ {-4, -4, -4,  5},
};

static int idx_dna_full(char a)
{
    switch (a) {
    case 'A':
        return 0;
    case 'T':
        return 1;
    case 'G':
        return 2;
    case 'C':
        return 3;
    default:
        fprintf(stderr, "unknown character for DNA_FULL: %c", a);
        exit(1);
    }
}

int scoring_function_nucleotides_dna_full(char a, char b)
{
    return DNA_FULL[idx_dna_full(a)][idx_dna_full(b)];
}

int scoring_function_default(char a, char b)
{
    if (a == b) {
        return 1;
    } else {
        return -1;
    }
}
