#include "fasta.h"
#include "needleman-wunsch.h"

#include <stdlib.h>

int main(int argc, char**argv)
{
    struct FASTA_DATA*fdata;
    unsigned fdata_len;
    unsigned i;

    fasta_file_read(argv[1], &fdata, &fdata_len);

    char*a_aligned;
    char*b_aligned;
    unsigned aligned_len;
    int score;
    needleman_wunsch_run("CAT", 3, "KATE", 4,
                         &a_aligned, &b_aligned, &aligned_len,
                         &score, scoring_function_default, -10);
    
    // fasta_file_write("stdout", fdata, fdata_len, 60);

    free(a_aligned);
    free(b_aligned);

    for (i = 0; i < fdata_len; i++) {
        fasta_data_clear(fdata + i);
    }
    
    free(fdata);
}
