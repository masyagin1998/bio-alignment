#include "fasta.h"
#include "needleman-wunsch.h"
#include "smith-waterman.h"
#include "hirschberg.h"
#include "scoring-functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>

#define MAX_STR 1024

#define SCORING_DEFAULT  "default"
#define SCORING_BLOSUM62 "blosum62"
#define SCORING_DNAFULL  "dnafull"

#define ALGO_NEEDLEMAN_WUNSCH "nw"
#define ALGO_SMITH_WATERMAN   "sw"
#define ALGO_HIRSCHBERG       "hirsch"

static void print_version(const char*program_name)
{
    fprintf(stderr, "%s\n", program_name);
    fprintf(stderr, "NW V0.1\n");
    fprintf(stderr, "NW is simple terminal tool for\n");
    fprintf(stderr, "processing FASTA data with\n");
    fprintf(stderr, "Needleman-Wunsch algorithm.\n");
    exit(0);
}

static void print_help(const char*program_name)
{
    fprintf(stderr, "Usage example: %s -i data/in.fasta -o data/out.fasta -s blosum62 -g -5\n", program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  --version (-v)\n");
    fprintf(stderr, "            Print version info.\n");
    fprintf(stderr, "  --help    (-h)\n");
    fprintf(stderr, "            Print this help info.\n");
    fprintf(stderr, "  --in      (-i)=<file1,file2>\n");
    fprintf(stderr, "            One or two input files, separated by commas.\n");
    fprintf(stderr, "  --out     (-o)=<file>\n");
    fprintf(stderr, "            Output file (default: stdout).\n");
    fprintf(stderr, "  --algo    (-a)=<nw|sw|hirsch>\n");
    fprintf(stderr, "            Algorithm for processing sequences\n");
    fprintf(stderr, "            (nw - Needleman-Wunsch)\n");
    fprintf(stderr, "            (sw - Smith-Waterman)\n");
    fprintf(stderr, "            (h - Hirschberg)\n");
    fprintf(stderr, "            (default: nw)\n");
    fprintf(stderr, "  --gap     (-g)=<int>\n");
    fprintf(stderr, "            Gap score (default: -2).\n");
    fprintf(stderr, "  --scoring (-s)=<blosum62|dnafull|default>\n");
    fprintf(stderr, "            Default: +1 for match, -1 for mismatch.\n");
    exit(0);
}

char in_0[MAX_STR] = "";
char in_1[MAX_STR] = "";
char out[MAX_STR] = "stdout";
char algo[MAX_STR] = "nw";
int  gap = -2;
char scoring[MAX_STR] = "default";

int parse_args(int argc, char**argv)
{
    struct option opts[] = {
        {"version", 0, 0, 'v'},
        {"help",    0, 0, 'h'},
        {"in",      1, 0, 'i'},
        {"out",     1, 0, 'o'},
        {"algo",    1, 0, 'a'},
        {"gap",     1, 0, 'g'},
        {"scoring", 1, 0, 's'},
        {0,0,0,0}
    };

    int c;
    int idx;
    while ((c = getopt_long(argc, argv, "i:o:a:g:s:vh", opts, &idx)) != -1) {
        switch (c) {
        case 'v':
            print_version(argv[0]);
            break;
        case 'h':
            print_help(argv[0]);
            break;
        case 'i': {
            unsigned len = strlen(optarg);            
            unsigned i;
            for (i = 0; i < len; i++) {
                if (optarg[i] == ',') {
                    break;
                }
            }
            if (i == len) {
                strncpy(in_0, optarg, sizeof(in_0));
            } else {
                optarg[i] = '\0';
                strncpy(in_0, optarg, sizeof(in_0));
                strncpy(in_1, optarg + (i + 1), sizeof(in_1));
            }
            break;   
        }
        case 'o':
            strncpy(out, optarg, sizeof(out));
            break;
        case 'a':
            strncpy(algo, optarg, sizeof(algo));
            break;
        case 'g':
            gap = atoi(optarg);
            break;
        case 's':
            strncpy(scoring, optarg, sizeof(scoring));
            break;
        default:
            /* do nothing. */
            break;
        }
    }

    return 0;
}

int main(int argc, char**argv)
{
    int r;

    struct FASTA_DATA*fdata_in_0 = NULL;
    unsigned fdata_in_0_len;
    struct FASTA_DATA*fdata_in_1 = NULL;
    unsigned fdata_in_1_len;

    int (*scoring_function)(char a, char b);
    void (*algo_run)(const char*a, unsigned a_len, const char*b, unsigned b_len,
                     char**a_aligned, char**b_aligned, unsigned*aligned_len,
                     int*score, int (*scoring_function)(char a, char b), int G);
    
    char*a_aligned;
    char*b_aligned;
    unsigned aligned_len;
    int score;
    
    struct FASTA_DATA*fdata_out;

    unsigned i;
    
    parse_args(argc, argv);

    if ((in_0[0] != '\0') && (in_1[0] != '\0')) {
        if ((r = fasta_file_read(in_0, &fdata_in_0, &fdata_in_0_len)) != 0) {
            fprintf(stderr, "Invdalid input file \"%s\"! Error code: %d.\n", in_0, r);
            exit(1);
        }
        if (fdata_in_0_len != 1) {
            fprintf(stderr, "Invdalid input file \"%s\" length (%d)!\n", in_0, fdata_in_0_len);
            exit(1);            
        }
        if ((r = fasta_file_read(in_1, &fdata_in_1, &fdata_in_1_len)) != 0) {
            fprintf(stderr, "Invdalid input file \"%s\"! Error code: %d.\n", in_1, r);
            exit(1);            
        }
        if (fdata_in_1_len != 1) {
            fprintf(stderr, "Invdalid input file \"%s\" length (%d)!\n", in_1, fdata_in_1_len);
            exit(1);            
        }        
    } else if (in_0 != '\0') {
        if ((r = fasta_file_read(in_0, &fdata_in_0, &fdata_in_0_len)) != 0) {
            fprintf(stderr, "Invdalid input file \"%s\"! Error code: %d.\n", in_0, r);
            exit(1);
        }
        if (fdata_in_0_len != 2) {
            fprintf(stderr, "Invdalid input file \"%s\" length (%d)!\n", in_0, fdata_in_0_len);
            exit(1);            
        }        
        fdata_in_1 = fdata_in_0 + 1;
    } else {
        fprintf(stderr, "No input FASTA-file(s) is(are) specified!\n");
        exit(1);
    }

    if (strcmp(scoring, SCORING_DEFAULT) == 0) {
        scoring_function = &scoring_function_default;
    } else if (strcmp(scoring, SCORING_BLOSUM62) == 0) {
        scoring_function = &scoring_function_amino_acids_blosum62;
    } else if (strcmp(scoring, SCORING_DNAFULL) == 0) {
        scoring_function = &scoring_function_nucleotides_dna_full;
    } else {
        fprintf(stderr, "Invalid scoring function \"%s\"!\n", scoring);
        exit(1);        
    }

    if (strcmp(algo, ALGO_NEEDLEMAN_WUNSCH) == 0) {
        algo_run = &needleman_wunsch_run;
    } else if (strcmp(algo, ALGO_SMITH_WATERMAN) == 0) {
        algo_run = &smith_waterman_run;
    } else if (strcmp(algo, ALGO_HIRSCHBERG) == 0) {
        algo_run = &hirschberg_run;
    } else {
        fprintf(stderr, "Invalid algorithm function \"%s\"!\n", algo);
        exit(1);        
    }

    algo_run(fdata_in_0->data, fdata_in_0->data_len,
             fdata_in_1->data, fdata_in_1->data_len,
             &a_aligned, &b_aligned, &aligned_len,
             &score, scoring_function, gap);

    if ((in_0[0] != '\0') && (in_1[0] != '\0')) {
        fasta_data_free(fdata_in_0);
        fasta_data_free(fdata_in_1);
    } else if (in_0 != '\0') {
        for (i = 0; i < 2; i++) {
            fasta_data_clear(fdata_in_0 + i);
        }
        free(fdata_in_0);
    }

    fdata_out = (struct FASTA_DATA*) malloc(sizeof(struct FASTA_DATA) * 2);
    fasta_data_conf(fdata_out,     "", 0, a_aligned, aligned_len);
    fasta_data_conf(fdata_out + 1, "", 0, b_aligned, aligned_len);

    free(a_aligned);
    free(b_aligned);
    
    fasta_file_write(out, fdata_out, 2, 80);

    for (i = 0; i < 2; i++) {
        fasta_data_clear(fdata_out + i);
    }
    free(fdata_out);

    if (strcmp(out, "stdout") == 0) {
        fprintf(stdout, "\nScore: %d\n", score);
    } else if (strcmp(out, "stderr") == 0) {
        fprintf(stderr, "\nScore: %d\n", score);
    } else {
        FILE*f = fopen(out, "a");
        fprintf(f, "\nScore: %d\n", score);
        fclose(f);
    }

    return 0;
}
