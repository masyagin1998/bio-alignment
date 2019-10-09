#include "fasta.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

struct FASTA_DATA*fasta_data_create()
{
    return (struct FASTA_DATA*) calloc(1, sizeof(struct FASTA_DATA));
}
    
enum FASTA_CODES fasta_data_conf(struct FASTA_DATA*fdata, const char*comment, unsigned comment_len, const char*data, unsigned data_len)
{
    fdata->comment_len = comment_len;
    fdata->comment_cap = comment_len;
    fdata->comment = (char*) malloc(comment_len * sizeof(char));
    strncpy(fdata->comment, comment, comment_len);
    fdata->data_len = data_len;
    fdata->data_cap = data_len;
    fdata->data = (char*) malloc(data_len * sizeof(char));
    strncpy(fdata->data, data, data_len);

    return FASTA_OK;
}

void fasta_data_clear(struct FASTA_DATA*fdata)
{
    free(fdata->comment);
    free(fdata->data);
    memset(fdata, 0, sizeof(struct FASTA_DATA));
}

void fasta_data_free(struct FASTA_DATA*fdata)
{
    fasta_data_clear(fdata);
    free(fdata);
}

#define DATA_BUF 8

enum FASTA_READER_STATES
{
    FASTA_INV        = -1,
    FASTA_IN_COMMENT =  0,
    FASTA_IN_DATA    =  1,
};

static char*buf_push_back(char*buf, unsigned*buf_len, unsigned*buf_cap, char sym)
{
    if ((*buf_len) == (*buf_cap)) {
        (*buf_cap) *= 2;
        buf = (char*) realloc(buf, (*buf_cap) * sizeof(char));
    }

    buf[(*buf_len)] = sym;
    (*buf_len)++;

    return buf;
}

enum FASTA_CODES fasta_file_read(const char*fname, struct FASTA_DATA**fdata, unsigned*n)
{
    enum FASTA_CODES r;
    
    FILE*f;

    char*comment = (char*) malloc(DATA_BUF * sizeof(char));
    unsigned comment_len = 0;
    unsigned comment_cap = DATA_BUF;

    char*data = (char*) malloc(DATA_BUF * sizeof(char));
    unsigned data_len = 0;
    unsigned data_cap = DATA_BUF;    

    int read;

    enum FASTA_READER_STATES frs = FASTA_INV;

    char buf[DATA_BUF];
    
    f = fopen(fname, "r");
    if (f == NULL) {
        r = FASTA_BAD_FILE;
        goto err0;
    }

    (*n) = 0;

    while ((read = fread(buf, sizeof(char), sizeof(buf) / sizeof(char), f)) > 0) {
        int i;
        for (i = 0; i < read; i++) {
            if ((frs == FASTA_INV) && (buf[i] == '>')) {
                frs = FASTA_IN_COMMENT;
            } else if ((frs == FASTA_IN_COMMENT) && (buf[i] != '\n')) {
                comment = buf_push_back(comment, &comment_len, &comment_cap, buf[i]);
            } else if ((frs == FASTA_IN_COMMENT) && (buf[i] == '\n')) {
                frs = FASTA_IN_DATA;
            } else if ((frs == FASTA_IN_DATA) && (((buf[i] >= 'A') && (buf[i] <= 'Z')) || (buf[i] == '-'))) {
                data = buf_push_back(data, &data_len, &data_cap, buf[i]);
            } else if ((frs == FASTA_IN_DATA) && (buf[i] == '\n')) {
                continue;
            } else if ((frs == FASTA_IN_DATA) && (buf[i] == '>')) {
                if ((*n) == 0) {
                    (*fdata) = (struct FASTA_DATA*) malloc(sizeof(struct FASTA_DATA));
                } else {
                    (*fdata) = (struct FASTA_DATA*) realloc((*fdata), ((*n) + 1) * sizeof(struct FASTA_DATA));
                }
                fasta_data_conf((*fdata) + (*n), comment, comment_len, data, data_len);
                (*n)++;
                frs = FASTA_IN_COMMENT;
                comment_len = 0;
                data_len = 0;                
            } else {
                r = FASTA_BAD_DATA;
                goto err1;
            }
        }
    }

    if ((frs == FASTA_IN_DATA) && (data_len != 0)) {
        if ((*n) == 0) {
            (*fdata) = (struct FASTA_DATA*) malloc(sizeof(struct FASTA_DATA));
        } else {
            (*fdata) = (struct FASTA_DATA*) realloc((*fdata), ((*n) + 1) * sizeof(struct FASTA_DATA));
        }
        fasta_data_conf((*fdata) + (*n), comment, comment_len, data, data_len);
        (*n)++;        
    }

    free(comment);
    free(data);
    fclose(f);

    return FASTA_OK;

 err1:
    {
        unsigned i;
        for (i = 0; i < (*n); i++) {
            fasta_data_clear((*fdata) + i);
        }
        free((*fdata));
    }
 err0:
    (*fdata) = NULL;
    (*n) = 0;
    return r;
}

enum FASTA_CODES fasta_file_write(const char*fname, const struct FASTA_DATA*fdata, unsigned n, unsigned max_data_len)
{
    enum FASTA_CODES r;

    unsigned i, j;
    int f_need_close = 0;
    
    FILE*f;

    if (strcmp(fname, "stdout") == 0) {
        f = stdout;
    } else if (strcmp(fname, "stderr") == 0) {
        f = stderr;
    } else {
        f = fopen(fname, "w");
        if (f == NULL) {
            r = FASTA_BAD_FILE;
            goto err0;
        }
        f_need_close = 1;
    }

    for (i = 0; i < n; i++) {
        fprintf(f, ">%.*s\n", fdata[i].comment_len, fdata[i].comment);
        for (j = 0; j < fdata[i].data_len / max_data_len; j++) {
            fprintf(f, "%.*s\n", max_data_len, fdata[i].data + j * max_data_len);
        }
        if (fdata[i].data_len % max_data_len != 0) {
            int offset = (fdata[i].data_len / max_data_len) * max_data_len;
            fprintf(f, "%.*s\n", fdata[i].data_len - offset, fdata[i].data + offset);
        }
        fprintf(f, "\n");
    }

    if (f_need_close) {
        fclose(f);
    }

    return FASTA_OK;

 err0:
    return r;
}
