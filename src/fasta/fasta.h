#ifndef FASTA_H
#define FASTA_H

enum FASTA_CODES
{
    FASTA_OK = 0,
    FASTA_BAD_FILE = -1,
    FASTA_BAD_DATA = -2,
};

struct FASTA_DATA
{
    char*comment;
    unsigned comment_len;
    unsigned comment_cap;
    char*data;
    unsigned data_len;
    unsigned data_cap;
};

struct FASTA_DATA*fasta_data_create();
enum FASTA_CODES fasta_data_conf(struct FASTA_DATA*fdata, const char*comment, unsigned comment_len, const char*data, unsigned data_len);
void fasta_data_clear(struct FASTA_DATA*fdata);
void fasta_data_free(struct FASTA_DATA*fdata);

enum FASTA_CODES fasta_file_read(const char*fname, struct FASTA_DATA**fdata, unsigned*n);
enum FASTA_CODES fasta_file_write(const char*fname, const struct FASTA_DATA*fdata, unsigned n, unsigned max_data_len);

#endif  /* FASTA_H */
