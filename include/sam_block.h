//
//  sam_line.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef XC_s2fastqIO_sam_line_h
#define XC_s2fastqIO_sam_line_h

#include <stdio.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "stream_model.h"
#include "pmf.h"
#include "qv_codebook.h"
#include "util.h"
#include "well.h"
#include "quantizer.h"


#define MAX_READ_LENGTH 1024

// This limits us to chunks that aren't too big to fit into a modest amount of memory at a time
#define MAX_LINES_PER_BLOCK			1000000
#define MAX_READS_PER_LINE			1022
#define READ_LINEBUF_LENGTH			(MAX_READS_PER_LINE+2)

// Error codes for reading a line block
#define LF_ERROR_NONE				0
#define LF_ERROR_NOT_FOUND			1
#define LF_ERROR_NO_MEMORY			2
#define LF_ERROR_TOO_LONG			4

#define get_qv_model_index(a,b) ((a & 0xff) << 8 | (b & 0xff))

#define MAX_ALPHA 5000000
#define MAX_CARDINALITY 50000000



typedef struct read_models_t{
    stream_model *flag;
    stream_model *pos;
    stream_model *pos_alpha;
    stream_model *match;
    stream_model *snps;
    stream_model *indels;
    stream_model *var;
    stream_model *chars;
    uint32_t read_length;
    char _readLength[4];
}*read_models;

// To store the model of the chars both in ref and target
enum BASEPAIR {
    A,
    C,
    G,
    T,
    N,
    O
};

/**
 *
 */
typedef struct read_line_t{
    char *cigar;
    char *edits;
    char *read;
    uint16_t invFlag;
    int pos;
    uint32_t read_length;
}*read_line;

/**
 *
 */
typedef struct read_block_t{
    struct read_line_t *lines;
    read_models models;
    uint32_t block_length;
} *read_block;


/**
 *
 */
typedef struct qv_line_t {
    symbol_t *data;
    uint32_t columns;
}*qv_line;


/**
 * Points to a file descriptor that includes important metadata about the qv
 */
typedef struct qv_block_t {
    struct alphabet_t *alphabet;
    uint32_t block_length;
    uint32_t columns;
    struct qv_line_t *qv_lines;
    struct cond_pmf_list_t *training_stats;
    struct cond_quantizer_list_t *qlist;
    struct distortion_t *dist;
    struct qv_options_t *opts;
    struct well_state_t well;
    stream_model *model;
    stream_model *codebook_model;
}*qv_block;

/**
 *
 */
typedef struct sam_block_t{
    qv_block QVs;
    read_block reads;
    FILE *fs;
    char *path;
    uint32_t read_length;
    uint32_t block_length;
    stream_model *codebook_model;
    FILE *fref;
    uint32_t current_line; // used for decompression
}*sam_block;


// Function Prototypes

int char2basepair(char c);
int basepair2char(enum BASEPAIR c);
char bp_complement(char c);

stream_model *initialize_stream_model_flag(uint32_t rescale);
stream_model* initialize_stream_model_pos(uint32_t rescale);
stream_model* initialize_stream_model_pos_alpha(uint32_t rescale);
stream_model* initialize_stream_model_match(uint32_t rescale);
stream_model* initialize_stream_model_snps(uint32_t readLength, uint32_t rescale);
stream_model* initialize_stream_model_indels(uint32_t readLength, uint32_t rescale);
stream_model* initialize_stream_model_var(uint32_t readLength, uint32_t rescale);
stream_model* initialize_stream_model_chars(uint32_t rescale);
stream_model* initialize_stream_model_qv(struct cond_quantizer_list_t *q_list, uint32_t rescale);
stream_model* initialize_stream_model_codebook(uint32_t rescale);

read_models alloc_read_models_t(uint32_t read_length);

void alloc_stream_model_qv(qv_block qvBlock);


sam_block alloc_sam_block_t(Arithmetic_stream as, FILE * fin, FILE *fref, struct qv_options_t *qv_opts, uint8_t decompression);
read_block alloc_read_block_t(uint32_t read_length);
qv_block alloc_qv_block_t(struct qv_options_t *opts, uint32_t read_length);
uint32_t get_read_length(FILE *f);
uint32_t load_sam_block(sam_block sb);

// Meat of the implementation
void calculate_statistics(struct qv_block_t *info);
double optimize_for_entropy(struct pmf_t *pmf, struct distortion_t *dist, double target, struct quantizer_t **lo, struct quantizer_t **hi, uint8_t verbose);
void generate_codebooks(struct qv_block_t *info);

// Master functions to handle codebooks in the output file
void write_codebooks(Arithmetic_stream as, struct qv_block_t *info);
void write_codebook(Arithmetic_stream as, struct cond_quantizer_list_t *quantizers, stream_model *model);
void read_codebooks(Arithmetic_stream as, struct qv_block_t *info);
struct cond_quantizer_list_t *read_codebook(Arithmetic_stream as, struct qv_block_t *info);

void initialize_qv_model(Arithmetic_stream as, qv_block qvBlock, uint8_t decompression);
void reset_QV_block(qv_block qvb, uint8_t direction);

stream_model *free_stream_model_qv(struct cond_quantizer_list_t *q_list, stream_model *s);


#endif
