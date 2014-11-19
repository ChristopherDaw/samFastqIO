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

#include "pmf.h"
#include "codebook.h"
#include "qv_compressor.h"
#include "util.h"

#define MAX_READ_LENGTH 1024
#define MAX_READS_PER_BLOCK 1000000

// This limits us to chunks that aren't too big to fit into a modest amount of memory at a time
#define MAX_LINES_PER_BLOCK			1000000
#define MAX_READS_PER_LINE			1022
#define READ_LINEBUF_LENGTH			(MAX_READS_PER_LINE+2)

// Error codes for reading a line block
#define LF_ERROR_NONE				0
#define LF_ERROR_NOT_FOUND			1
#define LF_ERROR_NO_MEMORY			2
#define LF_ERROR_TOO_LONG			4

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
    uint64_t block_length;
    uint32_t columns;
    struct qv_line_t *lines;
    struct distortion_t *dist;
    struct qv_options_t *opts;
    struct well_state_t well;
}*qv_block;

/**
 *
 */
typedef struct sam_file_t{
    qv_block QVs;
    read_block reads;
    FILE *fs;
    char *path;
    uint32_t read_length;
}*sam_file;


// Function Prototypes
sam_line alloc_sam_line_t();
read_file alloc_read_file_t();

uint8_t read_line_from_sam(struct sam_line_t *samLine, FILE *f);

#endif
