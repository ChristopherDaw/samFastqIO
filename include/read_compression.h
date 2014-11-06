//
//  reads_compression.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef XC_s2fastqIO_reads_compression_h
#define XC_s2fastqIO_reads_compression_h

#include "read_models.h"
#include "Arithmetic_stream.h"
#include "sam_line.h"

#define MAX_BP_CHR 300000000


//#define MAX_PATH_LENGTH 256
//#define MAX_HEADER 512

#define BITS_DELTA 7

uint32_t READ_LENGTH;
char _readLength[4];

char *reference;

uint8_t snpInRef[MAX_BP_CHR] = {0};
uint32_t cumsumP = 0;


// Stuctures to handle the reads
typedef struct ch_t{
    char refChar;
    char targetChar;
}ch_t;

typedef struct snp{
    uint32_t pos;
    enum BASEPAIR refChar;
    enum BASEPAIR targetChar;
    uint32_t ctr;
}snp;

typedef struct ins{
    uint32_t pos;
    enum BASEPAIR targetChar;
}ins;


// Protorypes for the compression functions
uint32_t compress_flag(Arithmetic_stream a, stream_model *F, uint16_t flag);
uint32_t compress_pos(Arithmetic_stream as, stream_model *P, uint32_t pos);
uint32_t compress_match(Arithmetic_stream a, stream_model *M, uint8_t match, uint32_t P);
uint32_t compress_snps(Arithmetic_stream a, stream_model *S, uint8_t numSnps);
uint32_t compress_indels(Arithmetic_stream a, stream_model *I, uint8_t numIndels);
uint32_t compress_chars(Arithmetic_stream a, stream_model *c, enum BASEPAIR ref, enum BASEPAIR target);
uint32_t compress_var(Arithmetic_stream a, stream_model *v, uint32_t pos, uint32_t prevPos, uint32_t flag);

// Prototypes for the functions to extract the information from the reads
uint32_t compress_edits(Arithmetic_stream as, read_models rs, char *edits, char *cigar, char *read, uint32_t P, uint8_t flag);
int add_snps_to_array(char* edits, snp* SNPs, unsigned int *numSnps, unsigned int insertionPos, char *read);
uint32_t compute_delta_to_first_snp(uint32_t prevPos, uint32_t readLen);

#endif
