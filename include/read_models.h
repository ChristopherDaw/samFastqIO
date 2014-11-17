//
//  reads_stream.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef __XC_s2fastqIO__reads_stream__
#define __XC_s2fastqIO__reads_stream__

#include "stream_model.h"

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
read_models alloc_read_models_t(uint32_t read_length);

#endif /* defined(__XC_s2fastqIO__reads_stream__) */