//
//  reads_stream.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef __XC_s2fastqIO__reads_stream__
#define __XC_s2fastqIO__reads_stream__

#include <stdio.h>
#include <ctype.h>
#include "stream_model.h"

#define MAX_ALPHA 5000000
#define MAX_CARDINALITY 50000000

typedef struct read_models_t{
    stream_model *flag;
    stream_model *pos;
    stream_model *match;
    stream_model *snps;
    stream_model *indels;
    stream_model *var;
    stream_model *chars;
    uint32_t read_length;
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



#endif /* defined(__XC_s2fastqIO__reads_stream__) */