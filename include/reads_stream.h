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
#include "Arithmetic_stream.h"

#define MAX_ALPHA 5000000
#define MAX_CARDINALITY 50000000

#endif /* defined(__XC_s2fastqIO__reads_stream__) */

typedef struct reads_stream_t{
    stream_stats flag;
    ppm0_stream_stats pos;
    stream_stats match;
    stream_stats snps;
    stream_stats indels;
    stream_stats var;
    stream_stats chars;
    uint32_t read_length;
}*reads_stream;