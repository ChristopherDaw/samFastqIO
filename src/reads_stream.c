//
//  reads_stream.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "reads_stream.h"


//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  INITIALIZATION                                      //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////

stream_stats *initialize_stream_stats_flag(uint32_t rescale){
    
    uint32_t i = 0;
    uint32_t context_size = 1;
    uint32_t counts_size = 8;
    
    stream_stats *s = (stream_stats*) calloc(context_size, sizeof(stream_stats));
    
    for (i = 0; i < context_size; i++) {
    
        s[i] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
        // Allocate memory
        s[i]->counts = (uint32_t*) calloc(counts_size, sizeof(uint32_t));
    
        // An extra for the cumcounts
        s[i]->counts += 1;
    
        s[i]->alphabetCard = 2;
    
        s[i]->counts[0] = 1;
        s[i]->counts[1] = 1;
    
        s[i]->n = 2;
    
        // STEP
        s[i]->step = 1;
    
        //rescale factor
        s[i]->rescale = rescale;
    }
    
    return s;
    
}

ppm0_stream_stats* initialize_stream_stats_pos(uint32_t rescale){
    
    ppm0_stream_stats *s;
    
    s = (ppm0_stream_stats*) calloc(1, sizeof(ppm0_stream_stats));
    
    s[0] = (ppm0_stream_stats) calloc(1, sizeof(struct ppm0_stream_stats_t));
    
    // Allocate memory
    s[0]->alphabet = (int32_t *) calloc(MAX_CARDINALITY + 2, sizeof(int32_t));
    s[0]->counts = (uint32_t*) calloc(MAX_CARDINALITY + 3, sizeof(uint32_t));
    s[0]->alphaExist = (uint8_t*) calloc(MAX_ALPHA + 2, sizeof(uint8_t));
    s[0]->alphaMap = (int32_t*) calloc(MAX_ALPHA + 2, sizeof(int32_t));
    
    // We shift the pointers to consider -1 (new alpha).
    s[0]->alphaMap += 1;
    s[0]->alphaExist += 1;
    s[0]->alphabet += 1;
    // An extra for the cumcounts
    s[0]->counts += 2;
    
    
    s[0]->alphabetCard = 0;
    s[0]->alphabet[-1] = -1;
    s[0]->alphaMap[-1] = -1;
    
    s[0]->alphaExist[-1] = 1;
    s[0]->counts[-1] = 1;
    s[0]->n = 1;
    
    // STEP
    s[0]->step = 10;
    
    s[0]->rescale = rescale;
    
    return s;
    
}

stream_stats* initialize_stream_stats_match(uint32_t rescale){
    
    uint32_t i = 0;
    uint32_t context_size = 256;
    
    stream_stats *s = (stream_stats*) calloc(context_size, sizeof(stream_stats));
    
    for (i = 0; i < context_size; i++) {
        
        s[i] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
        // Allocate memory
        s[i]->counts = (uint32_t*) calloc(8, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[i]->counts += 1;
        
        s[i]->alphabetCard = 2;
        
        s[i]->counts[0] = 1;
        s[i]->counts[1] = 1;
        
        s[i]->n = 2;
        
        // STEP
        s[i]->step = 1;
        
        //rescale bound
        s[i]->rescale = rescale;
    }
    
    
    
    
    
    return s;
    
}

stream_stats* initialize_stream_stats_S(uint32_t readLength, uint32_t rescale){
    
    stream_stats *s;
    
    uint32_t i = 0;
    
    s = (stream_stats*) calloc(1, sizeof(stream_stats));
    
    s[0] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
    // Allocate memory
    s[0]->counts = (uint32_t*) calloc(readLength + 2, sizeof(uint32_t));
    
    // An extra for the cumcounts
    s[0]->counts += 1;
    
    s[0]->alphabetCard = readLength;
    
    s[0]->n = 0;
    for (i = 0; i < readLength; i++) {
        s[0]->counts[i] = 1;
        s[0]->n += s[0]->counts[i];
    }
    
    // STEP
    s[0]->step = 10;
    
    //rescale bound
    s[0]->rescale = rescale;
    
    return s;
    
}

stream_stats* initialize_stream_stats_I(uint32_t readLength, uint32_t rescale){
    
    stream_stats *s;
    
    s = (stream_stats*) calloc(1, sizeof(stream_stats));
    
    uint32_t i = 0;
    
    s[0] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
    
    // Allocate memory
    s[0]->counts = (uint32_t*) calloc(readLength + 2, sizeof(uint32_t));
    
    // An extra for the cumcounts
    s[0]->counts += 1;
    
    s[0]->alphabetCard = readLength;
    
    s[0]->n = 0;
    for (i = 0; i < readLength; i++) {
        s[0]->counts[i] = 1;
        s[0]->n += s[0]->counts[i];
    }
    
    // STEP
    s[0]->step = 16;
    
    //rescale bound
    s[0]->rescale = rescale;
    
    return s;
    
}

stream_stats* initialize_stream_stats_var(uint32_t readLength, uint32_t rescale){
    
    stream_stats* s;
    
    uint32_t i = 0, j = 0;
    
    uint32_t num_models = 65536;
    
    s = (stream_stats*) calloc(num_models, sizeof(stream_stats));
    
    for (j = 0; j < num_models; j++) {
        
        s[j] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
        // Allocate memory
        s[j]->counts = (uint32_t*) calloc(readLength + 2, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[j]->counts += 1;
        
        s[j]->alphabetCard = readLength;
        
        s[j]->n = 0;
        for (i = 0; i < readLength; i++) {
            s[j]->counts[i] = 1;
            s[j]->n += s[0]->counts[i];
        }
        
        // STEP
        s[j]->step = 10;
        
        //rescale bound
        s[j]->rescale = rescale;
    }
    
    return s;
    
}

stream_stats* initialize_stream_stats_char(uint32_t rescale){
    
    stream_stats* s;
    
    s = (stream_stats*) calloc(6, sizeof(stream_stats));
    
    uint32_t i = 0, j;
    
    for (j = 0; j < 6; j++) {
        
        s[j] = (stream_stats) calloc(1, sizeof(struct stream_stats_t));
        
        // Allocate memory
        s[j]->counts = (uint32_t*) calloc(16, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[j]->counts += 1;
        
        s[j]->alphabetCard = 5;
        
        
        s[j]->n = 0;
        for (i = 0; i < 4; i++) {
            s[j]->counts[i] = (i == j)? 0:8;
            s[j]->n += s[j]->counts[i];
        }
        //{A,C,G,T} to N
        s[j]->counts[4] = 1;
        s[j]->n ++;
        
        // STEP
        s[j]->step = 8;
        
        //rescale bound
        s[j]->rescale = rescale;
    }
    
    s[0]->counts[1] += 8;
    s[0]->counts[2] += 8;
    s[0]->n += 16;
    
    s[1]->counts[0] += 8;
    s[1]->counts[3] += 8;
    s[1]->n += 16;
    
    s[2]->counts[0] += 8;
    s[2]->counts[3] += 8;
    s[2]->n += 16;
    
    s[3]->counts[1] += 8;
    s[3]->counts[2] += 8;
    s[3]->n += 16;
    
    // Allow N to N
    // s[4]->counts[4] = 1, s[4]->n++;
    
    // We need to increase the total count of the insertion case
    // s[5]->n++;
    
    return s;
    
}
