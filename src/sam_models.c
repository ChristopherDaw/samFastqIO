//
//  reads_stream.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "sam_block.h"


int char2basepair(char c)
{
    switch(c)
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

int basepair2char(enum BASEPAIR c)
{
    switch(c)
    {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

char bp_complement(char c){
    
    switch(c)
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return c;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////
//                                                                                      //
//                                                                                      //
//                                  INITIALIZATION                                      //
//                                                                                      //
//                                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////


stream_model* initialize_stream_model_id(uint32_t rescale, uint32_t context_size, uint32_t alphabet_card){
    
    uint32_t i = 0, j = 0;
    
    stream_model *s = (stream_model*) calloc(context_size, sizeof(stream_model));
    
    for (i = 0; i < context_size; i++) {
        
        s[i] = (stream_model) calloc(1, sizeof(struct stream_model_t));
        
        // Allocate memory
        s[i]->counts = (uint32_t*) calloc(alphabet_card + 1, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[i]->counts += 1;
        
        s[i]->alphabetCard = alphabet_card;
        
        s[i]->n = 0;
        for (j = 0; j < alphabet_card; j++) {
            s[i]->counts[j] = 1;
            s[i]->n += s[0]->counts[j];
        }
        
        // STEP
        s[i]->step = 10;
        
        //rescale bound
        s[i]->rescale = rescale;
    }
    
    
    
    
    
    return s;
    
}


stream_model *initialize_stream_model_flag(uint32_t rescale){
    
    uint32_t i = 0;
    uint32_t context_size = 1;
    uint32_t counts_size = 8;
    
    stream_model *s = (stream_model*) calloc(context_size, sizeof(stream_model));
    
    for (i = 0; i < context_size; i++) {
    
        s[i] = (stream_model) calloc(1, sizeof(struct stream_model_t));
    
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

stream_model* initialize_stream_model_pos(uint32_t rescale){
    
    stream_model *s;
    
    s = (stream_model*) calloc(1, sizeof(stream_model));
    
    s[0] = (stream_model) calloc(1, sizeof(struct stream_model_t));
    
    // Allocate memory
    s[0]->alphabet = (int32_t *) calloc(MAX_CARDINALITY, sizeof(int32_t));
    s[0]->counts = (uint32_t*) calloc(MAX_CARDINALITY, sizeof(uint32_t));
    s[0]->alphaExist = (uint8_t*) calloc(MAX_ALPHA, sizeof(uint8_t));
    s[0]->alphaMap = (int32_t*) calloc(MAX_ALPHA, sizeof(int32_t));
    
    // Initialize the alphabet assigning -1 to the bin 0 of the model
    s[0]->alphabetCard = 1;
    s[0]->alphabet[0] = -1;
    s[0]->alphaMap[0] = -1;
    
    s[0]->alphaExist[0] = 1;
    s[0]->counts[0] = 1;
    s[0]->n = 1;
    
    // STEP
    s[0]->step = 10;
    
    s[0]->rescale = rescale;
    
    return s;
    
}

stream_model* initialize_stream_model_pos_alpha(uint32_t rescale){
    
    uint32_t i = 0, j = 0;
    uint32_t context_size = 4;
    
    stream_model *s = (stream_model*) calloc(context_size, sizeof(stream_model));
    
    for (i = 0; i < context_size; i++) {
        
        s[i] = (stream_model) calloc(1, sizeof(struct stream_model_t));
        
        // Allocate memory
        s[i]->counts = (uint32_t*) calloc(257, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[i]->counts += 1;
        
        s[i]->alphabetCard = 256;
        
        s[i]->n = 0;
        for (j = 0; j < 256; j++) {
            s[i]->counts[j] = 1;
            s[i]->n += s[0]->counts[j];
        }
        
        // STEP
        s[i]->step = 10;
        
        //rescale bound
        s[i]->rescale = rescale;
    }
    
    
    
    
    
    return s;
    
}

stream_model* initialize_stream_model_match(uint32_t rescale){
    
    uint32_t i = 0;
    uint32_t context_size = 256;
    
    stream_model *s = (stream_model*) calloc(context_size, sizeof(stream_model));
    
    for (i = 0; i < context_size; i++) {
        
        s[i] = (stream_model) calloc(1, sizeof(struct stream_model_t));
        
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

stream_model* initialize_stream_model_snps(uint32_t readLength, uint32_t rescale){
    
    stream_model *s;
    
    uint32_t i = 0;
    
    s = (stream_model*) calloc(1, sizeof(stream_model));
    
    s[0] = (stream_model) calloc(1, sizeof(struct stream_model_t));
    
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

stream_model* initialize_stream_model_indels(uint32_t readLength, uint32_t rescale){
    
    stream_model *s;
    
    s = (stream_model*) calloc(1, sizeof(stream_model));
    
    uint32_t i = 0;
    
    s[0] = (stream_model) calloc(1, sizeof(struct stream_model_t));
    
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

stream_model* initialize_stream_model_var(uint32_t readLength, uint32_t rescale){
    
    stream_model* s;
    
    uint32_t i = 0, j = 0;
    
    uint32_t num_models = 0xffff;
    
    s = (stream_model*) calloc(num_models, sizeof(stream_model));
    
    for (j = 0; j < num_models; j++) {
        
        s[j] = (stream_model) calloc(1, sizeof(struct stream_model_t));
        
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

stream_model* initialize_stream_model_chars(uint32_t rescale){
    
    stream_model* s;
    
    s = (stream_model*) calloc(6, sizeof(stream_model));
    
    uint32_t i = 0, j;
    
    for (j = 0; j < 6; j++) {
        
        s[j] = (stream_model) calloc(1, sizeof(struct stream_model_t));
        
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

/**
 * Initialize stats structures used for adaptive arithmetic coding based on
 * the number of contexts required to handle the set of conditional quantizers
 * that we have (one context per quantizer)
 */
stream_model *initialize_stream_model_qv(struct cond_quantizer_list_t *q_list, uint32_t rescale) {
    
    stream_model *s;
    uint32_t i = 0, j = 0, k = 0, model_idx = 0;
    
    uint32_t num_models = 0xffff;
    
    s = (stream_model*) calloc(num_models, sizeof(stream_model));
    
    // Allocate jagged array, one set of stats per column
    for (i = 0; i < q_list->columns; ++i) {
        // And for each column, one set of stats per low/high quantizer per previous context
        for (j = 0; j < 2*q_list->input_alphabets[i]->size; ++j) {
            // Finally each individual stat structure needs to be filled in uniformly
            model_idx = get_qv_model_index(i, j);
            //model_idx = ((i & 0xff) << 8 | (j & 0xff));
            s[model_idx] = (stream_model) calloc(1, sizeof(struct stream_model_t));
            s[model_idx]->counts = (uint32_t *) calloc(q_list->q[i][j]->output_alphabet->size, sizeof(uint32_t));
            
            // Initialize the quantizer's stats uniformly
            for (k = 0; k < q_list->q[i][j]->output_alphabet->size; k++) {
                s[model_idx]->counts[k] = 1;
            }
            s[model_idx]->n = q_list->q[i][j]->output_alphabet->size;
            s[model_idx]->alphabetCard = q_list->q[i][j]->output_alphabet->size;
            
            // Step size is 8 counts per symbol seen to speed convergence
            s[model_idx]->step = 8;
            
            //rescale bound
            s[model_idx]->rescale = rescale;
            
        }
    }
    
    return s;
}

stream_model *free_stream_model_qv(struct cond_quantizer_list_t *q_list, stream_model *s) {
    
    uint32_t i = 0, j = 0, model_idx = 0;
    
    
    
    // Allocate jagged array, one set of stats per column
    for (i = 0; i < q_list->columns; ++i) {
        // And for each column, one set of stats per low/high quantizer per previous context
        for (j = 0; j < 2*q_list->input_alphabets[i]->size; ++j) {
            // Finally each individual stat structure needs to be filled in uniformly
            model_idx = get_qv_model_index(i, j);
            //model_idx = ((i & 0xff) << 8 | (j & 0xff));
            free_model(s[model_idx]);
        }
    }
    
    return s;
}


/**
 *
 */
read_models alloc_read_models_t(uint32_t read_length){
    
    uint32_t rescale = 1 << 20;
    
    read_models rtn = calloc(1, sizeof(struct read_models_t));
    
    rtn->read_length = read_length;
    sprintf(rtn->_readLength, "%d", rtn->read_length);
    
    rtn->flag = initialize_stream_model_flag(rescale);
    rtn->pos_alpha = initialize_stream_model_pos_alpha(rescale);
    rtn->pos = initialize_stream_model_pos(rescale);
    rtn->match = initialize_stream_model_match(rescale);
    rtn->snps = initialize_stream_model_snps(rtn->read_length, rescale);
    rtn->indels = initialize_stream_model_indels(rtn->read_length, rescale);
    rtn->var = initialize_stream_model_var(rtn->read_length, rescale);
    rtn->chars = initialize_stream_model_chars(rescale);
    
    
    return rtn;
}

/**
 *
 */
id_models alloc_id_models_t(){
    
    uint32_t rescale = 1 << 20;
    
    id_models rtn = calloc(1, sizeof(struct id_models_t));
    
    rtn->alpha_len = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
    rtn->alpha_value = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
    rtn->chars = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
    rtn->integer = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID * 4, 256);
    rtn->zero_run = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 256);
    rtn->token_type = initialize_stream_model_id(rescale, MAX_NUMBER_TOKENS_ID, 10);
    
    return rtn;
}


/**
 *
 */
stream_model* initialize_stream_model_codebook(uint32_t rescale){
    
    // It is a byte-based model with the previous byte as context.
    uint32_t i = 0, j = 0;
    uint32_t context_size = 256;
    
    stream_model *s = (stream_model*) calloc(context_size, sizeof(stream_model));
    
    for (i = 0; i < context_size; i++) {
        
        s[i] = (stream_model) calloc(1, sizeof(struct stream_model_t));
        
        // Allocate memory
        s[i]->counts = (uint32_t*) calloc(257, sizeof(uint32_t));
        
        // An extra for the cumcounts
        s[i]->counts += 1;
        
        s[i]->alphabetCard = 256;
        
        s[i]->n = 0;
        for (j = 0; j < 256; j++) {
            s[i]->counts[j] = 1;
            s[i]->n += s[i]->counts[j];
        }
        
        // STEP
        s[i]->step = 1;
        
        //rescale bound
        s[i]->rescale = rescale;
    }
    
    
    return s;
    
}



/**
 *
 */
void initialize_qv_model(Arithmetic_stream as, qv_block qvBlock, uint8_t decompression){
    
    uint32_t rescale = 1 << 20;
    
    // We need to generate the Codebooks of the QVs in order to initialize the stream_model
    if (decompression) {
        
        // Read the codebook from the input stream
        read_codebooks(as, qvBlock);
    }
    else{
        
        // Allocate the training set
        qvBlock->training_stats = alloc_conditional_pmf_list(qvBlock->alphabet, qvBlock->columns);
        
        // Calculate the statistics of the training set and generate the codebook
        calculate_statistics(qvBlock);
        generate_codebooks(qvBlock);
        
        // Write the codebook for this block in the output stream
        write_codebooks(as, qvBlock);

    }
    
    //print_codebook(qvBlock->qlist);
    
    // initialize the model of the qv using the codebooks
    qvBlock->model = initialize_stream_model_qv(qvBlock->qlist, rescale);
}


