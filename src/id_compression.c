//
//  id_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/10/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

// Compression of the IDs -- work based on the compression of ID in Samcomp by Mahonney and Bonfiled (2012)

#include <stdio.h>
#include "sam_block.h"

uint8_t decompress_uint8t(Arithmetic_stream as, stream_model model){
    
    // Send the value to the Arithmetic Stream
    uint8_t c = read_value_from_as(as, model);
    
    // Update model
    update_model(model, c);
    
    return c;
    
}


int compress_uint8t(Arithmetic_stream as, stream_model model, uint8_t c){
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(as, model, c);
    
    // Update model
    update_model(model, c);
    
    return 1;
    
}

int compress_id(Arithmetic_stream as, id_models models, char *id){
    
    static char prev_ID[1024] = {0};
    static uint32_t prev_tokens_ptr[1024] = {0};
    uint8_t token_len = 0, match_len = 0;
    uint32_t i = 0, k = 0, token_ctr = 0, digit_value = 0, digit_model = 0;
    
    char *id_ptr = id, *id_ptr_tok = NULL;
    
    while (*id_ptr != 0) {
        match_len += (*id_ptr == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++;
        id_ptr_tok = id_ptr + 1;
        
        // Check if the token is a alphabetic word
        if (isalpha(*id_ptr)) {
            
            while ( isalpha( *id_ptr_tok) ){
                
                // compare with the same token from previous ID
                match_len += (*id_ptr_tok == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++, id_ptr_tok++;
                
            }
            if (match_len == token_len && !isalpha(prev_ID[prev_tokens_ptr[token_ctr] + token_len])) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_ALPHA, the length of the string and the string
                compress_uint8t(as, models->token_type[token_ctr], ID_ALPHA);
                compress_uint8t(as, models->alpha_len[token_ctr], token_len);
                for (k = 0; k < token_len; k++) {
                    compress_uint8t(as, models->alpha_value[token_ctr], *(id_ptr+k));
                }
            }
            
        }
        // check if the token is a run of zeros
        else if (*id_ptr == '0') {
            
            while ( *id_ptr_tok == '0' ){
                
                // compare with the same token from previous ID
                match_len += ('0' == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++, id_ptr_tok++;
                
            }
            if (match_len == token_len && prev_ID[prev_tokens_ptr[token_ctr] + token_len] != '0') {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_ZEROS and the length of the zeros
                compress_uint8t(as, models->token_type[token_ctr], ID_ZEROS);
                compress_uint8t(as, models->zero_run[token_ctr], token_len);
            }
            
        }
        // Check if the token is a number smaller than
        else if (isdigit(*id_ptr)) {
            
            digit_value = digit_value * 10 + (*id_ptr - '0');
            while ( isdigit(*id_ptr_tok) ){
                digit_value = digit_value * 10 + (*id_ptr_tok - '0');
                // compare with the same token from previous ID
                match_len += (*id_ptr_tok == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++, id_ptr_tok++;
                
            }
            if ( match_len == token_len && !isdigit(prev_ID[prev_tokens_ptr[token_ctr] + token_len]) ) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_DIGIT, the value
                digit_model = (token_ctr << 2);
                compress_uint8t(as, models->integer[digit_model | 0], (digit_value >> 0) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 1], (digit_value >> 8) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 2], (digit_value >> 16) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 3], (digit_value >> 24) & 0xff);
            }
        }
        else {
            
            // compare with the same token from previous ID
            match_len += (*id_ptr == prev_ID[prev_tokens_ptr[token_ctr]]);
            
            if (match_len == token_len) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_CHAR
                compress_uint8t(as, models->chars[token_ctr], *id_ptr);
            }

        }
        
        prev_tokens_ptr[token_ctr] = i;
        i += token_len;
        id_ptr = id_ptr_tok;
        match_len = 0;
        token_len = 0;
        token_ctr++;
        
    }
    strcpy(prev_ID, id);
    compress_uint8t(as, models->token_type[token_ctr], ID_END);
    
    return 1;
}




int decompress_id(Arithmetic_stream as, id_models model, FILE *fs){
    
    static char prev_ID[1024] = {0};
    static uint32_t prev_tokens_ptr[1024] = {0};
    static uint32_t prev_tokens_len[1024] = {0};
    char id[1024] = {0};
    uint8_t token_len = 0, match_len = 0;
    uint32_t i = 0, k = 0, token_ctr = 0, digit_value = 0, digit_model = 0;
    
    enum token_type tok;
    
    while ( (tok = decompress_uint8t(as, model->token_type[token_ctr])) != ID_END ) {
        
        switch (tok) {
            case ID_MATCH:
                memcpy(id+i, &(prev_ID[prev_tokens_ptr[token_ctr]]), prev_tokens_len[token_ctr]);
                token_len = prev_tokens_len[token_ctr];
                break;
            case ID_ALPHA:
                token_len = decompress_uint8t(as, model->alpha_len[token_ctr]);
                for (k = 0; k < token_len; k++) {
                    id[i] = decompress_uint8t(as, model->alpha_value[token_ctr]);
                }
                break;
            case ID_DIGIT:
                digit_value = 0;
                digit_value |= ( decompress_uint8t(as, model->integer[(token_ctr << 2) | 0]) & 0xff ) << 0;
                digit_value |= ( decompress_uint8t(as, model->integer[(token_ctr << 2) | 1]) & 0xff ) << 4;
                digit_value |= ( decompress_uint8t(as, model->integer[(token_ctr << 2) | 2]) & 0xff ) << 8;
                digit_value |= ( decompress_uint8t(as, model->integer[(token_ctr << 2) | 3]) & 0xff ) << 24;
                sprintf(id+i, "%ud", digit_value);
                break;
            case ID_ZEROS:
                token_len = decompress_uint8t(as, model->zero_run[token_len]);
                memset(id+i, '0', token_len), i += token_len;
                break;
            case ID_CHAR:
                id[i] = decompress_uint8t(as, model->chars[token_len]);
                token_len = 1;
                
            default:
                break;
        }
        
        prev_tokens_ptr[token_ctr] = i;
        prev_tokens_len[token_ctr] = token_len;
        i += token_len;
        token_len = 0;
        token_ctr++;
    }
    
    strcpy(prev_ID, id);
    
    return 1;
}

