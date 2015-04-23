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

int compress_rname(Arithmetic_stream as, rname_models models, char *rname){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint32_t ctr = 0;
    
    if(strcmp(rname, prev_name) == 0){
        
        compress_uint8t(as, models->same_ref[0], 0);
        return 0;
        
    }
    
    else{
        compress_uint8t(as, models->same_ref[0], 1);
        while (*rname) {
            compress_uint8t(as, models->rname[prevChar], *rname);
            prev_name[ctr++] = *rname;
            prevChar = *rname++;
        }
        compress_uint8t(as, models->rname[prevChar], 0);
        prev_name[ctr] = 0;
        return 1;
    }
    
}

int compress_mapq(Arithmetic_stream as, mapq_models models, uint8_t mapq){
    
    compress_uint8t(as, models->mapq[0], mapq);
    
    return 0;
}

int decompress_mapq(Arithmetic_stream as, mapq_models models, uint8_t *mapq){
    
    *mapq = decompress_uint8t(as, models->mapq[0]);
    
    return 0;
}

int compress_rnext(Arithmetic_stream as, rnext_models models, char *rnext){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint32_t ctr = 0;
    
    switch (rnext[0]) {
        case '=':
            compress_uint8t(as, models->same_ref[0], 0);
            return 0;
        case '*':
            compress_uint8t(as, models->same_ref[0], 1);
            return 1;
        default:
            compress_uint8t(as, models->same_ref[0], 2);
            break;
    }
    
    
    while (*rnext) {
        compress_uint8t(as, models->rnext[prevChar], *rnext);
        prev_name[ctr++] = *rnext;
        prevChar = *rnext++;
    }
    compress_uint8t(as, models->rnext[prevChar], 0);
    prev_name[ctr] = 0;
    
    return 1;
}

int decompress_rnext(Arithmetic_stream as, rnext_models models, char *rnext){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint8_t equal_flag = 0, ch;
    
    uint32_t ctr = 0;
    
    equal_flag = decompress_uint8t(as, models->same_ref[0]);
    
    switch (equal_flag) {
        case 0:
            *rnext = '=', rnext++;
            *rnext = 0;
            return 0;
        case 1:
            *rnext = '*', rnext++;
            *rnext = 0;
            return 1;
        default:
            break;
    }
    
    while ( (ch = decompress_uint8t(as, models->rnext[prevChar])) ) {
        
        *rnext = ch, rnext++;
        prev_name[ctr++] = ch;
        prevChar = ch;
    }
    *rnext = 0;
    
    return 2;
    
}

int compress_pnext(Arithmetic_stream as, pnext_models models, uint32_t pos, int32_t tlen, uint32_t pnext, uint8_t rname_rnextDiff, char* cigar){
    
    uint32_t p = 0, pn = 0, cig_op_len = 0, readLength = 0;
    
    if (rname_rnextDiff) {
        
        compress_uint8t(as, models->raw_pnext[0], (pnext >> 0) & 0xff);
        compress_uint8t(as, models->raw_pnext[1], (pnext >> 8) & 0xff);
        compress_uint8t(as, models->raw_pnext[2], (pnext >> 16) & 0xff);
        compress_uint8t(as, models->raw_pnext[3], (pnext >> 24) & 0xff);
        
        return 0;
    }
    
    else
    
    if (tlen == 0) {
        if (pnext == pos) {
            compress_uint8t(as, models->assumption[0], 0);
        }
        
        else{
            compress_uint8t(as, models->assumption[0], 1);
            compress_pnext_raw(as, models, pos, pnext);
        }
        
        return 0;
    }
    else if (tlen > 0){
        
        p = pos;
        pn = pnext;
        
    }
    
    else{
        p = pnext;
        pn = pos;
        tlen = -tlen;
    }
    
    while (*cigar) {
        cig_op_len = atoi(cigar);
        cigar += compute_num_digits(cig_op_len);
        if (*cigar == 'M' || *cigar == 'D' || *cigar == 'N') {
            readLength += cig_op_len;
        }
        ++cigar;
    }
    
    if (pn == tlen + p - readLength) {
        compress_uint8t(as, models->assumption[0], 0);
    }
    else{
        compress_uint8t(as, models->assumption[0], 1);
        compress_pnext_raw(as, models, p, pn);
    }
    
        
    return 1;
}


int compress_pnext_raw(Arithmetic_stream as, pnext_models models, uint32_t pos, uint32_t pnext){
    
    uint32_t delta = 0;
    
    
    
    if (pnext == 0) {
        compress_uint8t(as, models->zero[0], 0);
        return 0;
    }
    else{
        compress_uint8t(as, models->zero[0], 1);
        
        if (pnext > pos){
            
            delta = pnext - pos;
            compress_uint8t(as, models->sign[0], 0);
        }
        else {
            delta = pos - pnext;
            compress_uint8t(as, models->sign[0], 1);
        }
        
        compress_uint8t(as, models->diff_pnext[0], (delta >> 0) & 0xff);
        compress_uint8t(as, models->diff_pnext[1], (delta >> 8) & 0xff);
        compress_uint8t(as, models->diff_pnext[2], (delta >> 16) & 0xff);
        compress_uint8t(as, models->diff_pnext[3], (delta >> 24) & 0xff);
    }
    
    return 1;
}


int decompress_pnext(Arithmetic_stream as, pnext_models models, uint32_t pos, int32_t tlen, uint32_t readLength, uint32_t *pnext, uint8_t rname_rnextDiff, char* cigar){
    
    uint32_t delta = 0, comp_flag = 0;
    
    
    
   
    comp_flag = decompress_uint8t(as, models->zero[0]);
    
    if ( comp_flag == 0) {
        *pnext = 0;
        return 0;
    }
    else{
        comp_flag = decompress_uint8t(as, models->sign[0]);
        
        delta |= (decompress_uint8t(as, models->raw_pnext[0]) & 0xff) << 0;
        delta |= (decompress_uint8t(as, models->raw_pnext[1]) & 0xff) << 8;
        delta |= (decompress_uint8t(as, models->raw_pnext[2]) & 0xff) << 16;
        delta |= (decompress_uint8t(as, models->raw_pnext[3]) & 0xff) << 24;
        
        if (comp_flag == 0){
            
            *pnext = delta + pos;
        }
        else {
            
            *pnext = pos - delta;
        }
    }
    
    return 1;
}

int compress_tlen(Arithmetic_stream as, tlen_models models, int32_t tlen){
    
    int32_t delta = 0;
    
    if (tlen == 0) {
        compress_uint8t(as, models->zero[0], 0);
        return 0;
    }
    else{
        compress_uint8t(as, models->zero[0], 1);
        
        if (tlen > 0){
            delta = tlen;
            compress_uint8t(as, models->sign[0], 0);
        }
        else {
            delta = -tlen;
            compress_uint8t(as, models->sign[0], 1);
            
        }
        compress_uint8t(as, models->tlen[0], (delta >> 0) & 0xff);
        compress_uint8t(as, models->tlen[1], (delta >> 8) & 0xff);
        compress_uint8t(as, models->tlen[2], (delta >> 16) & 0xff);
        compress_uint8t(as, models->tlen[3], (delta >> 24) & 0xff);
    }
    
    return 1;
}
int decompress_tlen(Arithmetic_stream as, tlen_models models, int32_t* tlen){
    
    int32_t delta = 0;
    uint32_t decomp_flag = 0;
    
    decomp_flag = decompress_uint8t(as, models->zero[0]);
    
    if ( decomp_flag == 0) {
        *tlen = 0;
        return 0;
    }
    else{
        
        decomp_flag = decompress_uint8t(as, models->sign[0]);
        
        delta |= (decompress_uint8t(as, models->tlen[0]) & 0xff) << 0;
        delta |= (decompress_uint8t(as, models->tlen[1]) & 0xff) << 8;
        delta |= (decompress_uint8t(as, models->tlen[2]) & 0xff) << 16;
        delta |= (decompress_uint8t(as, models->tlen[3]) & 0xff) << 24;
        
        if (decomp_flag == 0)
            *tlen = delta;
        else
            *tlen = -delta;
    }
    
    return 1;
}



int compress_id(Arithmetic_stream as, id_models models, char *id){
    
    static char prev_ID[1024] = {0};
    static uint32_t prev_tokens_ptr[1024] = {0};
    uint8_t token_len = 0, match_len = 0;
    uint32_t i = 0, k = 0, token_ctr = 0, digit_value = 0, digit_model = 0, prev_digit = 0;
    int delta = 0;
    
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
        // Check if the token is a number smaller than (1<<30)
        else if (isdigit(*id_ptr)) {
            
            digit_value = (*id_ptr - '0');
            prev_digit = prev_ID[prev_tokens_ptr[token_ctr] + token_len -1] - '0';
            
            while ( isdigit(*id_ptr_tok) && digit_value < (1<<30) ){
                digit_value = digit_value * 10 + (*id_ptr_tok - '0');
                prev_digit = prev_digit * 10 + (prev_ID[prev_tokens_ptr[token_ctr] + token_len] - '0');
                // compare with the same token from previous ID
                match_len += (*id_ptr_tok == prev_ID[prev_tokens_ptr[token_ctr] + token_len]), token_len++, id_ptr_tok++;
                
            }
            if ( match_len == token_len && !isdigit(prev_ID[prev_tokens_ptr[token_ctr] + token_len]) ) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else if ( (delta = (digit_value - prev_digit)) < 256 && delta > 0){
                compress_uint8t(as, models->token_type[token_ctr], ID_DELTA);
                compress_uint8t(as, models->delta[token_ctr], delta);
                
            }
            else {
                // Encode a token type ID_DIGIT and the value (byte-based)
                compress_uint8t(as, models->token_type[token_ctr], ID_DIGIT);
                digit_model = (token_ctr << 2);
                compress_uint8t(as, models->integer[digit_model | 0], (digit_value >> 0) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 1], (digit_value >> 8) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 2], (digit_value >> 16) & 0xff);
                compress_uint8t(as, models->integer[digit_model | 3], (digit_value >> 24) & 0xff);
            }
        }
        else {
            
            // compare with the same token from previous ID
            //match_len += (*id_ptr == prev_ID[prev_tokens_ptr[token_ctr]]);
            
            if (match_len == token_len) {
                // The token is the same as last ID
                // Encode a token_type ID_MATCH
                compress_uint8t(as, models->token_type[token_ctr], ID_MATCH);
                
            }
            else {
                // Encode a token type ID_CHAR and the char
                compress_uint8t(as, models->token_type[token_ctr], ID_CHAR);
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


int decompress_rname(Arithmetic_stream as, rname_models models, char *rname){
    
    static char prev_name[1024] = {0};
    static int prevChar = 0;
    
    uint8_t chr_change = 0, ch;
    
    uint32_t ctr = 0;
    
    chr_change = decompress_uint8t(as, models->same_ref[0]);
    
    if (chr_change) {
        
        while ( (ch = decompress_uint8t(as, models->rname[prevChar])) ) {
            
            if (ch == '\n') {
                return -1;
            }
            prev_name[ctr++] = ch;
            prevChar = ch;
            *rname = ch, rname++;
        }

    }
    
    return chr_change;
    
}




int decompress_id(Arithmetic_stream as, id_models model, char *id){
    
    static char prev_ID[1024] = {0};
    static uint32_t prev_tokens_ptr[1024] = {0};
    static uint32_t prev_tokens_len[1024] = {0};
    //char id[1024] = {0};
    uint8_t token_len = 0;
    uint32_t i = 0, k = 0, token_ctr = 0, digit_value = 0;
    uint32_t delta = 0;
    
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
                    id[i+k] = decompress_uint8t(as, model->alpha_value[token_ctr]);
                }
                break;
            case ID_DIGIT:
                digit_value = 0;
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 0]) & 0xff ) << 0);
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 1]) & 0xff ) << 8);
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 2]) & 0xff ) << 16);
                digit_value |= (( decompress_uint8t(as, model->integer[(token_ctr << 2) | 3]) & 0xff ) << 24);
                sprintf(id+i, "%u", digit_value);
                token_len = compute_num_digits(digit_value);
                break;
            case ID_DELTA:
                digit_value = 0;
                delta = decompress_uint8t(as, model->delta[token_ctr]);
                memcpy(id+i, &(prev_ID[prev_tokens_ptr[token_ctr]]), prev_tokens_len[token_ctr]);
                digit_value = atoi(id+i) + delta;
                sprintf(id+i, "%u", digit_value);
                token_len = compute_num_digits(digit_value);
                break;
            case ID_ZEROS:
                token_len = decompress_uint8t(as, model->zero_run[token_ctr]);
                memset(id+i, '0', token_len), i += token_len;
                break;
            case ID_CHAR:
                id[i] = decompress_uint8t(as, model->chars[token_ctr]);
                token_len = 1;
                break;
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
    //id[i++] = '\n';
    //putc('@', fs);
    //fwrite(id, i, sizeof(char), fs);
    
    return 1;
}

