//
//  sam_file_allocation.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/18/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "sam_line.h"


uint32_t get_read_length(FILE *f){
    
    int ch, header_bytes = 0;
    char buffer[1024]; // 1 KB buffer
    
    // We use this opportunity to remove the headers
    while ((ch = getc(f)) == '@') {
        fgets(buffer, 1024, f);
        header_bytes += strlen(buffer);
    }
    
    // rewind the file pointer to be at the beginning of the first read
    fseek(f, header_bytes, SEEK_SET);
    
    // Extract the first read
    fscanf(f, "%*s %*d %*s %*d %*d %s", buffer);
    
    // rewind the file pointer to be at the beginning of the first read
    fseek(f, header_bytes, SEEK_SET);
    
    
    return (uint32_t)strlen(buffer);
    
}


/**
 *
 */
read_block alloc_read_block_t(uint32_t read_length){
    
    uint32_t i = 0;
    
    read_block rf = (read_block) calloc(1, sizeof(struct read_block_t));
    
    rf->block_length = MAX_READS_PER_BLOCK;
    
    rf->lines = (read_line) calloc(rf->block_length, sizeof(struct read_line_t));
    

    
    // allocate the memory for each of the lines
    for (i = 0; i < rf->block_length; i++) {
        
        rf->lines[i].cigar = (char*) calloc(1, 2*read_length);
        rf->lines[i].edits = (char*) calloc(1, 2*read_length);
        rf->lines[i].read = (char*) calloc(1, read_length + 1);
    }
    
    return rf;
}

/**
 *
 */
qv_block alloc_qv_block_t(struct qv_options_t *opts, uint32_t read_length){
    
    qv_block qv_info;
    
    uint32_t i = 0;
    
    qv_info = (qv_block) calloc(1, sizeof(struct qv_block_t));
    
    // Allocate the QV alphabet and the distortion matrix
    struct distortion_t *dist = generate_distortion_matrix(41, opts->distortion);
    struct alphabet_t *alphabet = alloc_alphabet(41);
    
    // allocate the memory for each of the lines
    for (i = 0; i < qv_info->block_length; i++) {
        
        qv_info->lines[i].data = (symbol_t*) calloc(qv_info->columns, sizeof(symbol_t));
        qv_info->lines[i].columns = read_length;
    }
    
    // set the alphabet and distortion
    qv_info->alphabet = alphabet;
    qv_info->dist = dist;
    
    qv_info->opts = opts;
    
    return qv_info;
    
}

sam_file alloc_sam_file_t(FILE * fs, struct qv_options_t *qv_opts){
    
    sam_file sf = (sam_file) calloc(1, sizeof(struct sam_file_t));
    
    sf->fs = fs;
    
    // get the read length and move file pointer after headers
    sf->read_length = get_read_length(fs);
    
    // Allocate the memory for the three parts: READS, QVs, ID
    sf->reads = alloc_read_block_t(sf->read_length);
    sf->QVs = alloc_qv_block_t(qv_opts, sf->read_length);
    
    // TODO: IDs
    return sf;
    
}
