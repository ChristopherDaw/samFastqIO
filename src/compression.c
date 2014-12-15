//
//  compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include "sam_block.h"
#include "read_compression.h"

int compress_block(Arithmetic_stream as, sam_block samBlock){
    
    unsigned int i = 0;
    
    // Load the data from the first block
    //printf("Loading block of data into memory...\n");
    load_sam_block(samBlock);
    
    // Compute the codebook and initialize the QV model
    //printf("Computing the codebook for the QVs...\n");
//    initialize_qv_model(as, samBlock->QVs, COMPRESSION);
    
    
    //printf("Compressing the block...\n");
    // Loop over the lines of the sam block
    for (i = 0; i < samBlock->block_length; i++) {
        
        // Compress sam line
//        compress_read(as, samBlock->reads->models, &(samBlock->reads->lines[i]));
//        QVs_compress(as, samBlock->QVs, &(samBlock->QVs->qv_lines[i]));
        compress_id(as, samBlock->IDs->models, samBlock->IDs->IDs[i]);
        
    }
    
    // Check if we are in the last block
    if (i < MAX_LINES_PER_BLOCK){
        compress_pos(as, samBlock->reads->models->pos, samBlock->reads->models->pos_alpha, END_GENOME_FLAG);
        return 0;
    }
    
    return 1;
}

int decompress_block(Arithmetic_stream as, sam_block samBlock){
    
    static uint32_t chrCtr = 0;
    unsigned int i = 0;
    
    uint32_t decompression_flag = 0;
    
    // initialize the QV model
    printf("Computing the codebook for the QVs...\n");
//    initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
    
    printf("Decompressing the block...\n");
    // Loop over the lines of the sam block
    for (i = 0; i < samBlock->block_length; i++) {
        
//        decompression_flag = decompress_read(as,samBlock);
        
        if (decompression_flag == CHR_CHANGE_FLAG){
            
            printf("Chromosome %d decompressed.\n", ++chrCtr);
            
            // Store Ref sequence in memory
            store_reference_in_memory(samBlock->fref);
            
            // TODO clean this up
            // Clean snpInRef vector and reset cumsumP
            cumsumP = 0;
            for (i=0; i<MAX_BP_CHR; i++){
                snpInRef[i] = 0;
            }
            
            continue;
        }
        
        if (decompression_flag == END_GENOME_FLAG)
            return 0;
        
//        QVs_decompress(as, samBlock->QVs, samBlock->fs, decompression_flag);
        
        decompress_id(as, samBlock->IDs->models, samBlock->fs);
        
    }
    
    return 1;
}



void* compress(void *thread_info){
    
    uint64_t compress_file_size = 0, n = 0;
    time_t begin;
    time_t end;
    
    //printf("Compressing...\n");
    time(&begin);
    
    struct compressor_info_t info = *((struct compressor_info_t *)thread_info);
    
    struct qv_options_t opts = *(info.qv_opts);
    
    Arithmetic_stream as = alloc_arithmetic_stream(COMPRESSION);
    
    sam_block samBlock = alloc_sam_block_t(as, info.fsam, NULL, &opts, COMPRESSION);
    
    // Compress the blocks
    while(compress_block(as, samBlock)){
//        reset_QV_block(samBlock->QVs, COMPRESSION);
        n += samBlock->block_length;
    }
    
    n += samBlock->block_length;
    
    //end the compression
    compress_file_size = encoder_last_step(as);
    
    //printf("%f Million reads compressed using %f MB.\n", (double)n/1000000.0, (double)compress_file_size/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    fclose(info.fsam);
    
    time(&end);
    
    printf("Compression time elapsed: %ld seconds\n", end-begin);
    
    pthread_exit(NULL);
}


void* decompress(void *thread_info){
    
    uint64_t n = 0;
    uint32_t i = 0;
    uint8_t reference_flag = 0;
    clock_t begin = clock();
    clock_t ticks;
    
    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
    
    Arithmetic_stream as = alloc_arithmetic_stream(DECOMPRESSION);
    
    sam_block samBlock = alloc_sam_block_t(as, info->fsam, info->fref, info->qv_opts, DECOMPRESSION);
    
    // Start the decompression
    
    // Store Ref sequence in memory
    reference_flag = store_reference_in_memory(samBlock->fref);
    
    // TODO clean this up
    // Clean snpInRef vector and reset cumsumP
    cumsumP = 0;
    for (i=0; i<MAX_BP_CHR; i++){
        snpInRef[i] = 0;
    }
    
    // Decompress the blocks
    while(decompress_block(as, samBlock)){
//        reset_QV_block(samBlock->QVs, DECOMPRESSION);
        n += samBlock->block_length;
    }
    
    n += samBlock->block_length;
    
    //end the decompression
    //compress_file_size = encoder_last_step(as);
    
    ticks = clock() - begin;
    
    //printf("Decompression took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    
    //printf("%f Million reads decompressed.\n", (double)n/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    fclose(info->fsam);
    fclose(info->fref);
    
    return NULL;
}
