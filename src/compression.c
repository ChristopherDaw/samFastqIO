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
    uint8_t chr_change;
    
    float counts_rname = 0, counts_id = 0, counts_reads = 0, counts_qv = 0, counts_load = 0, counts_qv_model = 0, counts_quant = 0;
    
    clock_t begin;
    
    symbol_t *qArray;
    
    // Load the data from the first block
    //printf("Loading block of data into memory...\n");
    begin = clock();
    load_sam_block(samBlock);
    counts_load = (float)(clock() - begin)/CLOCKS_PER_SEC;
    // Compute the codebook and initialize the QV model
    //printf("Computing the codebook for the QVs...\n");
    begin = clock();
    initialize_qv_model(as, samBlock->QVs, COMPRESSION);
    //initialize_stream_model_qv_full(samBlock->QVs->model, samBlock->QVs->qlist);
    counts_qv_model = (float)(clock() - begin)/CLOCKS_PER_SEC;
    
    //begin = clock();
    //quantize_block(samBlock->QVs, samBlock->read_length);
    //printf("lala %f\n",  (float)(clock() - begin)/CLOCKS_PER_SEC);
    //qArray = copy_qlis_to_array(samBlock->QVs);
    
    //printf("Compressing the block...\n");
    // Loop over the lines of the sam block
    for (i = 0; i < samBlock->block_length; i++) {
        
        // Compress sam line
        begin = clock();
        chr_change = compress_rname(as, samBlock->rnames->models, samBlock->rnames->rnames[i]);
        counts_rname += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        if (chr_change == 1){
            
            //printf("Chromosome %d decompressed.\n", ++chrCtr);
            
            // Store Ref sequence in memory
            store_reference_in_memory(samBlock->fref);
            
            // Clean snpInRef vector and reset cumsumP
            cumsumP = 0;
            memset(snpInRef, 0, MAX_BP_CHR);
        }
        begin = clock();
        compress_id(as, samBlock->IDs->models, samBlock->IDs->IDs[i]);
        counts_id += (float)(clock() - begin)/CLOCKS_PER_SEC;
        begin = clock();
        compress_read(as, samBlock->reads->models, &(samBlock->reads->lines[i]), chr_change);
        counts_reads += (float)(clock() - begin)/CLOCKS_PER_SEC;
        begin = clock();
        quantize_line(samBlock->QVs, &(samBlock->QVs->qv_lines[i]), samBlock->read_length);
        counts_quant += (float)(clock() - begin)/CLOCKS_PER_SEC;
        begin = clock();
        //QVs_compress2(as, samBlock->QVs->qlist->input_alphabets, samBlock->QVs->qlist->qratio, &(samBlock->QVs->qv_lines[i]), qArray, samBlock->QVs->model, &(samBlock->QVs->well));
        QVs_compress3(as, samBlock->QVs->model, &(samBlock->QVs->qv_lines[i]));
        counts_qv += (float)(clock() - begin)/CLOCKS_PER_SEC;
    }
    
    printf("Time load \t %f\n", ((float)counts_load));
    printf("Time qv_model \t %f\n", ((float)counts_qv_model));
    printf("Time rname \t %f\n", ((float)counts_rname));
    printf("Time id \t %f\n", ((float)counts_id));
    printf("Time reads \t %f\n", ((float)counts_reads));
    
    printf("Time quant \t %f\n", ((float)counts_quant));
    printf("Time qv \t %f\n", ((float)counts_qv));
    
    
    
    // Check if we are in the last block
    if (i < MAX_LINES_PER_BLOCK){
        compress_rname(as, samBlock->rnames->models, "\n");
        return 0;
    }
    
    return 1;
}

int decompress_block(Arithmetic_stream as, sam_block samBlock){
    
    unsigned int i = 0;
    
    int32_t chr_change = 0;
    
    uint32_t decompression_flag = 0;
    
    // initialize the QV model
    printf("Computing the codebook for the QVs...\n");
    initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
    
    printf("Decompressing the block...\n");
    // Loop over the lines of the sam block
    for (i = 0; i < samBlock->block_length; i++) {
        
        chr_change = decompress_rname(as, samBlock->rnames->models, NULL);
        
        if (chr_change == -1)
            return 0;
        
        if (chr_change == 1){
            
            //printf("Chromosome %d decompressed.\n", ++chrCtr);
            
            // Store Ref sequence in memory
            store_reference_in_memory(samBlock->fref);
            
            // Clean snpInRef vector and reset cumsumP
            cumsumP = 0;
            memset(snpInRef, 0, MAX_BP_CHR);
        }
        
        decompress_id(as, samBlock->IDs->models, samBlock->fs);
        
        decompression_flag = decompress_read(as,samBlock, chr_change);
        
        QVs_decompress(as, samBlock->QVs, samBlock->fs, decompression_flag);
        
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
    
    Arithmetic_stream as = alloc_arithmetic_stream(info.mode);
    
    sam_block samBlock = alloc_sam_block_t(as, info.fsam, info.fref, &opts, info.mode);
    
    // Compress the blocks
    while(compress_block(as, samBlock)){
        reset_QV_block(samBlock->QVs, info.mode);
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
    clock_t begin = clock();
    clock_t ticks;
    
    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
    
    Arithmetic_stream as = alloc_arithmetic_stream(DECOMPRESSION);
    
    sam_block samBlock = alloc_sam_block_t(as, info->fsam, info->fref, info->qv_opts, DECOMPRESSION);
    
    // Start the decompression
    
    // Decompress the blocks
    while(decompress_block(as, samBlock)){
        reset_QV_block(samBlock->QVs, DECOMPRESSION);
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
