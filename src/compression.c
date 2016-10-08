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

int print_line(struct sam_line_t *sline, uint8_t print_mode, FILE *fs){
    
    char foo[] = "CIGAR";
    
    int32_t i = 0;
    
    switch (print_mode) {
        case 0:
            if ((sline->flag & 16) == 16) {
                // We need to inverse the read and QV
                fprintf(fs, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t",
                        sline->ID,
                        sline->flag,
                        sline->rname,
                        sline->pos,
                        sline->mapq,
                        //sline->cigar,
                        foo,
                        sline->rnext,
                        sline->pnext,
                        sline->tlen
                        );
                for (i = sline->readLength - 1; i >=0 ; --i)
                    fputc(bp_complement(sline->read[i]), fs);
                    
                fputc('\t', fs);
                
                for (i = sline->readLength - 1; i >=0 ; --i)
                    fputc(sline->quals[i], fs);
                
                fputc('\n', fs);
            }
            else{
                fprintf(fs, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
                   sline->ID,
                   sline->flag,
                   sline->rname,
                   sline->pos,
                   sline->mapq,
                   //sline->cigar,
                   foo,
                   sline->rnext,
                   sline->pnext,
                   sline->tlen,
                   sline->read,
                   sline->quals
                   );
            }
            break;
            
        default:
            break;
    }
    return 0;
}

int compress_line(Arithmetic_stream as, sam_block samBlock, uint8_t lossiness)  {
    
    uint8_t chr_change;
    
    // Load the data from the file
    if(load_sam_line(samBlock))
        return 0;
    
    if ( (samBlock->reads->lines->invFlag & 4) == 4) {
        return 1;
    }
        
    // Compress sam line
    
    chr_change = compress_rname(as, samBlock->rnames->models, *samBlock->rnames->rnames);
        
    if (chr_change == 1){
        // Store Ref sequence in memory
        store_reference_in_memory(samBlock->fref);
        // Clean snpInRef vector and reset cumsumP
        cumsumP = 0;
        memset(snpInRef, 0, MAX_BP_CHR);
    }
    
    compress_id(as, samBlock->IDs->models, *samBlock->IDs->IDs);
    
    compress_mapq(as, samBlock->mapq->models, *samBlock->mapq->mapq);
    
    compress_rnext(as, samBlock->rnext->models, *samBlock->rnext->rnext);
    
    compress_read(as, samBlock->reads->models, samBlock->reads->lines, chr_change);
    
    compress_tlen(as, samBlock->tlen->models, *samBlock->tlen->tlen);
    
    compress_pnext_raw(as, samBlock->pnext->models,  samBlock->reads->lines->pos, *samBlock->pnext->pnext);
    
    //compress_pnext(as, samBlock->pnext->models, samBlock->reads->lines->pos, *samBlock->tlen->tlen, *samBlock->pnext->pnext, (*samBlock->rnext->rnext[0] != '='), samBlock->reads->lines->cigar);
    
    if (lossiness == LOSSY)
        QVs_compress(as, samBlock->QVs, samBlock->QVs->qArray);
    else
        QVs_compress_lossless(as, samBlock->QVs->model, samBlock->QVs->qv_lines);
    
    
    return 1;
}

int decompress_line(Arithmetic_stream as, sam_block samBlock, uint8_t lossiness){
    
    int32_t chr_change = 0;
    
    uint32_t decompression_flag = 0;
    
    struct sam_line_t sline;
    
    sline.readLength = samBlock->read_length;
    
    //printf("Decompressing the block...\n");
    // Loop over the lines of the sam block
        
    chr_change = decompress_rname(as, samBlock->rnames->models, sline.rname);
        
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
        
    decompress_id(as, samBlock->IDs->models, sline.ID);
    
    decompress_mapq(as, samBlock->mapq->models, &sline.mapq);
    
    decompress_rnext(as, samBlock->rnext->models, sline.rnext);
        
    decompression_flag = decompress_read(as,samBlock, chr_change, &sline);
    
    decompress_tlen(as, samBlock->tlen->models, &sline.tlen);
    
    decompress_pnext(as, samBlock->pnext->models, sline.pos, sline.tlen, samBlock->read_length, &sline.pnext, sline.rnext[0] != '=', NULL);
        
    if (lossiness == LOSSY) {
            QVs_decompress(as, samBlock->QVs, decompression_flag, sline.quals);
    }
    else
        QVs_decompress_lossless(as, samBlock->QVs, decompression_flag, sline.quals);
    
    print_line(&sline, 0, samBlock->fs);
    
    return 1;
}




/*int compress_block(Arithmetic_stream as, sam_block samBlock){
    
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
    
    begin = clock();
    //quantize_block(samBlock->QVs, samBlock->read_length);
    qArray = copy_qlis_to_array(samBlock->QVs);
    printf("lala %f\n",  (float)(clock() - begin)/CLOCKS_PER_SEC);
    
    //printf("Compressing the block...\n");
    // Loop over the lines of the sam block
    for (i = 0; i < samBlock->block_length; i++) {
        
        // Compress sam line
        begin = clock();
        chr_change = compress_rname(as, samBlock->rnames->models, samBlock->rnames->rnames[i]);
        counts_rname += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        if (chr_change == 1){
            
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
        //begin = clock();
        //quantize_line(samBlock->QVs, &(samBlock->QVs->qv_lines[i]), samBlock->read_length);
        //counts_quant += (float)(clock() - begin)/CLOCKS_PER_SEC;
        begin = clock();
        //QVs_compress(as, samBlock->QVs, &(samBlock->QVs->qv_lines[i]), qArray);
//        QVs_compress(as, samBlock->QVs, i, qArray);
        //QVs_compress2(as, samBlock->QVs->qlist->input_alphabets, samBlock->QVs->qlist->qratio, &(samBlock->QVs->qv_lines[i]), qArray, samBlock->QVs->model, &(samBlock->QVs->well));
        //QVs_compress3(as, samBlock->QVs->model, &(samBlock->QVs->qv_lines[i]));
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
    
    //uint32_t decompression_flag = 0;
    
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
        
        //decompress_id(as, samBlock->IDs->models, samBlock->fs);
        
        //decompression_flag = decompress_read(as,samBlock, chr_change);
        
        //QVs_decompress(as, samBlock->QVs, samBlock->fs, decompression_flag);
        
    }
    
    return 1;
}*/



void* compress(void *thread_info){
    
    uint64_t compress_file_size = 0;
    clock_t begin;
    clock_t ticks;
    
    uint32_t lineCtr = 0;
    
    printf("Compressing...\n");
    begin = clock();
    
    struct compressor_info_t info = *((struct compressor_info_t *)thread_info);
    
    struct qv_options_t opts = *(info.qv_opts);
    
    // Allocs the Arithmetic and the I/O stream
    Arithmetic_stream as = alloc_arithmetic_stream(info.mode, info.fcomp);
    
    // Allocs the different blocks and all the models for the Arithmetic
    sam_block samBlock = alloc_sam_models(as, info.fsam, info.fref, &opts, info.mode);
    
    if (info.lossiness == LOSSY) {
        compress_int(as, samBlock->codebook_model, LOSSY);
        initialize_qv_model(as, samBlock->QVs, COMPRESSION);
    }
    else
        compress_int(as, samBlock->codebook_model, LOSSLESS);
    
    while (compress_line(as, samBlock, info.lossiness)) {
        ++lineCtr;
    }
    // Load and compress the blocks
    //while(compress_block(as, samBlock)){
    //    reset_QV_block(samBlock->QVs, info.mode);
    
    //   n += samBlock->block_length;
    //}
    
    // Check if we are in the last block
    compress_rname(as, samBlock->rnames->models, "\n");
    
    //end the compression
    compress_file_size = encoder_last_step(as);
    
    printf("Final Size: %lld\n", compress_file_size);
    
    //printf("%f Million reads compressed using %f MB.\n", (double)n/1000000.0, (double)compress_file_size/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    fclose(info.fsam);
    
    ticks = clock() - begin;
    
    printf("Compression took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    
    //pthread_exit(NULL);
    return NULL;
}


void* decompress(void *thread_info){
    
    uint64_t n = 0;
    clock_t begin = clock();
    clock_t ticks;
    
    struct compressor_info_t *info = (struct compressor_info_t *)thread_info;
    
    Arithmetic_stream as = alloc_arithmetic_stream(info->mode, info->fcomp);
    
    sam_block samBlock = alloc_sam_models(as, info->fsam, info->fref, info->qv_opts, DECOMPRESSION);
    
    info->lossiness = decompress_int(as, samBlock->codebook_model);
    
    // Start the decompression
    // initialize the QV model
    if (info->lossiness == LOSSY) {
        initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
    }
    
    // Decompress the blocks
    while(decompress_line(as, samBlock, info->lossiness)){
        //reset_QV_block(samBlock->QVs, DECOMPRESSION);
        n++;
    }
    
    n += samBlock->block_length;
    
    //end the decompression
    //compress_file_size = encoder_last_step(as);
    
    ticks = clock() - begin;
    
    printf("Decompression took %f\n", ((float)ticks)/CLOCKS_PER_SEC);
    
    //printf("%f Million reads decompressed.\n", (double)n/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    fclose(info->fsam);
    fclose(info->fref);
    
    return NULL;
}
