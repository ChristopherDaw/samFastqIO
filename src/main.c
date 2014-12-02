//
//  main.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/4/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include <stdio.h>
#include <time.h>
#include "sam_block.h"
#include "Arithmetic_stream.h"
#include "read_compression.h"

int compress_block(Arithmetic_stream as, sam_block samBlock){
    
    unsigned int i = 0;
    
    // Load the data from the first block
    printf("Loading block of data into memory...\n");
    load_sam_block(samBlock);
    
    // Compute the codebook and initialize the QV model
    printf("Computing the codebook for the QVs...\n");
    initialize_qv_model(as, samBlock->QVs, COMPRESSION);
    
    
    printf("Compressing the block...\n");
    // Loop over the lines of the sam block
    for (i = 0; i < samBlock->block_length; i++) {
        
        // Compress sam line
        compress_read(as, samBlock->reads->models, &(samBlock->reads->lines[i]));
        QVs_compress(as, samBlock->QVs, &(samBlock->QVs->qv_lines[i]));
        
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
    initialize_qv_model(as, samBlock->QVs, DECOMPRESSION);
    
    printf("Decompressing the block...\n");
    // Loop over the lines of the sam block
    for (i = 0; i < samBlock->block_length; i++) {
        
        decompression_flag = decompress_read(as,samBlock);
        
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
        
        QVs_decompress(as, samBlock->QVs, samBlock->fs, decompression_flag);
        
        
    }
    
    return 1;
}



int compress(FILE *fin, struct remote_file_info info, struct qv_options_t *qv_opts){
    
    uint64_t compress_file_size = 0, n = 0;
    
    Arithmetic_stream as = alloc_arithmetic_stream(info,ARITHMETIC_WORD_LENGTH, COMPRESSION);
    
    sam_block samBlock = alloc_sam_block_t(as, fin, NULL, qv_opts, COMPRESSION);
    
    // Compress the blocks
    while(compress_block(as, samBlock)){
        reset_QV_block(samBlock->QVs, COMPRESSION);
        n += samBlock->block_length;
    }
    
    n += samBlock->block_length;
    
    //end the compression
    compress_file_size = encoder_last_step(as);
    
    printf("%f Million reads compressed using %f MB.\n", (double)n/1000000.0, (double)compress_file_size/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    return 0;
}


int decompress(FILE *fout, struct remote_file_info info, FILE *fref, struct qv_options_t *qv_opts){
    
    uint64_t n = 0;
    uint32_t i = 0;
    uint8_t reference_flag = 0;
    
    Arithmetic_stream as = alloc_arithmetic_stream(info,ARITHMETIC_WORD_LENGTH, DECOMPRESSION);
    
    sam_block samBlock = alloc_sam_block_t(as, fout, fref, qv_opts, DECOMPRESSION);
    
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
        reset_QV_block(samBlock->QVs, DECOMPRESSION);
        n += samBlock->block_length;
    }
    
    n += samBlock->block_length;
    
    //end the decompression
    //compress_file_size = encoder_last_step(as);
    
    printf("%f Million reads decompressed.\n", (double)n/1000000.0);
    
    // free(samLine->cigar), free(samLine.edits), free(samLine.read_), free(samLine.identifier), free(samLine.refname);
    
    //fclose(fin);
    
    return 0;
}

/**
 * Displays a usage name
 * @param name Program name string
 */
void usage(const char *name) {
    printf("Usage: %s (options) [input file] [output file] [ref file]\n", name);
    printf("Options are:\n");
    //printf("\t-q\t\t\t: Store quality values in compressed file (default)\n");
    printf("\t-x\t\t: Regenerate FASTQ file from compressed file\n");
    printf("\t-c [ratio]\t: Compress using [ratio] bits per bit of input entropy per symbol\n");
    //printf("\t-r [rate]\t: Compress using fixed [rate] bits per symbol\n");
    printf("\t-d [M|L|A]\t: Optimize for MSE, Log(1+L1), L1 distortions, respectively (default: MSE)\n");
    //printf("\t-c [#]\t\t: Compress using [#] clusters (default: 3)\n");
    //printf("\t-u [FILE]\t: Write the uncompressed lossy values to FILE (default: off)\n");
    printf("\t-h\t\t: Print this help\n");
    printf("\t-s\t\t: Print summary stats\n");
    //printf("\t-t [lines]\t: Number of lines to use as training set (0 for all, 1000000 default)\n");
    printf("\t-v\t\t: Enable verbose output\n");
}



int main(int argc, const char * argv[]) {
    
    uint32_t extract, i = 0, file_idx = 0;
    
    struct qv_options_t opts;
    
    char input_name[1024], output_name[1024], ref_name[1024];
    
    char* ptr;
    
    FILE *fref = NULL;
    
    struct remote_file_info remote_info;
    
    clock_t begin = clock();
    
    opts.training_size = 1000000;
    opts.verbose = 0;
    opts.stats = 0;
    opts.ratio = 1;
    opts.uncompressed = 0;
    opts.distortion = DISTORTION_MSE;
    
    extract = COMPRESSION;
    
    // No dependency, cross-platform command line parsing means no getopt
    // So we need to settle for less than optimal flexibility (no combining short opts, maybe that will be added later)
    i = 1;
    while (i < argc) {
        // Handle file names and reject any other untagged arguments
        if (argv[i][0] != '-') {
            switch (file_idx) {
                case 0:
                    strcpy(input_name,argv[i]);
                    file_idx = 1;
                    break;
                case 1:
                    strcpy(output_name,argv[i]);
                    file_idx = 2;
                    break;
                case 2:
                    strcpy(ref_name, argv[i]);
                    file_idx = 3;
                    break;
                default:
                    printf("Garbage argument \"%s\" detected.\n", argv[i]);
                    usage(argv[0]);
                    exit(1);
            }
            i += 1;
            continue;
        }
        
        // Flags for options
        switch(argv[i][1]) {
            case 'x':
                extract = DECOMPRESSION;
                i += 1;
                break;
            //case 'q':
                //extract = COMPRESSION;
                //i += 1;
                //break;
            case 'c':
                extract = COMPRESSION;
                opts.ratio = atof(argv[i+1]);
                opts.mode = MODE_RATIO;
                i += 2;
                break;
            case 'v':
                opts.verbose = 1;
                i += 1;
                break;
            case 'h':
                usage(argv[0]);
                exit(0);
            case 's':
                opts.stats = 1;
                i += 1;
                break;
            //case 't':
                //opts.training_size = atoi(argv[i+1]);
                //i += 2;
                //break;
            case 'd':
                switch (argv[i+1][0]) {
                    case 'M':
                        opts.distortion = DISTORTION_MSE;
                        break;
                    case 'L':
                        opts.distortion = DISTORTION_LORENTZ;
                        break;
                    case 'A':
                        opts.distortion = DISTORTION_MANHATTAN;
                        break;
                    default:
                        printf("Distortion measure not supported, using MSE.\n");
                        break;
                }
                i += 2;
                break;
            default:
                printf("Unrecognized option -%c.\n", argv[i][1]);
                usage(argv[0]);
                exit(1);
        }
    }
    
    if (extract == DECOMPRESSION && file_idx != 3) {
        printf("Missing required filenames in extract mode.\n");
        usage(argv[0]);
        exit(1);
    }
    
    if (extract == COMPRESSION && file_idx != 2) {
        printf("Wrong required filenames in compress mode.\n");
        usage(argv[0]);
        exit(1);
    }
    
    if (opts.verbose) {
        if (extract) {
            printf("%s will be decoded to %s.\n", input_name, output_name);
        }
        else {
            printf("%s will be encoded as %s.\n", input_name, output_name);
            if (opts.mode == MODE_RATIO)
                printf("Ratio mode selected, targeting %f compression ratio\n", opts.ratio);
            else if (opts.mode == MODE_FIXED)
                printf("Fixed-rate mode selected, targeting %f bits per symbol\n", opts.ratio);
            else if (opts.mode == MODE_FIXED_MSE)
                printf("Fixed-MSE mode selected, targeting %f average MSE per context\n", opts.ratio);
            // @todo other modes?
        }
    }
    
    if (extract == COMPRESSION) {
        ptr = strtok((char*)output_name, "@");
        strcpy(remote_info.username, ptr);
        ptr = strtok(NULL, ":");
        strcpy(remote_info.host_name, ptr);
        ptr = strtok(NULL, "\0");
        strcpy(remote_info.filename, ptr);
    }
    
    if (extract == DECOMPRESSION) {
        ptr = strtok((char*)input_name, "@");
        strcpy(remote_info.username, ptr);
        ptr = strtok(NULL, ":");
        strcpy(remote_info.host_name, ptr);
        ptr = strtok(NULL, "\0");
        strcpy(remote_info.filename, ptr);
    }
    
    //FILE * fcomp = (extract == COMPRESSION)? fopen(output_name, "w"):  fopen(input_name, "r");
    
    FILE * fsam = (extract == COMPRESSION)? fopen( input_name, "r"): fopen(output_name, "w");
    
    // Open the Ref file
    if (extract == DECOMPRESSION) {
        if (( fref = fopen ( ref_name , "r" ) ) == NULL){
            fputs ("Chromosome (ref) File error\n",stderr); exit (1);
        }
    }
    
    
    if (extract == COMPRESSION)
        compress(fsam, remote_info, &opts);
    else
        decompress(fsam, remote_info, fref, &opts);
    
    fclose(fsam);
    
    
    
    if (extract == DECOMPRESSION) fclose(fref);
    
    clock_t ticks = clock() - begin;
    
    printf("time: %f\n", ((float)ticks)/CLOCKS_PER_SEC);
#ifdef _WIN32
    system("pause");
#endif
    
    return 0;
    
    
    
    
}

