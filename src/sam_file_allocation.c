//
//  sam_file_allocation.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/18/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "sam_block.h"

void reset_QV_block(qv_block qvb, uint8_t direction){
    
    free_stream_model_qv(qvb->qlist, qvb->model);
    free_cond_quantizer_list(qvb->qlist);
    if (direction == COMPRESSION)
        free_conditional_pmf_list(qvb->training_stats);
    
}

uint32_t get_read_length(FILE *f){
    
    int ch, header_bytes = 0;
    char buffer[2048]; // 2 KB buffer
    
    // We use this opportunity to remove the headers
    while ((ch = getc(f)) == '@') {
        fgets(buffer, 2048, f);
        header_bytes += strlen(buffer) + 1; // +1 to account for the @
    }
    
    // rewind the file pointer to be at the beginning of the first read
    fseek(f, header_bytes, SEEK_SET);
    
    // Extract the first read
    fscanf(f, "%*s %*d %*s %*d %*d %*s %*s %*d %*d %s", buffer);
    
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
    
    rf->block_length = MAX_LINES_PER_BLOCK;
    
    rf->lines = (read_line) calloc(rf->block_length, sizeof(struct read_line_t));
    
    // allocate the memory for each of the lines
    for (i = 0; i < rf->block_length; i++) {
        
        rf->lines[i].cigar = (char*) calloc(1, 2*read_length);
        rf->lines[i].edits = (char*) calloc(1, 2*read_length);
        rf->lines[i].read = (char*) calloc(1, read_length + 3);
    }
    
    // Allocate (and initialize) the models for the reads
    rf->models = alloc_read_models_t(read_length);
    
    return rf;
}

/**
 *
 */
qv_block alloc_qv_block_t(struct qv_options_t *opts, uint32_t read_length){
    
    qv_block qv_info;
    
    uint32_t i = 0;
    
    qv_info = (qv_block) calloc(1, sizeof(struct qv_block_t));
    
    qv_info->block_length = MAX_LINES_PER_BLOCK;
    
    qv_info->columns = read_length;
    
    // Allocate the QV alphabet and the distortion matrix
    struct distortion_t *dist = generate_distortion_matrix(41, opts->distortion);
    struct alphabet_t *alphabet = alloc_alphabet(41);
    
    qv_info->qv_lines = (struct qv_line_t *) calloc(qv_info->block_length, sizeof(struct qv_line_t));
    
    // allocate the memory for each of the lines
    for (i = 0; i < qv_info->block_length; i++) {
        
        qv_info->qv_lines[i].data = (symbol_t*) calloc(qv_info->columns, sizeof(symbol_t));
        qv_info->qv_lines[i].columns = read_length;
    }
    
    // set the alphabet and distortion
    qv_info->alphabet = alphabet;
    qv_info->dist = dist;
    
    qv_info->opts = opts;
    
    return qv_info;
    
}

sam_block alloc_sam_block_t(Arithmetic_stream as, FILE * fin, FILE *fref, struct qv_options_t *qv_opts, uint8_t decompression){
    
    sam_block sb = (sam_block) calloc(1, sizeof(struct sam_block_t));
    
    sb->fs = fin;
    
    sb->fref = fref;
    
    // initialize the codebook_model
    uint32_t rescale = 1 << 20;
    sb->codebook_model = initialize_stream_model_codebook(rescale);
    
    sb->block_length = MAX_LINES_PER_BLOCK;
    
    // Get the Read Length
    if (decompression) {
        // get the readLength from the ios buffer
        sb->read_length =  decompress_int(as, sb->codebook_model);
    }
    else{
        // get the read length from input file and move file pointer after headers
        sb->read_length = get_read_length(sb->fs);
        // write readLength directly to AS using the codebook model
        compress_int(as, sb->codebook_model, sb->read_length);
    }
    
    // Allocate the memory for the three parts:
    //READS,
    sb->reads = alloc_read_block_t(sb->read_length);
    //QVs,
    sb->QVs = alloc_qv_block_t(qv_opts, sb->read_length);
    sb->QVs->codebook_model = sb->codebook_model;
    //@TODO: IDs
    
    return sb;
    
}

uint32_t load_sam_block(sam_block sb){
    
    int32_t i = 0, j = 0;
    read_line rline = NULL;
    qv_line qvline = NULL;
    
    char buffer[1024];
    char *ptr;
    
    for (i = 0; i < sb->block_length; i++) {
        
        rline = &(sb->reads->lines[i]);
        qvline = &(sb->QVs->qv_lines[i]);
        
        rline->read_length = sb->read_length;
        qvline->columns = sb->read_length;
        // Read compulsory fields
        if (fgets(buffer, 1024, sb->fs)) {
            // ID
            ptr = strtok(buffer, "\t");
            // FLAG
            ptr = strtok(NULL, "\t");
            rline->invFlag = atoi(ptr);
            // RNAME
            ptr = strtok(NULL, "\t");
            // POS
            ptr = strtok(NULL, "\t");
            rline->pos = atoi(ptr);
            // MAPQ
            ptr = strtok(NULL, "\t");
            // CIGAR
            ptr = strtok(NULL, "\t");
            strcpy(rline->cigar, ptr);
            // RNEXT
            ptr = strtok(NULL, "\t");
            // PNEXT
            ptr = strtok(NULL, "\t");
            // TLEN
            ptr = strtok(NULL, "\t");
            // SEQ
            ptr = strtok(NULL, "\t");
            strcpy(rline->read, ptr);
            // QUAL
            ptr = strtok(NULL, "\t");
            // Read the QVs and translate them to a 0-based scale
            // Check if the read is inversed
            if (rline->invFlag & 16) { // The read is inverted
                for (j = sb->read_length - 1; j >= 0; j--) {
                    qvline->data[j] = (*ptr) - 33, ptr++;
                }
            }
            else{ // The read is not inversed
                for (j = 0; j < sb->read_length; j++) {
                    qvline->data[j] = (*ptr) - 33, ptr++;
                }
            }
            // Read the AUX fields until end of line, and store the MD field
            while( NULL != (ptr = strtok(NULL, "\t")) ){
                // Do something
                if (*ptr == 'M' && *(ptr+1) == 'D'){
                    // skip MD:Z:
                    ptr += 5;
                    strcpy(rline->edits, ptr);
                    break;
                }
                
            }
        }
        else
            break;
       /* int ch = 0;
        if (EOF!=fscanf(sb->fs, "%*s %"SCNu16" %*s %d %*d %s %*s %*d %*d %s", &(rline->invFlag), &(rline->pos), rline->cigar, rline->read)){
          fgetc(sb->fs); //remove the \t
            
            // Read the QVs and translate them to a 0-based scale
            // Check if the read is inversed
            if (rline->invFlag & 16) { // The read is inversed
                for (j = sb->read_length - 1; j >= 0; j--) {
                    qvline->data[j] = fgetc(sb->fs) - 33;
                }
            }
            else{ // The read is not inversed
                for (j = 0; j < sb->read_length; j++) {
                    qvline->data[j] = fgetc(sb->fs) - 33;
                }
            }
         
        
            // Read the AUX fields until end of line, and store the MD field
            while('\n'!=(ch=fgetc(sb->fs))){
                // Do something
                if (ch == 'M'){
                    ch = fgetc(sb->fs);
                    if (ch == 'D'){
                        // Read :Z:
                        fgetc(sb->fs);// :
                        fgetc(sb->fs);// Z
                        fgetc(sb->fs);// :
                        fscanf(sb->fs, "%s", rline->edits);
                    }
                }
                
            }
            
        }
        else
            break;
    */
        
    }

    // Update the block length
    sb->block_length = i;
    sb->reads->block_length = i;
    sb->QVs->block_length = i;
    
    return i;
}
