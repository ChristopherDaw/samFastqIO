//
//  io_functions.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 12/16/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "Arithmetic_stream.h"

int clean_compressed_dir(struct io_stream_t* ios){
    
    uint32_t rc=0;
    ios->fileCtr = 0;
    
    do{
        sprintf(ios->filePath, IDOFILE_PATH_ROOT "%010d", ios->fileCtr);
        rc = remove(ios->filePath);
        ios->fileCtr++;
    }while (rc == 0);
    
    ios->fileCtr = 0;
    
    return rc;
    
}

void open_new_iofile(struct io_stream_t* ios){
    
    sprintf(ios->filePath, IDOFILE_PATH_ROOT "%010d", ios->fileCtr);
    ios->fileCtr++;
    
    if (ios->direcction == DECOMPRESSION) {
        while (file_available == 0) ;
        ios->fp = fopen(ios->filePath, "r");
        fread(ios->buf, sizeof(uint8_t), IO_STREAM_BUF_LEN, ios->fp);
        fclose(ios->fp);
//        remove(ios->filePath);
        file_available--;
        
    }
    else ios->fp = fopen(ios->filePath, "w");
}

/**
 * Writes out the current stream buffer regardless of fill amount
 */
void stream_write_buffer(struct io_stream_t *os) {
    
    open_new_iofile(os);
    fwrite(os->buf, sizeof(uint8_t), os->bufPos, os->fp);
    fclose(os->fp);
    file_available++;
    memset(os->buf, 0, sizeof(uint8_t)*(os->bufPos));
    os->written += os->bufPos;
    os->bufPos = 0;
}

/**
 * Fills out the current stream buffer regardless of fill amount
 */
void stream_fill_buffer(struct io_stream_t *is) {
    
    fclose(is->fp);
    open_new_iofile(is);
    //fread(is->buf, sizeof(uint8_t), IO_STREAM_BUF_LEN, is->fp);
    is->bufPos = 0;
}