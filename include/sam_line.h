//
//  sam_line.h
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef XC_s2fastqIO_sam_line_h
#define XC_s2fastqIO_sam_line_h

#include <stdio.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#define MAX_READ_LENGTH 1024

typedef struct sam_line_t{
    char *identifier;
    char *refname;
    char *cigar;
    char *edits;
    char *read;
    uint16_t invFlag;
    // int invFlag;
    int pos;
}*sam_line;

// Function Prototypes
sam_line alloc_sam_line_t();
uint8_t read_line_from_sam(struct sam_line_t *samLine, FILE *f);

#endif
