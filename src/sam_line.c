//
//  sam_line.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "sam_line.h"

uint8_t read_line_from_sam(struct sam_line_t *samLine, FILE *f){
    
    uint8_t thereIsLine = 0;
    int ch;
    
    // Read compulsory fields
    if (EOF!=fscanf(f, "%*s %"SCNu16" %*s %d %*d %s %*s %*d %*d %s", &samLine->invFlag, &samLine->pos, samLine->cigar, samLine->read)){
        thereIsLine = 1;
        // Read the AUX fields until end of line, and store the MD field
        while('\n'!=(ch=fgetc(f))){
            // Do something
            if (ch == 'M'){
                ch = fgetc(f);
                if (ch == 'D'){
                    // Read :Z:
                    ch = fgetc(f);
                    ch = fgetc(f);
                    ch = fgetc(f);
                    fscanf(f, "%s", samLine->edits);
                }
            }
                
        }
            
    }
    
    return thereIsLine;
    
}

