//
//  sam_line.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "sam_line.h"

sam_line alloc_sam_line_t(){
    
    sam_line samLine;
    
    samLine = (sam_line) calloc(1, sizeof(struct sam_line_t));
    samLine->identifier = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->refname = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->cigar = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->edits = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->read = (char*) calloc(1, 4*MAX_READ_LENGTH);
    samLine->QV = (char*) calloc(1, 4*MAX_READ_LENGTH);
    
    return samLine;
    
}


uint8_t read_line_from_sam(struct sam_line_t *samLine, FILE *f){
    
    static uint8_t firstCall = 1;
    uint8_t thereIsLine = 0;
    char *token;
    int ch;
    
    if (firstCall == 1){
        char line[100000] = {0};
        
        while ( fgets(line, 100000, f) != NULL) {
            if (line[0] == '@'){
                // Line is a header
                
            }else{
                firstCall = 0;
                // Line belongs to a read
                thereIsLine = 1;
                
                // Identifier
                token = strtok(line, "\t");
                strcpy(samLine->identifier, token);
                
                // Flag
                token = strtok(NULL, "\t");
                samLine->invFlag = atoi(token);
                
                // Refname
                token = strtok(NULL, "\t");
                strcpy(samLine->refname, token);
                
                // Pos
                token = strtok(NULL, "\t");
                samLine->pos = atoi(token);
                
                // MapQ
                token = strtok(NULL, "\t");
                
                // Cigar
                token = strtok(NULL, "\t");
                strcpy(samLine->cigar, token);
                
                // Rnext, Pnext, Tlen
                token = strtok(NULL, "\t");
                token = strtok(NULL, "\t");
                token = strtok(NULL, "\t");
                
                // Read
                token = strtok(NULL, "\t");
                strcpy(samLine->read, token);
                
                // QV sequence
                token = strtok(NULL, "\t");
                strcpy(samLine->QV, token);
                
                /* walk through AUX fields */
                samLine->edits[0] = '$';
                token = strtok(NULL, "\t");
                while( token != NULL )
                {
                    if (strncmp(token, "MD:Z:", 5) == 0) {
                        strcpy(samLine->edits, token + 5);
                    }
                    token = strtok(NULL, "\t");
                }
                break;
            }
        }
    }else{
        // firstCall = 0
        // Read compulsory fields
        if (EOF!=fscanf(f, "%s %"SCNu16" %s %d %*d %s %*s %*d %*d %s %s", samLine->identifier, &samLine->invFlag, samLine->refname, &samLine->pos, samLine->cigar, samLine->read, samLine->QV)){
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
        
        
    }
    
    return thereIsLine;
    
}
