//
//  read_decompression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/13/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "read_compression.h"


//**************************************************************//
//                                                              //
//                  STORE REFERENCE IN MEMORY                   //
//                                                              //
//**************************************************************//
int store_reference_in_memory(FILE* refFile){
    uint32_t letterCount, endoffile = 1;
    char header[1024];
    
    reference = (char *) malloc(MAX_BP_CHR*sizeof(char));
    
    // ******* Read and Store Reference****** //
    letterCount = 0;
    
    // Remove the header
    fgets(header, sizeof(header), refFile);
    
    while (fgets(&reference[letterCount], 1024, refFile))
    {
        
        if(reference[letterCount] == '>' || reference[letterCount] == '@'){
            endoffile = 0;
            break;
        }
        
        while (reference[letterCount++] != '\n' ) ;
        letterCount--;
        
    }
        
    reference[letterCount] = '\0';
    
    reference = (char *) realloc(reference, letterCount);
    
    if (endoffile)
        return END_GENOME_FLAG;
    
    return letterCount;
    
}


/************************
 * Decompress the read
 **********************/
uint32_t decompress_read(Arithmetic_stream as, sam_block sb, uint8_t chr_change){
    
    int invFlag, tempP;
    
    read_models models = sb->reads->models;
    
    // Decompress the read
    tempP = decompress_pos(as, models->pos, models->pos_alpha, chr_change);
    
    // Check wether we are changing Chromosomes or we are at the end pf the last chromosome
    if (tempP == CHR_CHANGE_FLAG)
        return CHR_CHANGE_FLAG;
    
    else if (tempP == END_GENOME_FLAG)
        return END_GENOME_FLAG;
    
    invFlag = decompress_flag(as, models->flag);
    
    reconstruct_read(as, models, tempP, invFlag, sb->fs);
    
    return invFlag;
}


/***********************
 * Decompress the Flag
 **********************/
uint32_t decompress_flag(Arithmetic_stream a, stream_model *F){
    
    
    // In this case we are just compressing the binary information of whether the
    // read is in reverse or not. we use F[0] as there is no context for the flag.
    
    int isReversed = 0;
    
    // Read the value from the Arithmetic Stream
    isReversed = read_value_from_as(a, F[0]);
    
    // Update model
    update_model(F[0], isReversed);
    
    return isReversed;
    
}

/***********************************
 * Decompress the Alphabet of Position
 ***********************************/
uint32_t decompress_pos_alpha(Arithmetic_stream as, stream_model *PA){
    
    uint32_t Byte = 0, x = 0;
    
    // we encode byte per byte i.e. x = [B0 B1 B2 B3]
    
    // Read B0 from the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[0]);
    // Update model
    update_model(PA[0], Byte);
    // Reconstruct B0
    x |= (Byte << 24);
    
    // Read B1 from the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[1]);
    // Update model
    update_model(PA[1], Byte);
    // Reconstruct B1
    x |= (Byte << 16);
    
    // Send B2 to the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[2]);
    // Update model
    update_model(PA[2], Byte);
    // Reconstruct B2
    x |= (Byte << 8);
    
    // Send B3 to the Arithmetic Stream using the alphabet model
    Byte = read_value_from_as(as, PA[3]);
    // Update model
    update_model(PA[3], Byte);
    // Reconstruct B3
    x |= (Byte);
    
    return x;
    
    
}

/**************************
 * Decompress the Position
 *************************/
uint32_t decompress_pos(Arithmetic_stream as, stream_model *P, stream_model *PA, uint8_t chr_change){
    
    static uint32_t prevPos = 0;
    
    int32_t pos, alphaMapX = 0, x = 0;
    
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    
    // Check if we are changing chromosomes.
    if (chr_change)
        prevPos = 0;
    
    // Read from the AS and get the position
    alphaMapX = read_value_from_as(as, P[0]);
    
    x = P[0]->alphabet[alphaMapX];
    
    // Update the statistics
    update_model(P[0], alphaMapX);
    
    // A new value of pos
    if (x == -1) {
        
        // Read from the AS to get the unknown alphabet letter alpha
        x = decompress_pos_alpha(as, PA);
        
        // Update the statistics of the alphabet for x
        P[0]->alphaExist[x] = 1;
        P[0]->alphaMap[x] = P[0]->alphabetCard; // We reserve the bin 0 for the new symbol flag
        P[0]->alphabet[P[0]->alphabetCard] = x;
        
        update_model(P[0], P[0]->alphabetCard++);
    }
    
    // Decompress the position diference (+ 1 to reserve 0 for new symbols)
    pos = prevPos + x - 1;
    
    prevPos = pos;
    
    return pos;
}

/****************************
 * Decompress the match
 *****************************/
uint32_t decompress_match(Arithmetic_stream a, stream_model *M, uint32_t P){
    
    uint32_t ctx = 0;
    static uint8_t  prevM = 0;
    
    uint8_t match = 0;
    
    
    // Compute Context
    P = (P != 1)? 0:1;
    //prevP = (prevP > READ_LENGTH)? READ_LENGTH:prevP;
    //prevP = (prevP > READ_LENGTH/4)? READ_LENGTH:prevP;
    
    ctx = (P << 1) | prevM;
    
    //ctx = 0;
    
    // Read the value from the Arithmetic Stream
    match = read_value_from_as(a, M[ctx]);
    
    // Update model
    update_model(M[ctx], match);
    
    prevM = match;
    
    return match;
}

/*************************
 * Decompress the snps
 *************************/
uint32_t decompress_snps(Arithmetic_stream a, stream_model *S){
    
    uint8_t numSnps = 0;
    // No context is used for the numSnps for the moment.
    
    // Send the value to the Arithmetic Stream
    numSnps = read_value_from_as(a, S[0]);
    
    // Update model
    update_model(S[0], numSnps);
    
    return numSnps;
    
}


/********************************
 * Decompress the indels
 *******************************/
uint32_t decompress_indels(Arithmetic_stream a, stream_model *I){
    
    uint8_t numIndels = 0;
    // No context is used for the numIndels for the moment.
    
    // Read the value from the Arithmetic Stream
    numIndels = read_value_from_as(a, I[0]);
    
    // Update model
    update_model(I[0], numIndels);
    
    return numIndels;
    
}

/*******************************
 * Decompress the variations
 ********************************/
uint32_t decompress_var(Arithmetic_stream a, stream_model *v,  uint32_t prevPos, uint32_t flag){
    
    uint32_t ctx = 0;
    uint32_t pos = 0;
    
    //flag = 0;
    ctx = prevPos << 1 | flag;
    
    // Read the value from the Arithmetic Stream
    pos = read_value_from_as(a, v[ctx]);
    
    // Update model
    update_model(v[ctx], pos);
    
    return pos;
    
}

/*****************************************
 * Decompress the chars
 ******************************************/
uint8_t decompress_chars(Arithmetic_stream a, stream_model *c, enum BASEPAIR ref){
    
    uint32_t target = 0;
    
    // Read the value from the Arithmetic Stream
    target = read_value_from_as(a, c[ref]);
    
    // Update model
    update_model(c[ref], target);
    
    return basepair2char((enum BASEPAIR)target);
    
}

/*****************************************
 * Reconstruct the read
 ******************************************/
uint32_t reconstruct_read(Arithmetic_stream as, read_models models, uint32_t pos, uint8_t invFlag, FILE *fs){
    
    unsigned int numIns = 0, numDels = 0, numSnps = 0, delPos = 0, ctrPos = 0, snpPos = 0, insPos = 0;
    uint32_t currentPos = 0, prevIns = 0, prev_pos = 0, delta = 0, deltaPos = 0;
    
    static uint32_t prevPos = 0;
    
    unsigned int ctrDels = 0, readCtr = 0;
    int i = 0;
    
    uint8_t match;
    
    enum BASEPAIR refbp;
    
    char *tempRead = (char*)alloca(models->read_length*sizeof(char) + 2);
    char *read = (char*)alloca(models->read_length*sizeof(char) + 4);
    
    read[models->read_length] = '\n';
    read[models->read_length + 1] = '+';
    read[models->read_length + 2] = '\n';
    
    if (pos < prevPos){
        deltaPos = pos;
    }else{
        deltaPos = pos - prevPos + 1;// deltaPos is 1-based.
    }
    prevPos = pos;
    
    // The read matches perfectly.
    match = decompress_match(as, models->match, deltaPos);
    
    // cumsumP is equal to pos
    cumsumP = pos;
    
    // If there is a match, reconstruct the read
    if (match){
        
        switch (invFlag) {
            case 0:
                for (ctrPos=0; ctrPos<models->read_length; ctrPos++)
                    read[readCtr++] = reference[pos + ctrPos - 1];
                break;
                
            case 1:
                for (ctrPos=0; ctrPos<models->read_length; ctrPos++)
                    read[readCtr++] = bp_complement( reference[pos + models->read_length - 1 - ctrPos - 1] );
                
                break;
            default:
                printf("ERROR: invFlag must be 0 or 1\n");
                
        }
        
        fwrite(read, models->read_length + 3, sizeof(char), fs);
        return 1;
    }
    
    // There is no match, retreive the edits
    else{
        numSnps = decompress_snps(as, models->snps);
        
        
        if (numSnps == 0){
            numSnps = decompress_indels(as, models->indels);
            numDels = decompress_indels(as, models->indels);
            numIns = decompress_indels(as, models->indels);
        }
    }
    
    // Reconstruct the read
    
    // Deletions
    prev_pos = 0;
    for (ctrDels = 0; ctrDels < numDels; ctrDels++){
        
        delPos = decompress_var(as, models->var, prev_pos, invFlag);
        prev_pos += delPos;
        
        // Do not take the deleted characters from the reference
        for (ctrPos = 0; ctrPos<delPos; ctrPos++){
            tempRead[currentPos] = reference[pos + currentPos - 1 + ctrDels];
            currentPos++;
        }
    }
    
    
    // Fill up the rest of the read up to numIns
    for (ctrPos = currentPos; ctrPos<models->read_length - numIns; ctrPos++){
        tempRead[currentPos] = reference[pos + currentPos - 1 + numDels];
        currentPos++;
    }
    
    // SNPS
    currentPos = 0;
    prev_pos = 0;
    for (i = 0; i < numSnps; i++){
        
//        assert(currentPos < models->read_length);
        
        // compute delta to next snp
        delta = compute_delta_to_first_snp(prev_pos, models->read_length);
        delta = (delta << BITS_DELTA);
        
        snpPos = decompress_var(as, models->var, delta + prev_pos, invFlag);
        prev_pos += snpPos + 1;
        snpInRef[cumsumP + prev_pos - 1 - 1] = 1;
        
        refbp = char2basepair( tempRead[currentPos + snpPos] );
        tempRead[currentPos + snpPos] = decompress_chars(as, models->chars, refbp);
        currentPos = currentPos + snpPos + 1;
        
    }
    currentPos = 0;
    
    // Insertions and write the read
    switch (invFlag) {
            // There is NO inversion
        case 0:
            // Write the insertions
            prev_pos = 0;
            for (i = 0; i < numIns; i++){
                
                insPos = decompress_var(as, models->var, prev_pos, invFlag);
                prev_pos += insPos;
                
                // write the read up to the insertion
                for (ctrPos=0; ctrPos<insPos; ctrPos++)
                    read[readCtr++] = tempRead[currentPos], currentPos++;
                
                // write the insertion
                read[readCtr++] = decompress_chars(as, models->chars, O);
            }
            
            // write the rest of the read
            for (ctrPos=currentPos; ctrPos < models->read_length - numIns; ctrPos++)
                read[readCtr++] = tempRead[currentPos], currentPos++;
            
            fwrite(read, models->read_length + 3, sizeof(char), fs);
            return 0;
            
            // There is an inversion
        case 1:
            prevIns = 0;
            prev_pos = 0;
            for (i = 0; i < numIns; i++){
                
                insPos = decompress_var(as, models->var, prev_pos, invFlag);
                prev_pos += insPos;
                insPos += prevIns;
                // move the read one position to the left (make room for the insertion)
                for (ctrPos = models->read_length - 1; ctrPos > insPos; ctrPos--)
                    tempRead[ctrPos] = tempRead[ctrPos - 1];
                // Add the insertion
                tempRead[insPos] = decompress_chars(as, models->chars, O);
                prevIns = insPos + 1;
            }
            
            // write the read inverted
            for (ctrPos=0; ctrPos < models->read_length; ctrPos++)
                read[readCtr++] = bp_complement(tempRead[ models->read_length - 1 - ctrPos]);
            
            fwrite(read, models->read_length + 3, sizeof(char), fs);
            return 1;
            
            
        default:
            printf("ERROR: invFLAG different from 0 and 1\n");
            return 2;
    }
}
