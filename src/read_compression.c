//
//  reads_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/5/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//


#include "read_compression.h"

/************************
 * Compress the read
 **********************/
uint32_t compress_read(Arithmetic_stream as, read_models models, read_line samLine){
    
    int tempF, PosDiff, chrPos;
    // Compress sam line
    PosDiff = compress_pos(as, models->pos, models->pos_alpha, samLine->pos);
    //tempF = compress_flag(as, models->flag, samLine->invFlag);
    tempF = compress_flag(as, models->flag, 0);
    chrPos = compress_edits(as, models, samLine->edits, samLine->cigar, samLine->read, PosDiff, tempF);
    
    assert(samLine->pos  == chrPos);

    return 1;
}


/***********************
 * Compress the Flag
 **********************/
uint32_t compress_flag(Arithmetic_stream a, stream_model *F, uint16_t flag){
    
    
    // In this case we are just compressing the binary information of whether the
    // read is in reverse or not. we use F[0] as there is no context for the flag.
    
    uint16_t x = 0;
    
    x = flag << 11;
    x >>= 15;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, F[0], x);
    
    // Update model
    update_model(F[0], x);
    
    return x;
    
}

/***********************************
 * Compress the Alphabet of Position
 ***********************************/
uint32_t compress_pos_alpha(Arithmetic_stream as, stream_model *PA, uint32_t x){
    
    uint32_t Byte = 0;
    
    // we encode byte per byte i.e. x = [B0 B1 B2 B3]
    
    // Send B0 to the Arithmetic Stream using the alphabet model
    Byte = x >> 24;
    send_value_to_as(as, PA[0], Byte);
    // Update model
    update_model(PA[0], Byte);
    
    // Send B1 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x00ff0000) >> 16;
    send_value_to_as(as, PA[1], Byte);
    // Update model
    update_model(PA[1], Byte);
    
    // Send B2 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x0000ff00) >> 8;
    send_value_to_as(as, PA[2], Byte);
    // Update model
    update_model(PA[2], Byte);
    
    // Send B3 to the Arithmetic Stream using the alphabet model
    Byte = (x & 0x000000ff);
    send_value_to_as(as, PA[3], Byte);
    // Update model
    update_model(PA[3], Byte);
    
    return 1;
    
    
}

/***********************
 * Compress the Position
 **********************/
uint32_t compress_pos(Arithmetic_stream as, stream_model *P, stream_model *PA, uint32_t pos){
    
    static uint32_t prevPos = 0;
    enum {SMALL_STEP = 0, BIG_STEP = 1};
    int i = 0;
    int32_t x = 0;
    
    // TODO diferent update models for updating -1 and already seen symbols
    // i.e., SMALL_STEP and BIG_STEP
    
    // Check if we are changing chromosomes.
    if (pos < prevPos) {
        
        cumsumP = 0;
        for (i=0; i<MAX_BP_CHR; i++){
            snpInRef[i] = 0;
        }
        
        // Send 0 to the Arithmetic Stream
        send_value_to_as(as, P[0], 0);
        
        // Update model
        update_model(P[0], 0);
        
        // Send CHR_CHANGE_FLAG to the Arithmetic Stream using the alphabet model
        compress_pos_alpha(as, PA, CHR_CHANGE_FLAG);
        
        // Compress the position
        x = pos;
        if (P[0]->alphaExist[x]){
            // Send x to the Arithmetic Stream
            send_value_to_as(as, P[0], P[0]->alphaMap[x]);
            // Update model
            update_model(P[0], P[0]->alphaMap[x]);
        }
        else{
            
            // Send 0 to the Arithmetic Stream
            send_value_to_as(as, P[0], 0);
            
            // Update model
            update_model(P[0], 0);
            
            // Send the new letter to the Arithmetic Stream using the alphabet model
            compress_pos_alpha(as, PA, x);
            
            // Update the statistics of the alphabet for x
            P[0]->alphaExist[x] = 1;
            P[0]->alphaMap[x] = P[0]->alphabetCard; // We reserve the bin 0 for the new symbol flag
            P[0]->alphabet[P[0]->alphabetCard] = x;
            
            
            // Update model
            update_model(P[0], P[0]->alphabetCard++);
        }
        
        
        prevPos = pos;
        
        return x;
    }
    
    
    // Compress the position diference (+ 1 to reserve 0 for new symbols)
    x = pos - prevPos + 1;
    
    if (P[0]->alphaExist[x]){
        // Send x to the Arithmetic Stream
        send_value_to_as(as, P[0], P[0]->alphaMap[x]);
        // Update model
        update_model(P[0], P[0]->alphaMap[x]);
    }
    else{
        
        // Send 0 to the Arithmetic Stream
        send_value_to_as(as, P[0], 0);
        
        // Update model
        update_model(P[0], 0);
        
        // Send the new letter to the Arithmetic Stream using the alphabet model
        compress_pos_alpha(as, PA, x);
        
        // Update the statistics of the alphabet for x
        P[0]->alphaExist[x] = 1;
        P[0]->alphaMap[x] = P[0]->alphabetCard; // We reserve the bin 0 for the new symbol flag
        P[0]->alphabet[P[0]->alphabetCard] = x;
        
        // Update model
        update_model(P[0], P[0]->alphabetCard++);
    }
    
    prevPos = pos;
    
    return x;
}

/****************************
 * Compress the match
 *****************************/
uint32_t compress_match(Arithmetic_stream a, stream_model *M, uint8_t match, uint32_t P){
    
    uint32_t ctx = 0;
    static uint8_t  prevM = 0;
    
    
    // Compute Context
    P = (P != 1)? 0:1;
    //prevP = (prevP > READ_LENGTH)? READ_LENGTH:prevP;
    //prevP = (prevP > READ_LENGTH/4)? READ_LENGTH:prevP;
    
    ctx = (P << 1) | prevM;
    
    //ctx = 0;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, M[ctx], match);
    
    // Update model
    update_model(M[ctx], match);
    
    prevM = match;
    
    return 1;
}

/*************************
 * Compress the snps
 *************************/
uint32_t compress_snps(Arithmetic_stream a, stream_model *S, uint8_t numSnps){
    
    
    // No context is used for the numSnps for the moment.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, S[0], numSnps);
    
    // Update model
    update_model(S[0], numSnps);
    
    return 1;
    
}


/********************************
 * Compress the indels
 *******************************/
uint32_t compress_indels(Arithmetic_stream a, stream_model *I, uint8_t numIndels){
    
    
    // Nos context is used for the numIndels for the moment.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, I[0], numIndels);
    
    // Update model
    update_model(I[0], numIndels);
    
    return 1;
    
}

/*******************************
 * Compress the variations
 ********************************/
uint32_t compress_var(Arithmetic_stream a, stream_model *v, uint32_t pos, uint32_t prevPos, uint32_t flag){
    
    uint32_t ctx = 0;
    
    //flag = 0;
    ctx = prevPos << 1 | flag;
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, v[ctx], pos);
    
    // Update model
    update_model(v[ctx], pos);
    
    return 1;
    
}

/*****************************************
 * Compress the chars
 ******************************************/
uint32_t compress_chars(Arithmetic_stream a, stream_model *c, enum BASEPAIR ref, enum BASEPAIR target){
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, c[ref], target);
    
    // Update model
    update_model(c[ref], target);
    
    return 1;
    
}

/*****************************************
 * Compress the edits
 ******************************************/
uint32_t compress_edits(Arithmetic_stream as, read_models rs, char *edits, char *cigar, char *read, uint32_t deltaP, uint8_t flag){
    
    unsigned int numIns = 0, numDels = 0, numSnps = 0, lastSnp = 1;
    int i = 0, M = 0, I = 0, D = 0, pos = 0, ctr = 0, prevPosI = 0, prevPosD = 0, ctrS = 0, S = 0;
    uint32_t delta = 0;
    
    // pos in the reference
    cumsumP = cumsumP + deltaP - 1;// DeltaP is 1-based
    
    uint32_t Dels[MAX_READ_LENGTH];
    ins Insers[MAX_READ_LENGTH];
    snp SNPs[MAX_READ_LENGTH];
    
    uint8_t firstCase = 1;
    
    uint32_t prev_pos = 0;
    
    if(strcmp(edits, rs->_readLength) == 0){
        // The read matches perfectly.
        compress_match(as, rs->match, 1, deltaP);
        return cumsumP;
    }
    
    compress_match(as, rs->match, 0, deltaP);
    
    // The read does not match perfectly
    // Compute the edits
    while (*cigar != 0){
        if ( isdigit( *(cigar + i) ) == 0 )
            switch ( *(cigar + i) ){
                    // compute the position of the edit
                case 'M':
                    M += atoi(cigar);
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                    
                    // Store a insertion and all the previous snps
                case 'I':
                    I = atoi(cigar);
                    for (ctr = 0; ctr < I ; ctr++) {
                        pos = M;
                        if (lastSnp != 0)
                            lastSnp = add_snps_to_array(edits, SNPs, &numSnps, pos + numIns, read);
                        
                        Insers[numIns].pos = pos - prevPosI;
                        Insers[numIns].targetChar = char2basepair(read[pos+numIns]);
                        prevPosI = pos;
                        numIns++;
                    }
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                    
                    // Store the deletion
                case 'D':
                    D = atoi(cigar);
                    for (ctr = 0; ctr < D ; ctr++) {
                        pos = M;
                        Dels[numDels] = pos - prevPosD;
                        prevPosD = pos;
                        numDels++;
                    }
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                    
                case '*':
                    return 1;
                case 'S':
                    if (firstCase == 1) {
                        // S is the first thing we see
                        S = atoi(cigar);
                        for (ctrS = 0; ctrS < S; ctrS++){
                            if (lastSnp != 0)
                                lastSnp = add_snps_to_array(edits, SNPs, &numSnps, numIns, read);
                            Insers[numIns].pos = 0;
                            Insers[numIns].targetChar = char2basepair(read[ctrS]);
                            numIns++;
                        }
                    }else{
                        // We are at the end
                        S = atoi(cigar);
                        for (ctr = 0; ctr < S ; ctr++) {
                            pos = M;
                            Insers[numIns].pos = pos - prevPosI;
                            Insers[numIns].targetChar = char2basepair(read[pos+numIns]);
                            prevPosI = pos;
                            numIns++;
                        }
                    }
                    
                    firstCase = 0;
                    cigar = cigar + i + 1;
                    i = -1;
                    break;
                default:
                    break;
                    //printf("Something besides MIDS appeared in the cigar\n");
                    //printf("%c\n", *(cigar + i));
            }
        i++;
    }
    
    if (lastSnp != 0)
        add_snps_to_array(edits, SNPs, &numSnps, rs->read_length + 1, read);
    
    
    // Compress the edits
    if ((numDels | numIns) == 0) {
        compress_snps(as, rs->snps, numSnps);
    }
    else{
        compress_snps(as, rs->snps, 0);
        compress_indels(as, rs->indels, numSnps);
        compress_indels(as, rs->indels, numDels);
        compress_indels(as, rs->indels, numIns);
    }
    
    // Store the positions and Chars in the corresponding vector
    prev_pos = 0;
    for (i = 0; i < numDels; i++){
        compress_var(as, rs->var, Dels[i], prev_pos, flag);
        prev_pos += Dels[i];
    }
    prev_pos = 0;
    for (i = 0; i < numSnps; i++){
        
        // compute delta to next snp
        delta = compute_delta_to_first_snp(prev_pos, rs->read_length);
        /*delta = READ_LENGTH + 2;
         for (j=0;j<READ_LENGTH - prev_pos; j++){
         if (snpInRef[cumsumP - 1 + j] == 1){
         delta = j;
         break;
         }
         }*/
        
        delta = (delta << BITS_DELTA);
        compress_var(as, rs->var, SNPs[i].pos, delta + prev_pos, flag);
        prev_pos += SNPs[i].pos + 1;
        snpInRef[cumsumP + prev_pos - 1 - 1] = 1;
        
        compress_chars(as, rs->chars, SNPs[i].refChar, SNPs[i].targetChar);
        
    }
    prev_pos = 0;
    for (i = 0; i < numIns; i++){
        compress_var(as, rs->var, Insers[i].pos, prev_pos, flag);
        prev_pos += Insers[i].pos;
        
        compress_chars(as, rs->chars, O, Insers[i].targetChar);
    }
    
    return cumsumP;
    
}



/******************************************
 * Function to look for snps in the cigar
 ****************************************/
int add_snps_to_array(char* edits, snp* SNPs, unsigned int *numSnps, unsigned int insertionPos, char *read){
    
    static unsigned int prevEditPtr = 0, cumPos = 0;
    
    int pos = 0, tempPos = 0, ctr;
    char ch = 0;
    
    uint8_t flag = 0;
    
    edits += prevEditPtr;
    
    while (*edits != 0 ) {
        
        pos = atoi(edits);
        
        tempPos = pos;
        
        ctr = compute_num_digits(pos);
        ch = *(edits+ctr);
        ctr++;
        
        // if there are deletions after pos, we need to add those positions that come after the deletions
        while (ch == '^'){
            while (isnumber(*(edits+ctr)) == 0) {
                ctr++;
            }
            tempPos += atoi(edits + ctr);
            
            ctr += compute_num_digits(atoi(edits + ctr));
            
            ch = *(edits+ctr);
            ctr++;
            if (ch == '\0'){
                flag = 1;
                break;
            }
        }
        
        if (flag == 1){
            flag = 0;
            break;
        }
        
        
        if (cumPos + tempPos >= insertionPos){
            cumPos++;
            return cumPos;
        }
        
        tempPos = atoi(edits);
        edits += compute_num_digits(tempPos);
        prevEditPtr += compute_num_digits(tempPos);
        
        //if ( isnumber(*(edits+1)) ) edits += 2, prevEditPtr += 2;
        //else edits += 1, prevEditPtr += 1;
        
        ch = *edits++, prevEditPtr++;
        
        while (ch == '^'){
            while (isnumber(*edits) == 0)
                edits++, prevEditPtr++;
            pos += atoi(edits);
            
            tempPos = atoi(edits);
            edits += compute_num_digits(tempPos);
            prevEditPtr += compute_num_digits(tempPos);
            
            //if ( isnumber(*(edits+1)) ) edits += 2, prevEditPtr += 2;
            //else edits += 1, prevEditPtr += 1;
            ch = *edits++, prevEditPtr++;
        }
        
        if (ch == '\0')
            break;
        
        cumPos += pos;
        
        SNPs[*numSnps].pos = pos;
        SNPs[*numSnps].refChar = char2basepair(ch);
        SNPs[*numSnps].targetChar = char2basepair(read[cumPos]);
        (*numSnps)++;
        cumPos++;
        if (*edits == 0)
            break;
    }
    
    prevEditPtr = 0;
    cumPos = 0;
    return 0;
}

uint32_t compute_delta_to_first_snp(uint32_t prevPos, uint32_t readLen){
    
    uint32_t deltaOut;
    uint32_t j = 0;
    
    deltaOut = readLen + 2;
    
    for (j=0;j<readLen - prevPos; j++){
        if (snpInRef[cumsumP - 1 + j + prevPos] == 1){
            deltaOut = j;
            break;
        }
    }
    
    return deltaOut;
}

uint32_t compute_num_digits(uint32_t x){
    
    //Get the number of digits (We assume readLength < 1000)
    
    if (x < 10)
        return 1;
    else if (x < 100)
        return 2;
    else
        return 3;
    
}
