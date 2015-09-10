//
//  aux_compression.h
//  samFastqIO
//
//  Created by Carlos Navarro Astiasarán on 29/08/15.
//  Copyright (c) 2015 Carlos Navarro Astiasarán. All rights reserved.
//

#ifndef samFastqIO_aux_compression_h
#define samFastqIO_aux_compression_h

#include <stdio.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#define TAGTYPELUTLENGTH 43
#define TYPELUTLENGTH 6
#define NOTFOUNDLUTFLAG 0x8000
#define MASKTAGCHAR1 0x7E00
#define MASKTAGCHAR2 0x01F8
#define MASKTYPE 0x0007



/*
 // Aux. format: TAG:TYPE:VALUE
 
 //     tag: [a-zA-Z][a-zA-Z0-9],
 //    type: [A,i,f,Z,H,B],
 */


//lookup table corresponding to the most used tag-type pairs given in the specification (SAMv1.pdf)
const char tagTypeLUT[TAGTYPELUTLENGTH][4] = {"AMi","ASi","BCZ","BQZ","CCZ","CMi","COZ","CPi","CQZ","CSZ","CTZ","E2Z","FIi","FSZ","FZB","LBZ","H0i","H1i","H2i","HIi","IHi","MCZ","MDZ","MQi","NHi","NMi","OQZ","OPi","OCZ","PGZ","PQi","PTZ","PUZ","QTZ","Q2Z","R2Z","RGZ","RTZ","SAZ","SMi","TCi","U2Z","UQi"};

//available types
const char typeLUT[]="AifZHB";

uint16_t preprocessTagType(char *ptr);
uint8_t charMap(char c);
char inverseCharMap(uint8_t c);
uint8_t typeMap(char t);




#endif
