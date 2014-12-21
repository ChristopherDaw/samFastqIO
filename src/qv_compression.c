//
//  qv_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/19/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//


#include "util.h"

#include "qv_codebook.h"
#include "sam_block.h"

/*************************
 * Compress the snps
 *************************/
uint32_t compress_qv(Arithmetic_stream a, stream_model model, uint8_t qv){
    
    // The context is given by the index idx.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, model, qv);
    
    // Update model
    update_model(model, qv);
    
    return 1;
    
}

uint32_t decompress_qv(Arithmetic_stream a, stream_model *model, uint32_t idx){
    
    // The context is given by the index idx.
    
    uint32_t qv;
    
    // Send the value to the Arithmetic Stream
    qv = read_value_from_as(a, model[idx]);
    
    // Update model
    update_model(model[idx], qv);
    
    return qv;
    
}

float counts_1 = 0, counts_2 = 0, counts_3 = 0, counts_4 = 0, counts_5 = 0, counts_6 = 0;

double QVs_compress(Arithmetic_stream as, qv_block info, qv_line line) {
    
    
    
    clock_t begin;
    
    uint32_t s = 0, idx = 0, q_state = 0, model_idx = 0;
    //double error = 0.0;
    uint8_t qv = 0, prev_qv = 0;
    uint32_t columns = line->columns;
    
    struct quantizer_t *q;
    
    stream_model mod;
        
    // Select first column's codebook with no left context
    //q = choose_quantizer(info->qlist, &info->well, 0, 0, &idx);
        
    // Quantize, compress and calculate error simultaneously
    // Reading in the lines corrects the quality values to alphabet offsets
    //qv = q->q[line->data[0]];
    //q_state = get_symbol_index(q->output_alphabet, qv);
    
    //compress_qv(as, info->model, get_qv_model_index(0, idx), q_state);
        
    //error = compute_distortion(line->data[0], qv, info->opts->distortion);
        
    //prev_qv = qv;
        
    for (s = 0; s < columns; ++s) {
        
        begin = clock();
        q = choose_quantizer(info->qlist, &info->well, s, prev_qv, &idx);
        counts_1 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        begin = clock();
        qv = q->q[line->data[s]];
        //q_state = get_symbol_index(q->output_alphabet, qv);
        counts_2 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        begin = clock();
        q_state = q->output_alphabet->indexes[qv];
        counts_3 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        begin = clock();
        model_idx = get_qv_model_index(s, idx);
        counts_4 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        //model_idx = 0;
        begin = clock();
        mod = info->model[model_idx];
        counts_5 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        begin = clock();
        //compress_qv(as, info->model[model_idx], q_state);
        compress_qv(as, mod, q_state);
        counts_6 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        //error += compute_distortion(line->data[s], qv, info->opts->distortion);
        prev_qv = qv;
    }
    
    //return error / ((double) columns);
    
    return 0;
}

void foo(){
    
    printf("Time 1 \t %f\n", ((float)counts_1));
    printf("Time 2 \t %f\n", ((float)counts_2));
    printf("Time 3 \t %f\n", ((float)counts_3));
    printf("Time 4 \t %f\n", ((float)counts_4));
    printf("Time 5 \t %f\n", ((float)counts_5));
    printf("Time 5 \t %f\n", ((float)counts_6));
}

double QVs_decompress(Arithmetic_stream as, qv_block info, FILE *fout, uint8_t inv) {
    
    uint32_t idx = 0, q_state = 0;
    uint8_t prev_qv = 0;
    uint32_t columns = info->columns;
    
    int s = 0;
    char *line = (char *) _alloca(columns+1);
    line[columns] = '\n';
    
    struct quantizer_t *q;
    
    // Select first column's codebook with no left context
    q = choose_quantizer(info->qlist, &info->well, 0, 0, &idx);
    
    // Quantize, compress and calculate error simultaneously
    // Reading in the lines corrects the quality values to alphabet offsets
    q_state = decompress_qv(as, info->model, get_qv_model_index(0, idx));
    line[0] = q->output_alphabet->symbols[q_state] + 33;
    prev_qv = line[0] - 33;
    
    for (s = 1; s < columns; ++s) {
        q = choose_quantizer(info->qlist, &info->well, s, prev_qv, &idx);
        q_state = decompress_qv(as, info->model, get_qv_model_index(s, idx));
        line[s] = q->output_alphabet->symbols[q_state] + 33;
        prev_qv = line[s] - 33;
    }
    
    /*if (inv) {
        for (s = columns - 1; s >= 0; s--) {
            fputc(line[s],fout);
        }
        fputc('\n', fout);
        return 1;
    }
    else fwrite(line, 1, columns +1, fout);
    */
    fwrite(line, 1, columns +1, fout);
    return 1;
}

