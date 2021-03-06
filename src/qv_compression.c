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
uint32_t compress_qv(Arithmetic_stream a, stream_model *model, uint32_t idx, uint8_t qv){
    
    // The context is given by the index idx.
    
    // Send the value to the Arithmetic Stream
    send_value_to_as(a, model[idx], qv);
    
    // Update model
    update_model(model[idx], qv);
    
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

//float counts_1 = 0, counts_2 = 0, counts_3 = 0, counts_4 = 0, counts_5 = 0, counts_6 = 0;

double QVs_compress2(Arithmetic_stream as, qv_block info) {
    
    uint32_t s = 0, idx = 0, q_state = 0;
    //double error = 0.0;
    uint8_t prev_qv = 0, currentQV = 0, quantizer_type, lossyqv;
    uint32_t columns = info->columns;
    
    struct quantizer_t *q;
    
    struct qv_line_t qline = *info->qv_lines;
    
    //stream_model mod;
    
    for (s = 0; s < columns; ++s) {
        
        //begin = clock();
        q = choose_quantizer(info->qlist, &info->well, s, prev_qv, &idx, &quantizer_type);
        
        currentQV = qline.data[s];
        
//        idx1 = (s*3362 + quantizer_type*1681 + prev_qv*41 + currentQV) << 1;
        
//        qv1 = qArray[idx1];
        
//        q_state1 = qArray[idx1+1];
        
        lossyqv = q->q[currentQV];
        
        q_state = q->output_alphabet->indexes[lossyqv];
        
//        assert(q_state == q_state1);
        
        //mod = model[get_qv_model_index(s, idx)];
        //begin = clock();
        //compress_qv(as, info->model[model_idx], q_state);
//        q_state1 = currentQV;
//        q_state1 = lossyqv;
        compress_qv(as, info->model, get_qv_model_index(s, idx), q_state);
        
//        compress_qv(as, info->model, get_qv_model_index(0, 0), q_state);
        
        prev_qv = lossyqv;
        
//        idx1 = currentQV;
    }

    
    return 0;
}


double QVs_compress(Arithmetic_stream as, qv_block info, symbol_t *qArray) {
    
    
    uint32_t s = 0, idx = 0, q_state1 = 0, idx1 = 0;
    uint32_t a_head = QV_ALPHABET_SIZE;
    uint32_t q_head = QV_ALPHABET_SIZE*QV_ALPHABET_SIZE;
    uint32_t s_head = q_head<<1;
    //double error = 0.0;
    uint8_t prev_qv = 0, currentQV = 0, quantizer_type, qv1;
    uint32_t columns = info->columns;
    
    struct qv_line_t qline = *info->qv_lines;
    struct alphabet_t **inAlpha = info->qlist->input_alphabets;
    uint8_t **ratios = info->qlist->qratio;
        
    for (s = 0; s < columns; ++s) {
        
        //begin = clock();
        //q = choose_quantizer(info->qlist, &info->well, s, prev_qv, &idx, &quantizer_type);
        
        //idx = get_symbol_index(info->qlist->input_alphabets[s], prev_qv);
        //idx = info->qlist->input_alphabets[s]->indexes[prev_qv];
        idx = inAlpha[s]->indexes[prev_qv];
//        assert(idx != ALPHABET_SYMBOL_NOT_FOUND);
        //if (well_1024a_bits(&info->well, 7) >= info->qlist->qratio[s][idx]) {
        if (well_1024a_bits(&info->well, 7) >= ratios[s][idx]) {
            idx = 2*idx+1;
            quantizer_type = 1;
        }
        else{
            quantizer_type = 0;
            idx = 2*idx;
        }
       //counts_1 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        //begin = clock();
        currentQV = qline.data[s];
        //qv = q->q[currentQV];
        
        //q_state = get_symbol_index(q->output_alphabet, qv);
        
        //counts_2 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        //begin = clock();
        //q_state = q->output_alphabet->indexes[qv];
        //counts_3 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        
        
        idx1 = (s*s_head + quantizer_type*q_head + prev_qv*a_head + currentQV) << 1;
        
        qv1 = qArray[idx1++];
        
        q_state1 = qArray[idx1];
        
        //assert(idx == idx1);
        
        

        //assert(q_state == q_state1);
        //assert(qv1 == currentQV);
        
        //model_idx = 0;
        
        //mod = info->model[get_qv_model_index(s, idx)];
        //begin = clock();
        //compress_qv(as, info->model[model_idx], q_state);
        //compress_qv(as, mod, q_state1);
        //counts_6 += (float)(clock() - begin)/CLOCKS_PER_SEC;
        compress_qv(as, info->model, get_qv_model_index(s, idx), q_state1);
        //error += compute_distortion(line->data[s], qv, info->opts->distortion);
        prev_qv = qv1;
    }
    
    //return error / ((double) columns);
    
    return 0;
}

double QVs_compress_lossless(Arithmetic_stream as, stream_model* models, qv_line line) {
    
    uint32_t s = 0;
    uint8_t q_state = 0, prev_state = 0;
    uint32_t columns = line->columns;
    
    for (s = 0; s < columns; ++s) {
        
        q_state = line->data[s];
        compress_qv(as, models, get_qv_model_index(s, prev_state), q_state);
        prev_state = q_state;
    }
    
    return 0;
}

double QVs_decompress_lossless(Arithmetic_stream as, qv_block info, uint8_t inv, char* quals) {
    
    uint8_t prev_qv = 0;
    uint32_t columns = info->columns;
    
    int s = 0;
    //char *line = (char *) _alloca(columns+1);
    quals[columns] = 0;
    
    for (s = 0; s < columns; ++s) {
        quals[s] = decompress_qv(as, info->model, get_qv_model_index(s, prev_qv));
        prev_qv = quals[s];
        quals[s] += 33;
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
    //fwrite(line, 1, columns +1, fout);
    return 1;
}

double QVs_decompress(Arithmetic_stream as, qv_block info, uint8_t inv, char* quals) {
    
    uint32_t idx = 0, q_state = 0;
    uint8_t prev_qv = 0, foo;
    uint32_t columns = info->columns;
    
    int s = 0;
    //char *line = (char *) _alloca(columns+1);
    quals[columns] = 0;
    
    struct quantizer_t *q;
    
    // Select first column's codebook with no left context
    q = choose_quantizer(info->qlist, &info->well, 0, 0, &idx, &foo);
    
    // Quantize, compress and calculate error simultaneously
    // Reading in the lines corrects the quality values to alphabet offsets
    q_state = decompress_qv(as, info->model, get_qv_model_index(0, idx));
    quals[0] = q->output_alphabet->symbols[q_state] + 33;
    prev_qv = quals[0] - 33;
    
    for (s = 1; s < columns; ++s) {
        q = choose_quantizer(info->qlist, &info->well, s, prev_qv, &idx, &foo);
        q_state = decompress_qv(as, info->model, get_qv_model_index(s, idx));
        quals[s] = q->output_alphabet->symbols[q_state] + 33;
        prev_qv = quals[s] - 33;
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
    //fwrite(line, 1, columns +1, fout);
    return 1;
}

