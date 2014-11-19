//
//  qv_compression.c
//  XC_s2fastqIO
//
//  Created by Mikel Hernaez on 11/19/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "codebook.h"


void qv_compress(char *input_name, char *output_name, struct qv_options_t *opts) {
    
    struct qv_block_t qv_info;
    struct distortion_t *dist = generate_distortion_matrix(41, opts->distortion);
    struct alphabet_t *alphabet = alloc_alphabet(41);
    uint32_t status;
    struct hrtimer_t cluster_time, stats, encoding, total;
    FILE *fout, *funcompressed = NULL;
    uint64_t bytes_used;
    double distortion;
    
    start_timer(&total);
    
    qv_info.alphabet = alphabet;
    qv_info.dist = dist;
    
    // Load input file all at once
    status = load_file(input_name, &qv_info, 0);
    if (status != LF_ERROR_NONE) {
        printf("load_file returned error: %d\n", status);
        exit(1);
    }
    
    // Set up clustering data structures
    qv_info.opts = opts;
    
    // Then find stats and generate codebooks
    start_timer(&stats);
    calculate_statistics(&qv_info);
    generate_codebooks(&qv_info);
    stop_timer(&stats);
    
    if (opts->verbose) {
        printf("Stats and codebook generation took %.4f seconds\n", get_timer_interval(&stats));
        // @todo expected distortion is inaccurate due to lack of pmf
        //printf("Expected distortion: %f\n", opts->e_dist);
    }
    
    // Note that we want \r\n translation in the input
    // but we do not want it in the output
    fout = fopen(output_name, "wb");
    if (!fout) {
        perror("Unable to open output file");
        exit(1);
    }
    
    if (opts->uncompressed) {
        funcompressed = fopen(opts->uncompressed_name, "w");
        if (!funcompressed) {
            perror("Unable to open uncompressed file");
            exit(1);
        }
    }
    
    // @todo qv_compression should use quality_file structure with data in memory, now
    start_timer(&encoding);
    write_codebooks(fout, &qv_info);
    bytes_used = start_qv_compression(&qv_info, fout, &distortion, opts->distortion, funcompressed);
    stop_timer(&encoding);
    stop_timer(&total);
    
    fclose(fout);
    
    // Verbose stats
    if (opts->verbose) {
        // @todo add cluster info here
        switch (opts->distortion) {
            case DISTORTION_MANHATTAN:
                printf("L1 distortion: %f\n", distortion);
                break;
            case DISTORTION_MSE:
                printf("MSE distortion: %f\n", distortion);
                break;
            case DISTORTION_LORENTZ:
                printf("log(1+L1) distortion: %f\n", distortion);
                break;
            default:
                break;
        }
        printf("Lines: %llu\n", qv_info.lines);
        printf("Columns: %u\n", qv_info.columns);
        printf("Total bytes used: %llu\n", bytes_used);
        printf("Encoding took %.4f seconds.\n", get_timer_interval(&total));
        printf("Total time elapsed: %.4f seconds.\n", get_timer_interval(&total));
    }
    
    // Parse-able stats
    if (opts->stats) {
        printf("rate, %.4f, distortion, %.4f, time, %.4f, size, %llu \n", (bytes_used*8.)/((double)(qv_info.lines)*qv_info.columns), distortion, get_timer_interval(&total), bytes_used);
    }
}
