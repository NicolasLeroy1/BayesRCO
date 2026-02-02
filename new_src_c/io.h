#ifndef IO_H
#define IO_H

#include "bayesRCO.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * Safe memory allocation macro with NULL check and error handling.
 * Exits the program with an error message if allocation fails.
 * 
 * Usage: SAFE_CALLOC(ptr, count, type)
 * Example: SAFE_CALLOC(array, 100, double)
 */
#define SAFE_CALLOC(ptr, count, type) do { \
    (ptr) = (type*)calloc((count), sizeof(type)); \
    if ((ptr) == NULL && (count) > 0) { \
        fprintf(stderr, "Error: Memory allocation failed for %s (%zu bytes) at %s:%d\n", \
                #ptr, (size_t)(count) * sizeof(type), __FILE__, __LINE__); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

/**
 * Safe memory allocation macro for malloc with NULL check.
 */
#define SAFE_MALLOC(ptr, count, type) do { \
    (ptr) = (type*)malloc((count) * sizeof(type)); \
    if ((ptr) == NULL && (count) > 0) { \
        fprintf(stderr, "Error: Memory allocation failed for %s (%zu bytes) at %s:%d\n", \
                #ptr, (size_t)(count) * sizeof(type), __FILE__, __LINE__); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

/**
 * Safe free macro that also nullifies the pointer.
 */
#define SAFE_FREE(ptr) do { \
    if ((ptr) != NULL) { \
        free(ptr); \
        (ptr) = NULL; \
    } \
} while(0)

void get_size(ModelConfig *config, GenomicData *gdata);
void load_phenos_plink(ModelConfig *config, GenomicData *gdata);
void load_snp_binary(ModelConfig *config, GenomicData *gdata);
void init_random_seed_custom(ModelConfig *config, prng_state *rs); // Renamed to avoid confusion
void allocate_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
void load_param(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore, MCMCState *mstate);
void load_categories(ModelConfig *config, GenomicData *gdata);
void write_dgv(ModelConfig *config, GenomicData *gdata);
void output_model(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore);
void output_beta(ModelConfig *config, MCMCState *mstate, GenomicData *gdata);

/**
 * Free all dynamically allocated memory in data structures.
 * Should be called before program exit to prevent memory leaks.
 */
void cleanup_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);

#endif // IO_H
