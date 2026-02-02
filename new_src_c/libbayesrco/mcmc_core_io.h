#ifndef MCMC_CORE_IO_H
#define MCMC_CORE_IO_H

#include "bayesRCO.h"
#include <stdlib.h>
#include <stdio.h>

/**
 * Safe memory allocation macro with NULL check.
 * In a library, we should ideally return ERR_MEMORY, but for now we'll 
 * keep the check and allow the caller to handle success/failure via return codes.
 */
#define SAFE_CALLOC(ptr, count, type) do { \
    (ptr) = (type*)calloc((count), sizeof(type)); \
    if ((ptr) == NULL && (count) > 0) return ERR_MEMORY; \
} while(0)

#define SAFE_MALLOC(ptr, count, type) do { \
    (ptr) = (type*)malloc((count) * sizeof(type)); \
    if ((ptr) == NULL && (count) > 0) return ERR_MEMORY; \
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

void init_random_seed_custom(ModelConfig *config, prng_state *rs);
int allocate_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
int load_param(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore, MCMCState *mstate);
int load_categories(ModelConfig *config, GenomicData *gdata);
void write_dgv(ModelConfig *config, GenomicData *gdata);
void output_model(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore);
void output_beta(ModelConfig *config, MCMCState *mstate, GenomicData *gdata);
void cleanup_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);

#endif // MCMC_CORE_IO_H
