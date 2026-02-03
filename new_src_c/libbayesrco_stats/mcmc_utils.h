#ifndef MCMC_UTILS_H
#define MCMC_UTILS_H

#include "bayesrco_types.h"
#include "rng.h"
#include "mcmc_sampling.h"

// Common MCMC utility functions
void mcmc_save_samples_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
void mcmc_init_common(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore);
void mcmc_start_values_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs);
void mcmc_iteration_pre_common(ModelConfig *config, MCMCState *mstate, prng_state *rs);
void mcmc_update_hypers_common(int nc, ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs);

#endif // MCMC_UTILS_H
