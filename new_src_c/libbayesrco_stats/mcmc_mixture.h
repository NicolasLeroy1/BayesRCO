#ifndef MCMC_MIXTURE_H
#define MCMC_MIXTURE_H

#include "bayesrco_types.h"
#include "rng.h"

// Mixture model MCMC functions
void mcmc_mixture_init(ModelConfig *config, GenomicData *gdata, MCMCState *mstate);
void mcmc_mixture_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs);

#endif // MCMC_MIXTURE_H
