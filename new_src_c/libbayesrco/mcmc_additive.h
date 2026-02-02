#ifndef MCMC_ADDITIVE_H
#define MCMC_ADDITIVE_H

#include "bayesRCO.h"

// Additive model MCMC functions
void mcmc_additive_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs);

#endif // MCMC_ADDITIVE_H
