#ifndef MCMC_H
#define MCMC_H

#include "bayesRCO.h"
#include "rng.h"

void run_mcmc(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs);

#endif // MCMC_H
