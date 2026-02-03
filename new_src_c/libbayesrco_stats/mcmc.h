#ifndef MCMC_H
#define MCMC_H

#include "bayesrco_types.h"
#include "rng.h"
#include "rng.h"

int run_mcmc(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs);

#endif // MCMC_H
