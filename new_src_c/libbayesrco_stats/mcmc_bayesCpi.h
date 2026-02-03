#ifndef MCMC_BAYESCPI_H
#define MCMC_BAYESCPI_H

#include "bayesrco_types.h"
#include "rng.h"

// BayesCpi model MCMC functions
void mcmc_bayesCpi_init(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
void mcmc_bayesCpi_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs);

#endif // MCMC_BAYESCPI_H
