#ifndef UTILS_H
#define UTILS_H

#include "bayesRCO.h"

// From mod_stats.f90
void permutate(prng_state *rs, int n, int *p);
void compute_dgv(GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
void compute_residuals(GenomicData *gdata, MCMCState *mstate);

// From mod_standardize.f90
void xcenter(ModelConfig *config, GenomicData *gdata);

#endif // UTILS_H
