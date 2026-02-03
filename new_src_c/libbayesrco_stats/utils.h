#ifndef UTILS_H
#define UTILS_H

#include "bayesrco_types.h"
#include "rng.h"

// From mod_stats.f90
void permutate(prng_state *rs, int n, int *p);
void compute_dgv(GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
void compute_residuals(GenomicData *gdata, MCMCState *mstate);

// From mod_standardize.f90 - pure computation version (no file I/O)
void standardize_genotypes(GenomicData *gdata, bool compute_freqs);

#endif // UTILS_H
