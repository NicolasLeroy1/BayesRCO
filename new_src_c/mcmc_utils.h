#ifndef MCMC_UTILS_H
#define MCMC_UTILS_H

#include "bayesRCO.h"

// Helper functions for matrix operations
double dot_product_col(double *X, int col_idx, int n_rows, int n_cols, double *vec);
void add_col_scalar(double *vec, double *X, int col_idx, int n_rows, int n_cols, double scalar);

// Helper for random discrete sampling
int sample_discrete(double *probs, int n, prng_state *rs);

// Common MCMC utility functions
void mcmc_save_samples_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
void mcmc_init_common(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore);
void mcmc_start_values_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs);
void mcmc_iteration_pre_common(ModelConfig *config, MCMCState *mstate, prng_state *rs);
void mcmc_update_hypers_common(int nc, ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs);

#endif // MCMC_UTILS_H
