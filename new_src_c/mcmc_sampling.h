#ifndef MCMC_SAMPLING_H
#define MCMC_SAMPLING_H

/**
 * @file mcmc_sampling.h
 * @brief Shared MCMC sampling functions for distribution selection and effect sampling.
 * 
 * This module contains core math and sampling routines used across all MCMC kernels.
 */

#include "bayesRCO.h"
#include "rng.h"

/* =========================================================================
 * Matrix and Vector Operations
 * ========================================================================= */

/**
 * Compute dot product of a column with a vector.
 * Calculates: sum = X[:, col_idx]' * vec
 * 
 * @param X           Genotype matrix (column-major: n_rows x n_cols)
 * @param col_idx     Index of current column
 * @param n_rows      Number of individuals (rows)
 * @param n_cols      Number of loci (columns)
 * @param vec         Vector to dot with (e.g., adjusted phenotypes)
 * @return            Dot product result
 */
double dot_product_col(double *X, int col_idx, int n_rows, int n_cols, double *vec);

/**
 * Add column * scalar to vector: vec += X[:, col_idx] * scalar
 * 
 * @param vec         Target vector (modified in-place)
 * @param X           Genotype matrix (column-major)
 * @param col_idx     Index of current column
 * @param n_rows      Number of individuals
 * @param n_cols      Number of loci
 * @param scalar      Multiplier
 */
void add_col_scalar(double *vec, double *X, int col_idx, int n_rows, int n_cols, double scalar);

/* =========================================================================
 * MCMC Specific Sampling
 * ========================================================================= */

/**
 * Compute log selection probabilities for each distribution.
 * 
 * @param log_probs       Output array of log selection probabilities (size: num_dist)
 * @param rhs             Right-hand side (X' * y_adj)
 * @param ssq             Sum of squares (X' * X)
 * @param vare            Residual variance
 * @param dist_variances  Variance for each distribution
 * @param log_mix_probs   Log mixture proportions
 * @param num_dist        Number of distributions
 */
void compute_log_selection_probs(
    double *log_probs,
    double rhs,
    double ssq,
    double vare,
    double *dist_variances,
    double *log_mix_probs,
    int num_dist
);

/**
 * Stabilize and normalize log probabilities using the log-sum-exp trick.
 * 
 * @param probs      Output array of normalized probabilities (size: n)
 * @param log_probs  Input array of log probabilities (size: n)
 * @param n          Number of probabilities
 */
void stabilize_log_probs(double *probs, double *log_probs, int n);

/**
 * Sample a distribution index from (potentially unnormalized) probabilities.
 * 
 * @param probs   Probability array (size: n)
 * @param n       Number of items
 * @param rs      Random number generator state
 * @return        Sampled index (0-based)
 */
int sample_discrete(double *probs, int n, prng_state *rs);

/**
 * Sample SNP effect given the selected distribution.
 * 
 * @param dist_idx        Selected distribution index (0-based)
 * @param rhs             Right-hand side (X'y_adj)
 * @param ssq             Sum of squares (X'X)
 * @param vare            Residual variance
 * @param dist_variances  Variance for each distribution
 * @param rs              Random number generator state
 * @return                Sampled SNP effect
 */
double sample_snp_effect(
    int dist_idx,
    double rhs,
    double ssq,
    double vare,
    double *dist_variances,
    prng_state *rs
);

#endif // MCMC_SAMPLING_H
