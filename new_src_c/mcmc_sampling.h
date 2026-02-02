#ifndef MCMC_SAMPLING_H
#define MCMC_SAMPLING_H

/**
 * @file mcmc_sampling.h
 * @brief Shared MCMC sampling functions for distribution selection and effect sampling.
 * 
 * This module contains common functions extracted from the mixture, additive, and BayesCpi
 * MCMC kernels to reduce code duplication and improve maintainability.
 */

#include "bayesRCO.h"
#include "rng.h"

/**
 * Compute log selection probabilities for each distribution.
 * 
 * Calculates the log-likelihood of each distribution for a given SNP based on:
 * - The right-hand side of the normal equation (rhs = X'y_adj)
 * - The sum of squares for the SNP (ssq = X'X)
 * - The current variance estimates
 * 
 * @param log_probs       Output array of log selection probabilities (size: num_dist)
 * @param rhs             Right-hand side of normal equation (X' * y_adj)
 * @param ssq             Sum of squares for current SNP (X' * X)
 * @param vare            Residual variance
 * @param dist_variances  Variance for each distribution (size: num_dist)
 * @param log_mix_probs   Log mixture proportions for current category (size: num_dist)
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
 * Converts log probabilities to normalized probabilities while avoiding
 * numerical overflow/underflow issues.
 * 
 * @param probs      Output array of normalized probabilities (size: n)
 * @param log_probs  Input array of log probabilities (size: n)
 * @param n          Number of probabilities
 */
void stabilize_log_probs(double *probs, double *log_probs, int n);

/**
 * Sample a distribution index from normalized probabilities.
 * 
 * @param probs   Normalized probability array (size: n)
 * @param n       Number of distributions
 * @param rs      Random number generator state
 * @return        Sampled distribution index (0-based)
 */
int sample_distribution_index(double *probs, int n, prng_state *rs);

/**
 * Sample SNP effect given the selected distribution.
 * 
 * For the null distribution (dist_idx == 0), returns 0.
 * For other distributions, samples from the posterior normal distribution.
 * 
 * @param dist_idx        Selected distribution index (0-based)
 * @param rhs             Right-hand side of normal equation
 * @param ssq             Sum of squares for current SNP
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

/**
 * Update the adjusted phenotype vector after sampling a new SNP effect.
 * 
 * Subtracts the new effect contribution from y_adj:
 *   y_adj = y_adj - X[:, snp_idx] * effect
 * 
 * @param y_adj       Adjusted phenotype vector (modified in-place)
 * @param X           Genotype matrix (column-major: nloci x nt)
 * @param snp_idx     Index of current SNP (column index)
 * @param effect      Effect size to subtract
 * @param n_ind       Number of individuals
 * @param n_loci      Number of loci
 */
void subtract_snp_contribution(
    double *y_adj,
    double *X,
    int snp_idx,
    double effect,
    int n_ind,
    int n_loci
);

/**
 * Add a SNP contribution back to the adjusted phenotype vector.
 * 
 * Adds the effect contribution back to y_adj before resampling:
 *   y_adj = y_adj + X[:, snp_idx] * effect
 * 
 * @param y_adj       Adjusted phenotype vector (modified in-place)
 * @param X           Genotype matrix (column-major: nloci x nt)
 * @param snp_idx     Index of current SNP (column index)
 * @param effect      Effect size to add back
 * @param n_ind       Number of individuals
 * @param n_loci      Number of loci
 */
void add_snp_contribution(
    double *y_adj,
    double *X,
    int snp_idx,
    double effect,
    int n_ind,
    int n_loci
);

/**
 * Compute the right-hand side of the effect sampling equation.
 * 
 * Calculates: rhs = X[:, snp_idx]' * y_adj
 * 
 * @param X           Genotype matrix (column-major: nloci x nt)
 * @param snp_idx     Index of current SNP (column index)
 * @param y_adj       Adjusted phenotype vector
 * @param n_ind       Number of individuals
 * @param n_loci      Number of loci
 * @return            Right-hand side value
 */
double compute_rhs(
    double *X,
    int snp_idx,
    double *y_adj,
    int n_ind,
    int n_loci
);

#endif // MCMC_SAMPLING_H
