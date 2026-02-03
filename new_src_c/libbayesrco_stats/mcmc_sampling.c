/**
 * @file mcmc_sampling.c
 * @brief Implementation of shared MCMC sampling functions.
 */

#include "mcmc_sampling.h"
#include <math.h>

#ifndef LOG_UPPER_LIMIT
#define LOG_UPPER_LIMIT 700.0
#endif

/* =========================================================================
 * Matrix and Vector Operations
 * ========================================================================= */

double dot_product_col(double *X, int col_idx, int n_rows, int n_cols, double *vec) {
    double sum = 0.0;
    int offset = col_idx * n_rows;
    for (int i = 0; i < n_rows; i++) {
        sum += X[offset + i] * vec[i];
    }
    return sum;
}

void add_col_scalar(double *vec, double *X, int col_idx, int n_rows, int n_cols, double scalar) {
    int offset = col_idx * n_rows;
    for (int i = 0; i < n_rows; i++) {
        vec[i] += X[offset + i] * scalar;
    }
}

/* =========================================================================
 * MCMC Specific Sampling
 * ========================================================================= */

void compute_log_selection_probs(
    double *log_probs,
    double rhs,
    double ssq,
    double vare,
    double *dist_variances,
    double *log_mix_probs,
    int num_dist
) {
    double ssq_over_vare = ssq / vare;
    
    /* Distribution 0 is the null (zero effect) distribution */
    log_probs[0] = log_mix_probs[0];
    
    /* For other distributions, compute posterior likelihood */
    for (int i = 1; i < num_dist; i++) {
        double vare_over_gp = vare / dist_variances[i];
        double log_det_V = log(dist_variances[i] * ssq_over_vare + 1.0);
        double uhat = rhs / (ssq + vare_over_gp);
        
        /* Log posterior probability:
         * -0.5 * (log(det(V)) - (rhs * uhat / vare)) + log(mix_prob)
         */
        log_probs[i] = -0.5 * (log_det_V - (rhs * uhat / vare)) + log_mix_probs[i];
    }
}

void stabilize_log_probs(double *probs, double *log_probs, int n) {
    for (int k = 0; k < n; k++) {
        double skk = log_probs[k];
        double sum_exp = 0.0;
        int overflow = 0;
        
        for (int j = 0; j < n; j++) {
            if (j == k) continue;
            double diff = log_probs[j] - skk;
            
            if (diff > LOG_UPPER_LIMIT) {
                overflow = 1;
                break;
            }
            if (diff < -LOG_UPPER_LIMIT) continue;
            
            sum_exp += exp(diff);
        }
        
        if (overflow) {
            probs[k] = 0.0;
        } else {
            probs[k] = 1.0 / (1.0 + sum_exp);
        }
    }
}

int sample_discrete(double *probs, int n, prng_state *rs) {
    double r = rng_uniform(rs, 0.0, 1.0);
    double cumsum = 0.0;
    
    for (int i = 0; i < n; i++) {
        cumsum += probs[i];
        if (r < cumsum) {
            return i;
        }
    }
    return n - 1;
}

double sample_snp_effect(
    int dist_idx,
    double rhs,
    double ssq,
    double vare,
    double *dist_variances,
    prng_state *rs
) {
    if (dist_idx == 0) return 0.0;
    
    double v1 = ssq + vare / dist_variances[dist_idx];
    double mean = rhs / v1;
    double sd = sqrt(vare / v1);
    
    return rng_normal(rs, mean, sd);
}
