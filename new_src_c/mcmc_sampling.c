/**
 * @file mcmc_sampling.c
 * @brief Implementation of shared MCMC sampling functions.
 * 
 * Contains common sampling routines used across mixture, additive, and BayesCpi kernels.
 */

#include "mcmc_sampling.h"
#include <math.h>

/* Log upper limit for numerical stability */
#ifndef LOG_UPPER_LIMIT
#define LOG_UPPER_LIMIT 700.0
#endif

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
         * -0.5 * log(det(V)) + 0.5 * (rhs * uhat / vare) + log(mix_prob)
         */
        log_probs[i] = -0.5 * (log_det_V - (rhs * uhat / vare)) + log_mix_probs[i];
    }
}

void stabilize_log_probs(double *probs, double *log_probs, int n) {
    /*
     * Use log-sum-exp trick to convert log probabilities to normalized probabilities.
     * For each k, compute: P(k) = 1 / (1 + sum_{j!=k} exp(log_probs[j] - log_probs[k]))
     * 
     * This formulation is numerically stable because we're computing relative probabilities.
     */
    for (int k = 0; k < n; k++) {
        double skk = log_probs[k];
        double sum_exp = 0.0;
        int overflow = 0;
        
        for (int j = 0; j < n; j++) {
            if (j == k) continue;
            
            double diff = log_probs[j] - skk;
            
            if (diff > LOG_UPPER_LIMIT) {
                /* This distribution dominates - probability is essentially 0 for k */
                overflow = 1;
                break;
            }
            
            if (diff < -LOG_UPPER_LIMIT) {
                /* This distribution contributes negligibly */
                continue;
            }
            
            sum_exp += exp(diff);
        }
        
        if (overflow) {
            probs[k] = 0.0;
        } else {
            probs[k] = 1.0 / (1.0 + sum_exp);
        }
    }
}

int sample_distribution_index(double *probs, int n, prng_state *rs) {
    double r = rng_uniform(rs, 0.0, 1.0);
    double cumsum = 0.0;
    
    for (int i = 0; i < n; i++) {
        cumsum += probs[i];
        if (r < cumsum) {
            return i;
        }
    }
    
    /* Fallback to last distribution (handles floating point rounding) */
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
    /* Null distribution (index 0) has zero effect */
    if (dist_idx == 0) {
        return 0.0;
    }
    
    /* Sample from posterior normal distribution:
     * mean = rhs / (ssq + vare/gp)
     * variance = vare / (ssq + vare/gp)
     */
    double v1 = ssq + vare / dist_variances[dist_idx];
    double mean = rhs / v1;
    double sd = sqrt(vare / v1);
    
    return rng_normal(rs, mean, sd);
}

void subtract_snp_contribution(
    double *y_adj,
    double *X,
    int snp_idx,
    double effect,
    int n_ind,
    int n_loci
) {
    /* Column-major access: X[snp_idx * n_ind + i] for individual i */
    for (int i = 0; i < n_ind; i++) {
        y_adj[i] -= X[snp_idx * n_ind + i] * effect;
    }
}

void add_snp_contribution(
    double *y_adj,
    double *X,
    int snp_idx,
    double effect,
    int n_ind,
    int n_loci
) {
    /* Column-major access: X[snp_idx * n_ind + i] for individual i */
    for (int i = 0; i < n_ind; i++) {
        y_adj[i] += X[snp_idx * n_ind + i] * effect;
    }
}

double compute_rhs(
    double *X,
    int snp_idx,
    double *y_adj,
    int n_ind,
    int n_loci
) {
    double sum = 0.0;
    
    /* Column-major access: X[snp_idx * n_ind + i] for individual i */
    for (int i = 0; i < n_ind; i++) {
        sum += X[snp_idx * n_ind + i] * y_adj[i];
    }
    
    return sum;
}
