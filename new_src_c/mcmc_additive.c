/**
 * @file mcmc_additive.c
 * @brief Additive model MCMC kernel implementation.
 * 
 * Implements an additive model where each SNP's effect is the sum of
 * its effects across all annotation categories it belongs to.
 */

#include "mcmc_additive.h"
#include "mcmc_utils.h"
#include "mcmc_sampling.h"
#include <math.h>
#include <stdbool.h>

/**
 * Run one iteration of the additive model MCMC kernel.
 * 
 * For each SNP and each category it belongs to:
 * 1. Sample a distribution assignment
 * 2. Sample an effect for that category
 * 
 * The total SNP effect is the sum across all categories.
 */
void mcmc_additive_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;
    int ndist = config->num_distributions;
    int nt = gdata->num_phenotyped_individuals;
    
    /* Process each locus in permuted order */
    for (int k = 0; k < nloci; k++) {
        int snploc = gdata->permvec[k];
        double current_effect = mstate->snp_effects[snploc];
        
        double ssq = gdata->snp_correlations[snploc];
        double ssq_over_vare = ssq / mstate->variance_residual;
        
        /* If SNP has non-zero effect, add it back */
        if (gdata->current_distribution[snploc] > 1) {
            add_snp_contribution(mstate->adjusted_phenotypes, gdata->genotypes,
                                snploc, current_effect, nt, nloci);
        }
        
        /* Compute RHS: X'y_adj */
        double rhs = compute_rhs(gdata->genotypes, snploc, 
                                mstate->adjusted_phenotypes, nt, nloci);
        
        /* Process each category in permuted order */
        for (int kcat = 0; kcat < ncat; kcat++) {
            int j = gdata->permannot[kcat];
            
            /* Update log mixture probabilities */
            for (int i = 0; i < ndist; i++) {
                mstate->log_p[i][j] = log(mstate->p[i][j]);
            }
            
            /* Only process if this SNP belongs to this category */
            if (gdata->categories[snploc][j] != 1) continue;
            
            /* Compute log selection probabilities */
            double *log_probs = mstate->log_likelihoods;
            log_probs[0] = mstate->log_p[0][j];
            for (int kk = 1; kk < ndist; kk++) {
                double logdetV = log(mstate->genomic_values[kk] * ssq_over_vare + 1.0);
                double uhat = rhs / (ssq + mstate->residual_variance_over_distribution_variances[kk]);
                log_probs[kk] = -0.5 * (logdetV - (rhs * uhat / mstate->variance_residual)) + mstate->log_p[kk][j];
            }
            
            /* Stabilize and normalize */
            stabilize_log_probs(mstate->selection_probs, log_probs, ndist);
            
            /* Sample distribution */
            int dist_idx = sample_distribution_index(mstate->selection_probs, ndist, rs);
            
            /* Store assignment (1-based for compatibility) */
            gdata->distribution_per_category[snploc][j] = dist_idx + 1;
            
            /* Sample effect */
            double cat_effect;
            if (dist_idx == 0) {
                cat_effect = 0.0;
            } else {
                cat_effect = sample_snp_effect(dist_idx, rhs, ssq,
                                              mstate->variance_residual,
                                              mstate->genomic_values, rs);
                /* Subtract effect from adjusted phenotype */
                subtract_snp_contribution(mstate->adjusted_phenotypes, gdata->genotypes,
                                         snploc, cat_effect, nt, nloci);
            }
            gdata->effects_per_category[snploc][j] = cat_effect;
        }
    }
    
    /* Sum effects across categories and determine max distribution */
    for (int i = 0; i < nloci; i++) {
        double sum_effect = 0.0;
        int max_track = 0;
        for (int j = 0; j < ncat; j++) {
            sum_effect += gdata->effects_per_category[i][j];
            if (gdata->distribution_per_category[i][j] > max_track) {
                max_track = gdata->distribution_per_category[i][j];
            }
        }
        mstate->snp_effects[i] = sum_effect;
        gdata->current_distribution[i] = max_track;
    }
    
    /* Compute per-distribution statistics */
    for (int j = 0; j < ncat; j++) {
        for (int d = 0; d < ndist; d++) {
            int count = 0;
            double sum_sq = 0.0;
            for (int k = 0; k < nloci; k++) {
                if (gdata->distribution_per_category[k][j] == d + 1) {
                    count++;
                    double eff = gdata->effects_per_category[k][j];
                    sum_sq += eff * eff;
                }
            }
            mstate->snps_per_distribution[d][j] = count;
            mstate->variance_per_distribution[d][j] = sum_sq;
        }
    }
    
    /* Count included loci (any distribution > 1 means included) */
    double included_sum = 0.0;
    for (int i = 0; i < nloci; i++) {
        if (gdata->current_distribution[i] > 1) {
            gdata->included_loci[i] = 1.0;
            included_sum += 1.0;
        } else {
            gdata->included_loci[i] = 0.0;
        }
    }
    mstate->included = (int)included_sum;
}
