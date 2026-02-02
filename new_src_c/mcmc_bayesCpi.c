/**
 * @file mcmc_bayesCpi.c
 * @brief BayesCpi MCMC kernel implementation.
 * 
 * Implements the BayesCpi model where the mixture proportion (pi) is estimated
 * from the data rather than fixed. This allows adaptive estimation of the
 * proportion of SNPs with non-zero effects.
 */

#include "mcmc_bayesCpi.h"
#include "mcmc_utils.h"
#include "mcmc_sampling.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

/**
 * Initialize BayesCpi-specific data structures.
 * 
 * Sets up the annotation tracking matrices and initializes all SNPs
 * to distribution 2 (first non-null distribution) for SNPs in each category.
 */
void mcmc_bayesCpi_init(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;
    
    /* Initialize effects per category to zero */
    for (int i = 0; i < nloci; i++) {
        for (int j = 0; j < ncat; j++) {
            gdata->effects_per_category[i][j] = 0.0;
        }
    }
    
    /* Initialize distribution assignment: 
     * SNPs in a category start in distribution 2 (1-based indexing for output compatibility) 
     * Distribution 1 = null (zero effect), Distribution 2+ = non-null
     */
    for (int j = 0; j < ncat; j++) {
        for (int i = 0; i < nloci; i++) {
            if (gdata->categories[i][j] == 1) {
                gdata->distribution_per_category[i][j] = 2;
            }
        }
    }
}

/**
 * Run one iteration of the BayesCpi MCMC kernel.
 * 
 * For each category and each SNP in that category:
 * 1. Add back the current effect to the adjusted phenotype
 * 2. Compute log selection probabilities for each distribution
 * 3. Sample a new distribution assignment
 * 4. Sample a new effect (0 for null distribution)
 * 5. Subtract the new effect from the adjusted phenotype
 * 
 * After processing all SNPs, compute summary statistics.
 */
void mcmc_bayesCpi_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;
    int ndist = config->num_distributions;
    int nt = gdata->num_phenotyped_individuals;
    
    /* Temporary arrays for selection probabilities */
    double *log_probs = mstate->log_likelihoods;
    double *probs = mstate->selection_probs;
    
    /* Process each category */
    for (int j = 0; j < ncat; j++) {
        /* Compute log mixture proportions for this category */
        double *log_mix_probs = mstate->log_p[0];  /* Note: this should index into the j-th column */
        for (int i = 0; i < ndist; i++) {
            /* Access p[i][j] and compute log */
            log_mix_probs[i] = log(mstate->p[i][j]);
        }

        /* Process each locus in permuted order */
        for (int k = 0; k < nloci; k++) {
            int snploc = gdata->permvec[k];
            
            /* Only process if this locus belongs to this category */
            if (gdata->categories[snploc][j] != 1) continue;
            
            double ssq = gdata->snp_correlations[snploc];
            double current_effect = gdata->effects_per_category[snploc][j];
            int current_dist = gdata->distribution_per_category[snploc][j];
            
            /* If previously included (non-null), add effect back to adjusted phenotype */
            if (current_dist > 1) {
                add_snp_contribution(mstate->adjusted_phenotypes, gdata->genotypes,
                                    snploc, current_effect, nt, nloci);
            }
            
            /* Compute right-hand side: X'y_adj */
            double rhs = compute_rhs(gdata->genotypes, snploc, 
                                    mstate->adjusted_phenotypes, nt, nloci);
            
            /* Compute log selection probabilities */
            /* Note: need to build log_mix_probs array for column j */
            double log_mix_for_cat[DEFAULT_NUM_DISTRIBUTIONS];
            for (int d = 0; d < ndist; d++) {
                log_mix_for_cat[d] = log(mstate->p[d][j]);
            }
            
            compute_log_selection_probs(log_probs, rhs, ssq, 
                                       mstate->variance_residual,
                                       mstate->genomic_values,  /* dist_variances */
                                       log_mix_for_cat, ndist);
            
            /* Stabilize and normalize probabilities */
            stabilize_log_probs(probs, log_probs, ndist);
            
            /* Sample distribution index (0-based) */
            int dist_idx = sample_distribution_index(probs, ndist, rs);
            
            /* Store distribution assignment (1-based for output compatibility) */
            gdata->distribution_per_category[snploc][j] = dist_idx + 1;
            gdata->current_distribution[snploc] = dist_idx + 1;
            
            /* Sample effect */
            double new_effect;
            if (dist_idx == 0) {
                new_effect = 0.0;
            } else {
                new_effect = sample_snp_effect(dist_idx, rhs, ssq,
                                              mstate->variance_residual,
                                              mstate->genomic_values, rs);
                /* Subtract effect from adjusted phenotype */
                subtract_snp_contribution(mstate->adjusted_phenotypes, gdata->genotypes,
                                         snploc, new_effect, nt, nloci);
                mstate->included++;
            }
            
            gdata->effects_per_category[snploc][j] = new_effect;
            
            /* Early exit if marker set size reached */
            if (config->marker_set_size > 0 && mstate->rep > config->marker_replicates) {
                if (mstate->included >= config->marker_set_size) break;
            }
        }
    }
    
    /* Sum effects across categories: g[i] = sum_j(gannot[i][j]) */
    for (int i = 0; i < nloci; i++) {
        double sum = 0.0;
        for (int j = 0; j < ncat; j++) {
            sum += gdata->effects_per_category[i][j];
        }
        mstate->snp_effects[i] = sum;
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
    
    /* Count included loci (any category has non-null assignment) */
    for (int i = 0; i < nloci; i++) {
        gdata->included_loci[i] = 0.0;
    }
    for (int j = 0; j < ncat; j++) {
        for (int i = 0; i < nloci; i++) {
            if (gdata->distribution_per_category[i][j] > 1) {
                gdata->included_loci[i] = 1.0;
            }
        }
    }
    
    double included_sum = 0.0;
    for (int i = 0; i < nloci; i++) {
        included_sum += gdata->included_loci[i];
    }
    mstate->included = (int)included_sum;
}
