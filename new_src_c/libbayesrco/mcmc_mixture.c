/**
 * @file mcmc_mixture.c
 * @brief Mixture model MCMC kernel implementation.
 * 
 * Implements the BayesRC mixture model where SNPs can belong to multiple
 * annotation categories, and the category assignment is sampled during MCMC.
 */

#include "mcmc_mixture.h"
#include "mcmc_utils.h"
#include "mcmc_sampling.h"
#include <math.h>
#include <stdbool.h>

/**
 * Initialize mixture model-specific data structures.
 */
void mcmc_mixture_init(ModelConfig *config, GenomicData *gdata, MCMCState *mstate) {
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;
    
    /* Initialize distribution tracking to 2 (first non-null) */
    for (int k = 0; k < nloci; k++) {
        gdata->current_distribution[k] = 2;
    }
    
    /* For SNPs with single annotation, set current category */
    for (int k = 0; k < nloci; k++) {
        if (gdata->annotations_per_locus[k] == 1) {
            for (int j = 0; j < ncat; j++) {
                if (gdata->categories[k][j] == 1) {
                    gdata->current_category[k] = j;
                }
            }
        }
    }
}

/**
 * Run one iteration of the mixture model MCMC kernel.
 * 
 * Two-stage sampling:
 * 1. Sample annotation category for multi-annotated SNPs
 * 2. Sample distribution and effect for each SNP
 */
void mcmc_mixture_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;
    int ndist = config->num_distributions;
    int nt = gdata->num_phenotyped_individuals;
    
    /* Reset distribution per category tracking */
    for (int k = 0; k < nloci; k++) {
        for (int j = 0; j < ncat; j++) {
            gdata->distribution_per_category[k][j] = 0;
        }
    }

    /* =====================================================================
     * Stage 1: Sample annotation category for multi-annotated SNPs
     * ===================================================================== */
    for (int k = 0; k < nloci; k++) {
        int snploc = gdata->permvec[k];
        
        if (gdata->annotations_per_locus[snploc] > 1) {
            double ssq = gdata->snp_correlations[snploc];
            double current_effect = mstate->snp_effects[snploc];
            int current_dist = gdata->current_distribution[snploc];
            
            /* Prepare Dirichlet posterior for annotation sampling */
            for (int j = 0; j < ncat; j++) {
                mstate->category_dirichlet_scratch[j] = (double)gdata->categories[snploc][j];
            }
            if (mstate->rep != 1) {
                /* Add count for current category */
                mstate->category_dirichlet_scratch[gdata->current_category[snploc]] += 1.0;
            }
            
            /* Sample annotation probabilities from Dirichlet */
            rng_dirichlet2(rs, ncat, mstate->category_dirichlet_scratch, mstate->category_probabilities);
            
            /* Compute adjusted phenotype with current effect added back */
            double *ytemp = mstate->ytemp;
            if (current_dist > 1) {
                for (int row = 0; row < nt; row++) {
                    ytemp[row] = mstate->adjusted_phenotypes[row] + 
                                 gdata->genotypes[snploc * nt + row] * current_effect;
                }
            } else {
                for (int row = 0; row < nt; row++) {
                    ytemp[row] = mstate->adjusted_phenotypes[row];
                }
            }
            
            /* Compute RHS */
            double rhs = dot_product_col(gdata->genotypes, snploc, nt, nloci, ytemp);
            
            /* Find maximum log probability for stability */
            double ssq_over_vare = ssq / mstate->variance_residual;
            double maxs = 0.0;
            for (int i = 1; i < ndist; i++) {
                double uhat = rhs / (ssq + mstate->residual_variance_over_distribution_variances[i]);
                double temp = 0.5 * uhat * rhs / mstate->variance_residual;
                if (temp > maxs) maxs = temp;
            }
            
            /* Compute annotation selection scores */
            for (int j = 0; j < ncat; j++) {
                if (gdata->categories[snploc][j] == 1) {
                    double ss_val = mstate->p[0][j] * exp(-maxs);
                    for (int kk = 1; kk < ndist; kk++) {
                        double detV = mstate->genomic_values[kk] * ssq_over_vare + 1.0;
                        double uhat = rhs / (ssq + mstate->residual_variance_over_distribution_variances[kk]);
                        ss_val += mstate->p[kk][j] * pow(detV, -0.5) * 
                                  exp(0.5 * uhat * rhs / mstate->variance_residual - maxs);
                    }
                    mstate->ss[j] = log(mstate->category_probabilities[j]) + log(ss_val);
                }
            }
            
            /* Stabilize and normalize */
            for (int kk = 0; kk < ncat; kk++) {
                if (gdata->categories[snploc][kk] == 1) {
                    double skk = mstate->ss[kk];
                    double sk = 0.0;
                    bool overflow = false;
                    for (int l = 0; l < ncat; l++) {
                        if (gdata->categories[snploc][l] == 1 && l != kk) {
                            double clike = mstate->ss[l] - skk;
                            if (clike > LOG_UPPER_LIMIT) { overflow = true; break; }
                            if (clike < -LOG_UPPER_LIMIT) continue;
                            sk += exp(clike);
                        }
                    }
                    mstate->sstemp[kk] = overflow ? 0.0 : 1.0 / (1.0 + sk);
                } else {
                    mstate->sstemp[kk] = 0.0;
                }
            }
            
            /* Sample annotation */
            gdata->current_category[snploc] = sample_discrete(mstate->sstemp, ncat, rs);
        }
    }
    
    /* =====================================================================
     * Stage 2: Sample distribution and effect for each SNP
     * ===================================================================== */
    for (int k = 0; k < nloci; k++) {
        int snploc = gdata->permvec[k];
        int j = gdata->current_category[snploc];
        
        double ssq = gdata->snp_correlations[snploc];
        double current_effect = mstate->snp_effects[snploc];
        int current_dist = gdata->current_distribution[snploc];
        
        /* Add back current effect to adjusted phenotype */
        if (current_dist > 1) {
            add_col_scalar(mstate->adjusted_phenotypes, gdata->genotypes,
                           snploc, nt, nloci, current_effect);
        }
        
        /* Compute RHS */
        double rhs = dot_product_col(gdata->genotypes, snploc, nt, nloci, 
                                    mstate->adjusted_phenotypes);
        
        /* Compute log selection probabilities */
        double *log_probs = mstate->log_likelihoods;
        log_probs[0] = mstate->log_p[0][j];
        for (int kk = 1; kk < ndist; kk++) {
            double ssq_over_vare = ssq / mstate->variance_residual;
            double logdetV = log(mstate->genomic_values[kk] * ssq_over_vare + 1.0);
            double uhat = rhs / (ssq + mstate->residual_variance_over_distribution_variances[kk]);
            log_probs[kk] = -0.5 * (logdetV - (rhs * uhat / mstate->variance_residual)) + mstate->log_p[kk][j];
        }
        
        /* Stabilize and normalize */
        stabilize_log_probs(mstate->selection_probs, log_probs, ndist);
        
        /* Sample distribution */
        int dist_idx = sample_discrete(mstate->selection_probs, ndist, rs);
        
        /* Store assignments (1-based for compatibility) */
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
            add_col_scalar(mstate->adjusted_phenotypes, gdata->genotypes,
                           snploc, nt, nloci, -new_effect);
            mstate->included++;
        }
        mstate->snp_effects[snploc] = new_effect;
    }
    
    /* =====================================================================
     * Compute per-distribution statistics
     * ===================================================================== */
    for (int j = 0; j < ncat; j++) {
        for (int i = 0; i < ndist; i++) {
            int cnt = 0;
            double sum_g2 = 0.0;
            for (int k = 0; k < nloci; k++) {
                if (gdata->distribution_per_category[k][j] == i + 1) {
                    cnt++;
                    sum_g2 += mstate->snp_effects[k] * mstate->snp_effects[k];
                }
            }
            mstate->snps_per_distribution[i][j] = cnt;
            mstate->variance_per_distribution[i][j] = sum_g2;
        }
    }
    
    /* Count included SNPs (distribution > 1) */
    int count_zero = 0;
    for (int j = 0; j < ncat; j++) {
        count_zero += mstate->snps_per_distribution[0][j];
    }
    mstate->included = nloci - count_zero;
}
