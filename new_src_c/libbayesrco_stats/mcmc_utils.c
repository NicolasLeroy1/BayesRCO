/**
 * @file mcmc_utils.c
 * @brief Common MCMC utility functions.
 * 
 * Contains shared functions for MCMC initialization, iteration updates,
 * hyperparameter sampling, and sample accumulation.
 */

#include "mcmc_utils.h"
#include "utils.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
 * MCMC Common Functions
 * ========================================================================= */

/**
 * Save samples for posterior summary.
 * Accumulates SNP effects, variances, and hyperparameters.
 * 
 * Note: File I/O (writing to hyp file) is handled by the calling layer.
 */
void mcmc_save_samples_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    int nloci = gdata->num_loci;
    int ndist = config->num_distributions;
    int ncat = config->num_categories;
    
    mstate->counter++;
    
    /* Accumulate SNP effects */
    for (int i = 0; i < nloci; i++) {
        mstore->sum_snp_effects[i] += mstate->snp_effects[i];
        mstore->varistore[i] += mstate->snp_effects[i] * mstate->snp_effects[i];
    }
    
    /* Accumulate distribution statistics */
    for (int j = 0; j < ncat; j++) {
        for (int i = 0; i < ndist; i++) {
            mstore->sum_mixture_proportions[IDX2(i, j, ncat)] += mstate->p[IDX2(i, j, ncat)];
            mstore->sum_snps_per_distribution[IDX2(i, j, ncat)] += (double)mstate->snps_per_distribution[IDX2(i, j, ncat)];
            mstore->sum_variance_per_distribution[IDX2(i, j, ncat)] += mstate->variance_per_distribution[IDX2(i, j, ncat)];
        }
    }
    
    /* Accumulate hyperparameters: mu, nsnp, vara, vare */
    mstore->mu_vare_store[0] += mstate->mu;
    mstore->mu_vare_store[1] += mstate->included;
    mstore->mu_vare_store[2] += mstate->variance_genetic;
    mstore->mu_vare_store[3] += mstate->variance_residual;
    
    /* Online variance calculation (Welford's algorithm) */
    if (mstate->counter > 1) {
        for (int i = 0; i < nloci; i++) {
            double term = ((double)mstate->counter * mstate->snp_effects[i] - mstore->sum_snp_effects[i]);
            mstore->varustore[i] += term * term / ((double)mstate->counter * ((double)mstate->counter - 1.0));
        }
    }
    
    /* Note: File I/O (fprintf to config->fp_hyp) removed. 
     * I/O is handled by the calling layer (libbayesrco_io). */
}

/**
 * Initialize common MCMC data structures.
 * Computes genotype sum of squares and annotation counts.
 */
void mcmc_init_common(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore) {
    int nloci = gdata->num_loci;
    int nt = gdata->num_phenotyped_individuals;
    int ncat = config->num_categories;
    
    /* Compute X'X for each SNP (genotype sum of squares) */
    for (int i = 0; i < nloci; i++) {
        double sum = 0.0;
        for (int r = 0; r < nt; r++) {
            double val = gdata->genotypes[i * nt + r];
            sum += val * val;
        }
        gdata->snp_correlations[i] = sum;
    }
    
    /* Count annotations per locus */
    for (int i = 0; i < nloci; i++) {
        int sum_c = 0;
        for (int j = 0; j < ncat; j++) {
            sum_c += gdata->categories[IDX2(i, j, ncat)];
        }
        gdata->annotations_per_locus[i] = sum_c;
    }
}

/**
 * Set starting values for MCMC state.
 */
void mcmc_start_values_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs) {
    int nloci = gdata->num_loci;
    int nind = gdata->num_individuals;
    int nt = gdata->num_phenotyped_individuals;
    int ndist = config->num_distributions;
    int ncat = config->num_categories;
    
    mstate->mu = 1.0;
    for (int i = 0; i < nt; i++) {
        mstate->adjusted_phenotypes[i] = 0.0;
    }
    
    /* Compute mean phenotype for training individuals */
    double sum_y = 0.0;
    for (int i = 0; i < nind; i++) {
        if (gdata->trains[i] == 0) {
            sum_y += gdata->phenotypes[i];
        }
    }
    mstate->yhat = sum_y / mstate->nnind;
    
    /* Compute phenotype variance */
    double sum_sq = 0.0;
    for (int i = 0; i < nind; i++) {
        if (gdata->trains[i] == 0) {
            double diff = gdata->phenotypes[i] - mstate->yhat;
            sum_sq += diff * diff;
        }
    }
    mstate->vary = sum_sq / (mstate->nnind - 1.0);
    
    /* Initialize distribution variances */
    for (int i = 0; i < ndist; i++) {
        mstate->genomic_values[i] = mstate->variance_scaling_factors[i] * mstate->variance_genetic;
    }
    
    /* Initialize mixture proportions */
    mstate->scale = 0.0;
    for (int j = 0; j < ncat; j++) {
        mstate->p[IDX2(0, j, ncat)] = 0.5;
        double sum_rest = 0.0;
        for (int i = 1; i < ndist; i++) {
            mstate->p[IDX2(i, j, ncat)] = 1.0 / mstate->variance_scaling_factors[i];
            sum_rest += mstate->p[IDX2(i, j, ncat)];
        }
        for (int i = 1; i < ndist; i++) {
            mstate->p[IDX2(i, j, ncat)] = 0.5 * mstate->p[IDX2(i, j, ncat)] / sum_rest;
        }
    }
    
    /* Initialize SNP effects */
    double g_val = sqrt(mstate->variance_genetic / (0.5 * (double)nloci));
    for (int i = 0; i < nloci; i++) {
        mstate->snp_effects[i] = g_val;
    }
    
    /* Initialize permutation vector */
    for (int k = 0; k < nloci; k++) {
        gdata->permvec[k] = k;
    }
    
    compute_residuals(gdata, mstate);
}

/**
 * Pre-iteration setup: sample residual variance and update intercept.
 */
void mcmc_iteration_pre_common(ModelConfig *config, MCMCState *mstate, prng_state *rs) {
    int nt = (int)mstate->nnind;
    
    mstate->included = 0;
    
    /* Sample residual variance (if not using VCE) */
    if (!config->VCE) {
        double dp = 0.0;
        for (int i = 0; i < nt; i++) {
            dp += mstate->adjusted_phenotypes[i] * mstate->adjusted_phenotypes[i];
        }
        mstate->variance_residual = dp / rng_chi_square(rs, mstate->nnind + 3.0);
    }
    
    /* Update intercept: add back old mu, sample new mu, subtract new mu */
    double sum_yadj = 0.0;
    for (int i = 0; i < nt; i++) {
        mstate->adjusted_phenotypes[i] += mstate->mu;
        sum_yadj += mstate->adjusted_phenotypes[i];
    }
    
    mstate->mu = rng_normal(rs, sum_yadj / mstate->nnind, sqrt(mstate->variance_residual / mstate->nnind));
    
    for (int i = 0; i < nt; i++) {
        mstate->adjusted_phenotypes[i] -= mstate->mu;
    }
    
    /* Precompute log distribution variances and vare/gp ratios */
    for (int i = 1; i < config->num_distributions; i++) {
        mstate->log_distribution_variances[i] = log(mstate->genomic_values[i]);
        mstate->residual_variance_over_distribution_variances[i] = mstate->variance_residual / mstate->genomic_values[i];
    }
}

/**
 * Update hyperparameters after processing all SNPs.
 */
void mcmc_update_hypers_common(int nc, ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs) {
    int nloci = gdata->num_loci;
    int nt = (int)mstate->nnind;
    int ndist = config->num_distributions;
    int ncat = config->num_categories;
    
    if (config->VCE) {
        /* Compute sum of squared effects */
        double sum_g2 = 0.0;
        for (int i = 0; i < nloci; i++) {
            sum_g2 += mstate->snp_effects[i] * mstate->snp_effects[i];
        }
        
        /* Sample genetic variance */
        mstate->scale = ((double)mstate->included * sum_g2 + config->vara_ap * config->dfvara) / 
                        (config->dfvara + (double)mstate->included);
        
        mstate->variance_genetic = rng_scaled_inverse_chi_square(rs, 
            (double)mstate->included + config->dfvara, mstate->scale);
        
        /* Update distribution variances */
        if (nc == 2) {
            /* BayesCpi: single non-null distribution */
            double denom = (mstate->included > 0) ? (double)mstate->included : 1.0;
            mstate->genomic_values[1] = mstate->variance_genetic / denom;
            if (mstate->included == 0) mstate->genomic_values[1] = 0.0;
        } else {
            for (int i = 0; i < ndist; i++) {
                mstate->genomic_values[i] = mstate->variance_scaling_factors[i] * mstate->variance_genetic;
            }
        }
        
        /* Sample residual variance */
        double dp = 0.0;
        for (int i = 0; i < nt; i++) {
            dp += mstate->adjusted_phenotypes[i] * mstate->adjusted_phenotypes[i];
        }
        
        double vare_num = dp + config->vare_ap * config->dfvare;
        mstate->variance_residual = vare_num / (mstate->nnind + config->dfvare);
        mstate->variance_residual = rng_scaled_inverse_chi_square(rs, 
            mstate->nnind + config->dfvare, mstate->variance_residual);
    }
    
    /* Sample mixture proportions from Dirichlet posterior */
    for (int j = 0; j < ncat; j++) {
        for (int i = 0; i < ndist; i++) {
            mstate->dirichlet_scratch[i] = (double)mstate->snps_per_distribution[IDX2(i, j, ncat)] + mstate->dirichlet_priors[i];
        }
        rng_dirichlet(rs, ndist, mstate->dirichlet_scratch, mstate->ytemp);
        
        for (int i = 0; i < ndist; i++) {
            mstate->p[IDX2(i, j, ncat)] = mstate->ytemp[i];
            mstate->log_p[IDX2(i, j, ncat)] = log(mstate->p[IDX2(i, j, ncat)]);
        }
    }
}
