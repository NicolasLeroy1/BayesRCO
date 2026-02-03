/**
 * @file bayesrco_stats.c
 * @brief Public API implementation for BayesRCO statistical library.
 * 
 * This file provides the main entry points that wrap the internal MCMC functions.
 */

#include "bayesrco_stats.h"
#include "mcmc.h"
#include "rng.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* =========================================================================
 * Helper Functions
 * ========================================================================= */

/**
 * Allocate MCMCResults arrays.
 */
int allocate_results(MCMCResults *results, int nloci, int ndist, int ncat, int nind) {
    if (!results) return ERR_INVALID_ARG;
    
    results->posterior_means = (double*)calloc(nloci, sizeof(double));
    if (!results->posterior_means) return ERR_MEMORY;
    
    results->posterior_vars = (double*)calloc(nloci, sizeof(double));
    if (!results->posterior_vars) { free(results->posterior_means); return ERR_MEMORY; }
    
    results->distribution_probs = (double*)calloc(nloci * ndist, sizeof(double));
    if (!results->distribution_probs) { 
        free(results->posterior_means); 
        free(results->posterior_vars); 
        return ERR_MEMORY; 
    }
    
    results->category_probs = (double*)calloc(nloci * ncat, sizeof(double));
    if (!results->category_probs) { 
        free(results->posterior_means); 
        free(results->posterior_vars); 
        free(results->distribution_probs);
        return ERR_MEMORY; 
    }
    
    results->predicted_values = (double*)calloc(nind, sizeof(double));
    if (!results->predicted_values) { 
        free(results->posterior_means); 
        free(results->posterior_vars); 
        free(results->distribution_probs);
        free(results->category_probs);
        return ERR_MEMORY; 
    }
    
    return SUCCESS;
}

/**
 * Free MCMCResults arrays.
 */
void free_results(MCMCResults *results) {
    if (!results) return;
    
    if (results->posterior_means) { free(results->posterior_means); results->posterior_means = NULL; }
    if (results->posterior_vars) { free(results->posterior_vars); results->posterior_vars = NULL; }
    if (results->distribution_probs) { free(results->distribution_probs); results->distribution_probs = NULL; }
    if (results->category_probs) { free(results->category_probs); results->category_probs = NULL; }
    if (results->predicted_values) { free(results->predicted_values); results->predicted_values = NULL; }
}

/**
 * Internal structure allocation (used by all model entry points).
 */
static int allocate_internal_data(const MCMCParams *params, int nloci, int nt, int ncat, int nind,
                                   GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    int ndist = params->num_distributions;
    
    /* Allocate GenomicData arrays */
    gdata->snp_correlations = (double*)calloc(nloci, sizeof(double));
    gdata->allele_frequencies = (double*)calloc(nloci, sizeof(double));
    gdata->annotations_per_locus = (int*)calloc(nloci, sizeof(int));
    gdata->current_category = (int*)calloc(nloci, sizeof(int));
    gdata->current_distribution = (int*)calloc(nloci, sizeof(int));
    gdata->distribution_per_category = (int*)calloc(nloci * ncat, sizeof(int));
    gdata->effects_per_category = (double*)calloc(nloci * ncat, sizeof(double));
    gdata->permvec = (int*)calloc(nloci, sizeof(int));
    gdata->permannot = (int*)calloc(ncat, sizeof(int));
    gdata->atemp = (int*)calloc(ncat, sizeof(int));
    gdata->predicted_values = (double*)calloc(nind, sizeof(double));
    
    if (!gdata->snp_correlations || !gdata->annotations_per_locus || 
        !gdata->current_category || !gdata->current_distribution ||
        !gdata->distribution_per_category || !gdata->effects_per_category ||
        !gdata->permvec || !gdata->permannot || !gdata->atemp ||
        !gdata->predicted_values) {
        return ERR_MEMORY;
    }
    
    /* Allocate MCMCState arrays */
    mstate->snp_effects = (double*)calloc(nloci, sizeof(double));
    mstate->adjusted_phenotypes = (double*)calloc(nt, sizeof(double));
    mstate->distribution_variances = (double*)calloc(ndist, sizeof(double));
    mstate->variance_scaling_factors = (double*)calloc(ndist, sizeof(double));
    mstate->dirichlet_priors = (double*)calloc(ndist, sizeof(double));
    mstate->log_distribution_variances = (double*)calloc(ndist, sizeof(double));
    mstate->residual_variance_over_distribution_variances = (double*)calloc(ndist, sizeof(double));
    mstate->p = (double*)calloc(ndist * ncat, sizeof(double));
    mstate->log_p = (double*)calloc(ndist * ncat, sizeof(double));
    mstate->variance_per_distribution = (double*)calloc(ndist * ncat, sizeof(double));
    mstate->snps_per_distribution = (int*)calloc(ndist * ncat, sizeof(int));
    mstate->dirichlet_scratch = (double*)calloc(ndist, sizeof(double));
    mstate->category_probabilities = (double*)calloc(ncat, sizeof(double));
    mstate->ytemp = (double*)calloc(nloci, sizeof(double));
    mstate->category_dirichlet_scratch = (double*)calloc(ncat, sizeof(double));
    mstate->log_likelihoods = (double*)calloc(ndist, sizeof(double));
    mstate->selection_probs = (double*)calloc(ndist, sizeof(double));
    mstate->sstemp = (double*)calloc(ncat, sizeof(double));
    mstate->ss = (double*)calloc(ncat, sizeof(double));
    mstate->z = (double*)calloc(nt, sizeof(double));
    
    if (!mstate->snp_effects || !mstate->adjusted_phenotypes ||
        !mstate->distribution_variances || !mstate->variance_scaling_factors ||
        !mstate->dirichlet_priors || !mstate->log_distribution_variances ||
        !mstate->residual_variance_over_distribution_variances ||
        !mstate->p || !mstate->log_p || !mstate->variance_per_distribution ||
        !mstate->snps_per_distribution || !mstate->dirichlet_scratch ||
        !mstate->category_probabilities || !mstate->ytemp ||
        !mstate->category_dirichlet_scratch || !mstate->log_likelihoods ||
        !mstate->selection_probs || !mstate->sstemp || !mstate->ss || !mstate->z) {
        return ERR_MEMORY;
    }
    
    /* Allocate MCMCStorage arrays */
    mstore->sum_snp_effects = (double*)calloc(nloci, sizeof(double));
    mstore->varustore = (double*)calloc(nloci, sizeof(double));
    mstore->varistore = (double*)calloc(nloci, sizeof(double));
    mstore->mu_vare_store = (double*)calloc(NUM_HYPERPARAMETER_STATS, sizeof(double));
    mstore->sum_snps_per_distribution = (double*)calloc(ndist * ncat, sizeof(double));
    mstore->sum_variance_per_distribution = (double*)calloc(ndist * ncat, sizeof(double));
    mstore->sum_mixture_proportions = (double*)calloc(ndist * ncat, sizeof(double));
    mstore->sum_distribution_counts = (double*)calloc(nloci * ndist, sizeof(double));
    mstore->sum_category_counts = (double*)calloc(nloci * ncat, sizeof(double));
    
    if (!mstore->sum_snp_effects || !mstore->varustore || !mstore->varistore ||
        !mstore->mu_vare_store || !mstore->sum_snps_per_distribution ||
        !mstore->sum_variance_per_distribution || !mstore->sum_mixture_proportions ||
        !mstore->sum_distribution_counts || !mstore->sum_category_counts) {
        return ERR_MEMORY;
    }
    
    return SUCCESS;
}

/**
 * Free internal data structures.
 */
static void free_internal_data(GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    /* Free GenomicData */
    if (gdata->snp_correlations) free(gdata->snp_correlations);
    if (gdata->allele_frequencies) free(gdata->allele_frequencies);
    if (gdata->annotations_per_locus) free(gdata->annotations_per_locus);
    if (gdata->current_category) free(gdata->current_category);
    if (gdata->current_distribution) free(gdata->current_distribution);
    if (gdata->distribution_per_category) free(gdata->distribution_per_category);
    if (gdata->effects_per_category) free(gdata->effects_per_category);
    if (gdata->permvec) free(gdata->permvec);
    if (gdata->permannot) free(gdata->permannot);
    if (gdata->atemp) free(gdata->atemp);
    if (gdata->predicted_values) free(gdata->predicted_values);
    
    /* Free MCMCState */
    if (mstate->snp_effects) free(mstate->snp_effects);
    if (mstate->adjusted_phenotypes) free(mstate->adjusted_phenotypes);
    if (mstate->distribution_variances) free(mstate->distribution_variances);
    if (mstate->variance_scaling_factors) free(mstate->variance_scaling_factors);
    if (mstate->dirichlet_priors) free(mstate->dirichlet_priors);
    if (mstate->log_distribution_variances) free(mstate->log_distribution_variances);
    if (mstate->residual_variance_over_distribution_variances) free(mstate->residual_variance_over_distribution_variances);
    if (mstate->p) free(mstate->p);
    if (mstate->log_p) free(mstate->log_p);
    if (mstate->variance_per_distribution) free(mstate->variance_per_distribution);
    if (mstate->snps_per_distribution) free(mstate->snps_per_distribution);
    if (mstate->dirichlet_scratch) free(mstate->dirichlet_scratch);
    if (mstate->category_probabilities) free(mstate->category_probabilities);
    if (mstate->ytemp) free(mstate->ytemp);
    if (mstate->category_dirichlet_scratch) free(mstate->category_dirichlet_scratch);
    if (mstate->log_likelihoods) free(mstate->log_likelihoods);
    if (mstate->selection_probs) free(mstate->selection_probs);
    if (mstate->sstemp) free(mstate->sstemp);
    if (mstate->ss) free(mstate->ss);
    if (mstate->z) free(mstate->z);
    
    /* Free MCMCStorage */
    if (mstore->sum_snp_effects) free(mstore->sum_snp_effects);
    if (mstore->varustore) free(mstore->varustore);
    if (mstore->varistore) free(mstore->varistore);
    if (mstore->mu_vare_store) free(mstore->mu_vare_store);
    if (mstore->sum_snps_per_distribution) free(mstore->sum_snps_per_distribution);
    if (mstore->sum_variance_per_distribution) free(mstore->sum_variance_per_distribution);
    if (mstore->sum_mixture_proportions) free(mstore->sum_mixture_proportions);
    if (mstore->sum_distribution_counts) free(mstore->sum_distribution_counts);
    if (mstore->sum_category_counts) free(mstore->sum_category_counts);
}

/**
 * Copy results from internal structures to public results struct.
 */
static void copy_results(MCMCResults *results, GenomicData *gdata, MCMCState *mstate, 
                         MCMCStorage *mstore, int nloci, int ndist, int ncat, int nind) {
    /* Copy per-SNP results */
    for (int i = 0; i < nloci; i++) {
        results->posterior_means[i] = mstore->sum_snp_effects[i];
        results->posterior_vars[i] = mstore->varustore[i];
    }
    
    /* Copy distribution probabilities */
    for (int i = 0; i < nloci * ndist; i++) {
        results->distribution_probs[i] = mstore->sum_distribution_counts[i];
    }
    
    /* Copy category probabilities */
    for (int i = 0; i < nloci * ncat; i++) {
        results->category_probs[i] = mstore->sum_category_counts[i];
    }
    
    /* Copy hyperparameters */
    results->mu = mstore->mu_vare_store[0];
    results->nsnp = mstore->mu_vare_store[1];
    results->vara = mstore->mu_vare_store[2];
    results->vare = mstore->mu_vare_store[3];
    
    /* Copy predicted values */
    for (int i = 0; i < nind; i++) {
        results->predicted_values[i] = gdata->predicted_values[i];
    }
}

/* =========================================================================
 * Public API Implementation
 * ========================================================================= */

/**
 * Run BayesRCO mixture model.
 */
int run_bayesrco_mixture(
    const MCMCParams *params,
    const double *genotypes,
    const double *phenotypes,
    const int *categories,
    int nloci, int nt, int ncat, int nind,
    const int *trains,
    uint64_t seed,
    MCMCResults *results
) {
    if (!params || !genotypes || !phenotypes || !results) {
        return ERR_INVALID_ARG;
    }
    
    int ret;
    int ndist = params->num_distributions;
    
    /* Initialize internal structures */
    ModelConfig config;
    GenomicData gdata;
    MCMCState mstate;
    MCMCStorage mstore;
    prng_state rs;
    
    memset(&config, 0, sizeof(ModelConfig));
    memset(&gdata, 0, sizeof(GenomicData));
    memset(&mstate, 0, sizeof(MCMCState));
    memset(&mstore, 0, sizeof(MCMCStorage));
    
    /* Setup config from params */
    config.num_iterations = params->num_iterations;
    config.burnin_iterations = params->burnin_iterations;
    config.thinning_interval = params->thinning_interval;
    config.num_distributions = params->num_distributions;
    config.num_categories = ncat;
    config.permute = params->permute;
    config.mixture = true;
    config.nobayesCpi = true;
    config.VCE = params->VCE;
    config.dfvara = params->dfvara;
    config.dfvare = params->dfvare;
    config.vara_ap = params->vara_ap;
    config.vare_ap = params->vare_ap;
    
    /* Setup gdata */
    gdata.num_loci = nloci;
    gdata.num_individuals = nind;
    gdata.num_phenotyped_individuals = nt;
    gdata.genotypes = (double*)genotypes;  /* Cast away const - internal use only */
    gdata.phenotypes = (double*)phenotypes;
    gdata.categories = (int*)categories;
    gdata.trains = (int*)trains;
    
    /* Allocate internal arrays */
    ret = allocate_internal_data(params, nloci, nt, ncat, nind, &gdata, &mstate, &mstore);
    if (ret != SUCCESS) {
        free_internal_data(&gdata, &mstate, &mstore);
        return ret;
    }
    
    /* Copy variance scaling factors and dirichlet priors */
    for (int i = 0; i < ndist; i++) {
        mstate.variance_scaling_factors[i] = params->variance_scaling_factors[i];
        mstate.dirichlet_priors[i] = params->dirichlet_priors[i];
    }
    mstate.variance_genetic = params->variance_genetic;
    mstate.variance_residual = params->variance_residual;
    mstate.nnind = (double)nt;
    
    /* Initialize RNG */
    if (seed == 0) {
        seed = (uint64_t)time(NULL);
    }
    uint64_t s_packed[4] = {seed, seed ^ 0xDEADBEEF, seed ^ 0x12345678, seed ^ 0x87654321};
    rng_seed(&rs, s_packed);
    
    /* Run MCMC */
    ret = run_mcmc(&config, &gdata, &mstate, &mstore, &rs);
    
    if (ret == SUCCESS) {
        /* Copy results */
        copy_results(results, &gdata, &mstate, &mstore, nloci, ndist, ncat, nind);
    }
    
    /* Note: We don't free gdata.genotypes, gdata.phenotypes, gdata.categories, gdata.trains
     * because they point to caller-provided data */
    gdata.genotypes = NULL;
    gdata.phenotypes = NULL;
    gdata.categories = NULL;
    gdata.trains = NULL;
    
    free_internal_data(&gdata, &mstate, &mstore);
    
    return ret;
}

/**
 * Run BayesRCO additive model.
 */
int run_bayesrco_additive(
    const MCMCParams *params,
    const double *genotypes,
    const double *phenotypes,
    const int *categories,
    int nloci, int nt, int ncat, int nind,
    const int *trains,
    uint64_t seed,
    MCMCResults *results
) {
    if (!params || !genotypes || !phenotypes || !results) {
        return ERR_INVALID_ARG;
    }
    
    int ret;
    int ndist = params->num_distributions;
    
    /* Initialize internal structures */
    ModelConfig config;
    GenomicData gdata;
    MCMCState mstate;
    MCMCStorage mstore;
    prng_state rs;
    
    memset(&config, 0, sizeof(ModelConfig));
    memset(&gdata, 0, sizeof(GenomicData));
    memset(&mstate, 0, sizeof(MCMCState));
    memset(&mstore, 0, sizeof(MCMCStorage));
    
    /* Setup config from params - additive model */
    config.num_iterations = params->num_iterations;
    config.burnin_iterations = params->burnin_iterations;
    config.thinning_interval = params->thinning_interval;
    config.num_distributions = params->num_distributions;
    config.num_categories = ncat;
    config.permute = params->permute;
    config.mixture = false;  /* Additive model */
    config.nobayesCpi = true;
    config.VCE = params->VCE;
    config.dfvara = params->dfvara;
    config.dfvare = params->dfvare;
    config.vara_ap = params->vara_ap;
    config.vare_ap = params->vare_ap;
    
    /* Setup gdata */
    gdata.num_loci = nloci;
    gdata.num_individuals = nind;
    gdata.num_phenotyped_individuals = nt;
    gdata.genotypes = (double*)genotypes;
    gdata.phenotypes = (double*)phenotypes;
    gdata.categories = (int*)categories;
    gdata.trains = (int*)trains;
    
    /* Allocate internal arrays */
    ret = allocate_internal_data(params, nloci, nt, ncat, nind, &gdata, &mstate, &mstore);
    if (ret != SUCCESS) {
        free_internal_data(&gdata, &mstate, &mstore);
        return ret;
    }
    
    /* Copy variance scaling factors and dirichlet priors */
    for (int i = 0; i < ndist; i++) {
        mstate.variance_scaling_factors[i] = params->variance_scaling_factors[i];
        mstate.dirichlet_priors[i] = params->dirichlet_priors[i];
    }
    mstate.variance_genetic = params->variance_genetic;
    mstate.variance_residual = params->variance_residual;
    mstate.nnind = (double)nt;
    
    /* Initialize RNG */
    if (seed == 0) {
        seed = (uint64_t)time(NULL);
    }
    uint64_t s_packed[4] = {seed, seed ^ 0xDEADBEEF, seed ^ 0x12345678, seed ^ 0x87654321};
    rng_seed(&rs, s_packed);
    
    /* Run MCMC */
    ret = run_mcmc(&config, &gdata, &mstate, &mstore, &rs);
    
    if (ret == SUCCESS) {
        copy_results(results, &gdata, &mstate, &mstore, nloci, ndist, ncat, nind);
    }
    
    gdata.genotypes = NULL;
    gdata.phenotypes = NULL;
    gdata.categories = NULL;
    gdata.trains = NULL;
    
    free_internal_data(&gdata, &mstate, &mstore);
    
    return ret;
}

/**
 * Run BayesCpi model.
 */
int run_bayesrco_bayesCpi(
    const MCMCParams *params,
    const double *genotypes,
    const double *phenotypes,
    const int *categories,
    int nloci, int nt, int ncat, int nind,
    const int *trains,
    uint64_t seed,
    MCMCResults *results
) {
    if (!params || !genotypes || !phenotypes || !results) {
        return ERR_INVALID_ARG;
    }
    
    int ret;
    int ndist = params->num_distributions;
    
    /* Initialize internal structures */
    ModelConfig config;
    GenomicData gdata;
    MCMCState mstate;
    MCMCStorage mstore;
    prng_state rs;
    
    memset(&config, 0, sizeof(ModelConfig));
    memset(&gdata, 0, sizeof(GenomicData));
    memset(&mstate, 0, sizeof(MCMCState));
    memset(&mstore, 0, sizeof(MCMCStorage));
    
    /* Setup config from params - BayesCpi model */
    config.num_iterations = params->num_iterations;
    config.burnin_iterations = params->burnin_iterations;
    config.thinning_interval = params->thinning_interval;
    config.num_distributions = params->num_distributions;
    config.num_categories = ncat;
    config.permute = params->permute;
    config.mixture = false;
    config.nobayesCpi = false;  /* Enable BayesCpi */
    config.VCE = params->VCE;
    config.dfvara = params->dfvara;
    config.dfvare = params->dfvare;
    config.vara_ap = params->vara_ap;
    config.vare_ap = params->vare_ap;
    
    /* Setup gdata */
    gdata.num_loci = nloci;
    gdata.num_individuals = nind;
    gdata.num_phenotyped_individuals = nt;
    gdata.genotypes = (double*)genotypes;
    gdata.phenotypes = (double*)phenotypes;
    gdata.categories = (int*)categories;
    gdata.trains = (int*)trains;
    
    /* Allocate internal arrays */
    ret = allocate_internal_data(params, nloci, nt, ncat, nind, &gdata, &mstate, &mstore);
    if (ret != SUCCESS) {
        free_internal_data(&gdata, &mstate, &mstore);
        return ret;
    }
    
    /* Copy variance scaling factors and dirichlet priors */
    for (int i = 0; i < ndist; i++) {
        mstate.variance_scaling_factors[i] = params->variance_scaling_factors[i];
        mstate.dirichlet_priors[i] = params->dirichlet_priors[i];
    }
    mstate.variance_genetic = params->variance_genetic;
    mstate.variance_residual = params->variance_residual;
    mstate.nnind = (double)nt;
    
    /* Initialize RNG */
    if (seed == 0) {
        seed = (uint64_t)time(NULL);
    }
    uint64_t s_packed[4] = {seed, seed ^ 0xDEADBEEF, seed ^ 0x12345678, seed ^ 0x87654321};
    rng_seed(&rs, s_packed);
    
    /* Run MCMC */
    ret = run_mcmc(&config, &gdata, &mstate, &mstore, &rs);
    
    if (ret == SUCCESS) {
        copy_results(results, &gdata, &mstate, &mstore, nloci, ndist, ncat, nind);
    }
    
    gdata.genotypes = NULL;
    gdata.phenotypes = NULL;
    gdata.categories = NULL;
    gdata.trains = NULL;
    
    free_internal_data(&gdata, &mstate, &mstore);
    
    return ret;
}
