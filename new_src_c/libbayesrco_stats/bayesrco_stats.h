/**
 * @file bayesrco_stats.h
 * @brief Public API for BayesRCO statistical library.
 * 
 * This header provides the main entry points for running BayesRCO models.
 * Users should call one of the run_bayesrco_* functions after loading data.
 */

#ifndef BAYESRCO_STATS_H
#define BAYESRCO_STATS_H

#include "bayesrco_types.h"

/* =========================================================================
 * Public API - Model Entry Points
 * ========================================================================= */

/**
 * Run BayesRCO mixture model.
 * 
 * The mixture model allows SNPs to belong to multiple annotation categories,
 * and samples the category assignment during MCMC.
 * 
 * @param params     MCMC parameters (iterations, hyperparameters, etc.)
 * @param genotypes  Standardized genotype matrix, column-major (nloci x nt)
 * @param phenotypes Phenotype values for training individuals (size: nt)
 * @param categories Annotation matrix, flattened row-major (nloci x ncat)
 * @param nloci      Number of SNP loci
 * @param nt         Number of training individuals
 * @param ncat       Number of annotation categories
 * @param nind       Total number of individuals (for predictions)
 * @param trains     Training indicator array (0=train, 1=test), size: nind
 * @param seed       Random seed for reproducibility (0 for time-based)
 * @param results    Output structure (caller allocates arrays)
 * @return           SUCCESS (0) on success, error code otherwise
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
);

/**
 * Run BayesRCO additive model.
 * 
 * The additive model sums effects across all annotation categories
 * rather than selecting a single category per SNP.
 * 
 * @param params     MCMC parameters
 * @param genotypes  Standardized genotype matrix, column-major (nloci x nt)
 * @param phenotypes Phenotype values for training individuals (size: nt)
 * @param categories Annotation matrix, flattened row-major (nloci x ncat)
 * @param nloci      Number of SNP loci
 * @param nt         Number of training individuals
 * @param ncat       Number of annotation categories
 * @param nind       Total number of individuals
 * @param trains     Training indicator array (0=train, 1=test), size: nind
 * @param seed       Random seed
 * @param results    Output structure
 * @return           SUCCESS (0) on success, error code otherwise
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
);

/**
 * Run BayesCpi model.
 * 
 * The BayesCpi model estimates the mixture proportion (pi) from the data
 * rather than using fixed priors.
 * 
 * @param params     MCMC parameters
 * @param genotypes  Standardized genotype matrix, column-major (nloci x nt)
 * @param phenotypes Phenotype values for training individuals (size: nt)
 * @param categories Annotation matrix, flattened row-major (nloci x ncat)
 * @param nloci      Number of SNP loci
 * @param nt         Number of training individuals
 * @param ncat       Number of annotation categories
 * @param nind       Total number of individuals
 * @param trains     Training indicator array (0=train, 1=test), size: nind
 * @param seed       Random seed
 * @param results    Output structure
 * @return           SUCCESS (0) on success, error code otherwise
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
);

/* =========================================================================
 * Helper Functions
 * ========================================================================= */

/**
 * Allocate MCMCResults arrays.
 * 
 * @param results Pointer to results structure
 * @param nloci   Number of SNP loci
 * @param ndist   Number of distributions
 * @param ncat    Number of categories
 * @param nind    Total number of individuals
 * @return        SUCCESS on success, ERR_MEMORY on allocation failure
 */
int allocate_results(MCMCResults *results, int nloci, int ndist, int ncat, int nind);

/**
 * Free MCMCResults arrays.
 * 
 * @param results Pointer to results structure
 */
void free_results(MCMCResults *results);

#endif /* BAYESRCO_STATS_H */
