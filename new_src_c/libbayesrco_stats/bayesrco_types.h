/**
 * @file bayesrco_types.h
 * @brief Core data structures for BayesRCO.
 * 
 * This header defines all data types used by the statistical library.
 * It contains no I/O dependencies.
 */

#ifndef BAYESRCO_TYPES_H
#define BAYESRCO_TYPES_H

#include <stdbool.h>
#include <stdint.h>

/* =========================================================================
 * Constants & Error Codes
 * ========================================================================= */

#define PI 3.141592653589793238462
#define LOG_UPPER_LIMIT 700.0
#define MISSING_VALUE -9999.0
#define GENOTYPE_MISSING_THRESHOLD 3.0

/* Error Codes */
#define SUCCESS 0
#define ERR_FILE_IO 1
#define ERR_MEMORY 2
#define ERR_INVALID_ARG 3
#define ERR_RUNTIME 4

/* Named constants for magic numbers */
#define DEFAULT_NUM_DISTRIBUTIONS 4
#define DEFAULT_NUM_ITERATIONS 50000
#define DEFAULT_BURNIN 20000
#define DEFAULT_THINNING 10
#define DEFAULT_MARKER_REPLICATES 5000
#define NUM_HYPERPARAMETER_STATS 4
#define PERMUTE_BATCH_SIZE 100
#define PATH_MAX_LENGTH 4096
#define LINE_BUFFER_SIZE 4096
#define SMALL_BUFFER_SIZE 1024

/* Indexing Macros */
#define IDX2(i, j, nj) ((i)*(nj) + (j))

/* =========================================================================
 * Internal Configuration (used by MCMC internals)
 * ========================================================================= */

/**
 * @struct ModelConfig
 * @brief Internal configuration for MCMC model (no I/O fields).
 * 
 * Used internally by MCMC kernels. Does not contain file handles or paths.
 */
typedef struct {
    /* Algorithm parameters */
    int num_iterations;      /**< Total number of MCMC iterations */
    int burnin_iterations;   /**< Number of burn-in iterations */
    int thinning_interval;   /**< Thinning interval for output */
    int num_distributions;   /**< Number of effect size distributions */
    int num_categories;      /**< Number of annotation categories */
    int marker_set_size;     /**< Size of marker set for permutation */
    int marker_replicates;   /**< Number of marker replicates */
    
    /* Flags */
    bool permute;            /**< Permute marker order each iteration */
    bool mixture;            /**< Use mixture model (vs additive) */
    bool nobayesCpi;         /**< Disable BayesCpi (use fixed pi) */
    bool VCE;                /**< Variance component estimation mode */
    
    /* Hyperparameters */
    double dfvara;           /**< Degrees of freedom for genetic variance prior */
    double dfvare;           /**< Degrees of freedom for residual variance prior */
    double vara_ap;          /**< A priori genetic variance */
    double vare_ap;          /**< A priori residual variance */
} ModelConfig;

/* =========================================================================
 * MCMC Parameters (public API input)
 * ========================================================================= */

/**
 * @struct MCMCParams
 * @brief Parameters for running MCMC models.
 * 
 * This is the primary input structure for the public API.
 */
typedef struct {
    /* Algorithm parameters */
    int num_iterations;      /**< Total number of MCMC iterations */
    int burnin_iterations;   /**< Number of burn-in iterations */
    int thinning_interval;   /**< Thinning interval for output */
    int num_distributions;   /**< Number of effect size distributions */
    int num_categories;      /**< Number of annotation categories */
    
    /* Flags */
    bool permute;            /**< Permute marker order each iteration */
    bool VCE;                /**< Variance component estimation mode */
    
    /* Hyperparameters */
    double dfvara;           /**< Degrees of freedom for genetic variance prior */
    double dfvare;           /**< Degrees of freedom for residual variance prior */
    double vara_ap;          /**< A priori genetic variance */
    double vare_ap;          /**< A priori residual variance */
    double variance_genetic; /**< Initial genetic variance */
    double variance_residual;/**< Initial residual variance */
    
    /* Distribution parameters (caller-allocated, size: num_distributions) */
    double *variance_scaling_factors;  /**< Variance scaling per distribution (gpin) */
    double *dirichlet_priors;          /**< Dirichlet prior parameters (delta) */
} MCMCParams;

/* =========================================================================
 * MCMC Results (public API output)
 * ========================================================================= */

/**
 * @struct MCMCResults
 * @brief Results from MCMC run.
 * 
 * This is the primary output structure from the public API.
 * Caller is responsible for allocating arrays with appropriate sizes.
 */
typedef struct {
    /* Per-SNP results (caller-allocated, size: nloci) */
    double *posterior_means;       /**< Posterior mean of SNP effects */
    double *posterior_vars;        /**< Posterior variance of SNP effects */
    
    /* Per-SNP distribution probabilities (caller-allocated, size: nloci * ndist) */
    double *distribution_probs;    /**< Probability of each distribution per SNP */
    
    /* Per-SNP category probabilities (caller-allocated, size: nloci * ncat) */
    double *category_probs;        /**< Probability of each category per SNP */
    
    /* Hyperparameter estimates */
    double mu;                     /**< Posterior mean of intercept */
    double vara;                   /**< Posterior mean of genetic variance */
    double vare;                   /**< Posterior mean of residual variance */
    double nsnp;                   /**< Posterior mean of included SNPs */
    
    /* Per-individual predictions (caller-allocated, size: nind) */
    double *predicted_values;      /**< Direct genomic values */
} MCMCResults;

/* =========================================================================
 * Internal Data Structures
 * ========================================================================= */

/**
 * @struct GenomicData
 * @brief Genomic data including genotypes, phenotypes, and annotations.
 * 
 * Used internally by MCMC kernels.
 */
typedef struct {
    /* Dimensions */
    int num_loci;                    /**< Number of SNP loci */
    int num_individuals;             /**< Total number of individuals */
    int num_phenotyped_individuals;  /**< Number with phenotype data */
    int num_training;                /**< Number in training set */
    int num_testing;                 /**< Number in test set */
    int nref;
    
    /* Phenotype data */
    double *phenotypes;              /**< Phenotype values, size: nind */
    double *predicted_values;        /**< Predicted genomic values, size: nind */
    
    /* Genotype data (column-major: nloci x nt) */
    double *genotypes;               /**< Standardized genotype matrix, size: nloci * nt */
    double *allele_frequencies;      /**< Allele frequencies, size: nloci */
    double *snp_correlations;        /**< SNP correlations (X'X diagonal), size: nloci */
    double *included_loci;           /**< Inclusion indicator, size: nloci */
    
    /* Annotation data */
    int *categories;                 /**< Annotation matrix (flattened), size: nloci * ncat */
    int *annotations_per_locus;      /**< Number of annotations per SNP, size: nloci */
    int *current_category;           /**< Current category assignment, size: nloci */
    int *current_distribution;       /**< Current distribution assignment, size: nloci */
    int *distribution_per_category;  /**< Distribution per category (flattened), size: nloci * ncat */
    double *effects_per_category;    /**< Effects per category (flattened), size: nloci * ncat */
    
    /* Permutation and training data */
    int *trains;                     /**< Training indicator (0=train, 1=test), size: nind */
    int *permvec;                    /**< Permutation vector for SNP order, size: nloci */
    int *permannot;                  /**< Permutation vector for categories, size: ncat */
    int *atemp;                      /**< Temporary annotation buffer, size: ncat */
    
    /* Prediction metrics */
    double msep, bhat, ahat, corr;
} GenomicData;

/**
 * @struct MCMCState
 * @brief Current state of the MCMC sampler.
 */
typedef struct {
    /* Variance components */
    double mu;                       /**< Intercept term */
    double variance_genetic;         /**< Genetic variance (vara) */
    double variance_residual;        /**< Residual variance (vare) */
    double scale;                    /**< Scale parameter */
    double yhat;                     /**< Mean phenotype */
    double vary;                     /**< Phenotype variance */
    double nnind;                    /**< Number of phenotyped individuals (as double) */
    
    /* Iteration tracking */
    int included;                    /**< Number of included SNPs (non-zero effect) */
    int counter;                     /**< Sample counter (post-burnin) */
    int rep;                         /**< Current iteration number */
    
    /* Effect estimates */
    double *snp_effects;             /**< SNP effect sizes, size: nloci */
    double *adjusted_phenotypes;     /**< Adjusted phenotypes (yadj), size: nt */
    double *genomic_values;          /**< Genomic values, size: ndist */
    
    /* Prior parameters */
    double *variance_scaling_factors;/**< Variance scaling per distribution, size: ndist */
    double *dirichlet_priors;        /**< Dirichlet prior parameters, size: ndist */
    
    /* Distribution state */
    double *log_distribution_variances;                  /**< Log of distribution variances, size: ndist */
    double *residual_variance_over_distribution_variances; /**< vare/gp ratio, size: ndist */
    double *p;                       /**< Mixture proportions (flattened), size: ndist * ncat */
    double *log_p;                   /**< Log mixture proportions (flattened), size: ndist * ncat */
    double *variance_per_distribution; /**< Variance per distribution (flattened), size: ndist * ncat */
    int *snps_per_distribution;      /**< SNP counts per distribution (flattened), size: ndist * ncat */
    
    /* Scratch space for computations */
    double *dirichlet_scratch;       /**< Dirichlet sampling scratch, size: ndist */
    double *category_probabilities;  /**< Category probabilities, size: ncat */
    double *ytemp;                   /**< Temporary phenotype vector, size: nloci */
    double *category_dirichlet_scratch; /**< Category Dirichlet scratch, size: ncat */
    double *log_likelihoods;         /**< Log-likelihoods for distributions, size: ndist */
    double *selection_probs;         /**< Selection probabilities, size: ndist */
    double *sstemp;                  /**< Category selection scratch, size: ncat */
    double *ss;                      /**< Category scores, size: ncat */
    double *z;                       /**< Working vector, size: nt */
    
    /* Per-iteration scratch variables */
    double zz;                       /**< Current SNP sum of squares (X'X) */
    double zz_vare;                  /**< zz / vare ratio */
    double rhs;                      /**< Right-hand side (X'y_adj) */
    double lhs;                      /**< Left-hand side */
    double v1;                       /**< Posterior variance denominator */
    double gk;                       /**< Current SNP effect being sampled */
    double logdetV;                  /**< Log determinant of V */
    double uhat;                     /**< Posterior mean */
    double total_ssq;                /**< Total sum of squares */
    double detV;                     /**< Determinant of V */
    double maxs, maxtemp;            /**< Maximum values for stability */
    double xhat, sk, skk, r, ssculm, clike;  /**< Sampling scratch */
} MCMCState;

/**
 * @struct MCMCStorage
 * @brief Accumulated statistics for posterior summary.
 */
typedef struct {
    /* Per-SNP accumulators */
    double *sum_snp_effects;         /**< Sum of SNP effects, size: nloci */
    double *varustore;               /**< Variance accumulator, size: nloci */
    double *varistore;               /**< Sum of squared effects, size: nloci */
    
    /* Hyperparameter accumulators: [0]=mu, [1]=nsnp, [2]=vara, [3]=vare */
    double *mu_vare_store;           /**< Size: NUM_HYPERPARAMETER_STATS (4) */
    
    /* Per-distribution accumulators (flattened, size: ndist * ncat) */
    double *sum_snps_per_distribution;      /**< Sum of SNP counts */
    double *sum_variance_per_distribution;  /**< Sum of variances */
    double *sum_mixture_proportions;        /**< Sum of mixture proportions */
    
    /* Per-SNP distribution counts */
    double *sum_distribution_counts; /**< Size: nloci * ndist (flattened) */
    double *sum_category_counts;     /**< Size: nloci * ncat (flattened) */
} MCMCStorage;

#endif /* BAYESRCO_TYPES_H */
