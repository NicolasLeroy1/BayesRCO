#ifndef BAYESRCO_H
#define BAYESRCO_H

#include <stdio.h>
#include <stdbool.h>
#include "rng.h"

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
#define PATH_MAX_LENGTH 4096 /* Increased for safety */
#define LINE_BUFFER_SIZE 4096
#define SMALL_BUFFER_SIZE 1024

/* Indexing Macros */
#define IDX2(i, j, nj) ((i)*(nj) + (j))

/* =========================================================================
 * Data Structures
 * ========================================================================= */

/**
 * @struct ModelConfig
 * @brief Configuration parameters for the MCMC model.
 */
typedef struct {
    /* Algorithm parameters */
    int num_iterations;      /* Total number of MCMC iterations (numit) */
    int burnin_iterations;   /* Number of burn-in iterations (burnin) */
    int thinning_interval;   /* Thinning interval for output (thin) */
    int num_distributions;   /* Number of effect size distributions (ndist) */
    int random_seed;         /* Random seed for reproducibility */
    int trait_column_index;  /* Column index for phenotype trait (n) */
    int num_categories;      /* Number of annotation categories (ncat) */
    int marker_set_size;     /* Size of marker set for permutation (msize) */
    int marker_replicates;   /* Number of marker replicates (mrep) */
    
    /* Flags */
    int use_indistinguishable_distributions;
    int burn;
    int use_annotations;
    int unit_log, unit_hyp, unit_loc, unit_cat, unit_beta;
    bool mcmc;               /* True for MCMC, false for prediction */
    bool snpout;             /* Output SNP-level results */
    bool permute;            /* Permute marker order */
    bool cat;                /* Use annotation categories */
    bool beta;               /* Output beta values */
    bool mixture;            /* Use mixture model (vs additive) */
    bool nobayesCpi;         /* Disable BayesCpi (use fixed pi) */
    bool VCE;                /* Variance component estimation */
    
    /* Hyperparameters */
    double dfvara;           /* Degrees of freedom for genetic variance prior */
    double dfvare;           /* Degrees of freedom for residual variance prior */
    double vara_ap;          /* A priori genetic variance */
    double vare_ap;          /* A priori residual variance */
    
    /* File paths */
    char genotype_file_path[PATH_MAX_LENGTH];
    char phenotype_file_path[PATH_MAX_LENGTH];
    char bim_file_path[PATH_MAX_LENGTH];
    char input_prefix[PATH_MAX_LENGTH];
    char output_prefix[PATH_MAX_LENGTH];
    char log_file_path[PATH_MAX_LENGTH];
    char freq_file_path[PATH_MAX_LENGTH];
    char gv_file_path[PATH_MAX_LENGTH];
    char hyp_file_path[PATH_MAX_LENGTH];
    char snp_file_path[PATH_MAX_LENGTH];
    char model_file_path[PATH_MAX_LENGTH];
    char param_file_path[PATH_MAX_LENGTH];
    char beta_file_path[PATH_MAX_LENGTH];
    char cat_file_path[PATH_MAX_LENGTH];
    char cat_input_file_path[PATH_MAX_LENGTH];
    
    /* File handles */
    FILE *fp_log, *fp_hyp, *fp_snp, *fp_cat, *fp_beta;
} ModelConfig;

/**
 * @struct GenomicData
 * @brief Genomic data including genotypes, phenotypes, and annotations.
 */
typedef struct {
    /* Dimensions */
    int num_loci;                    /* Number of SNP loci (nloci) */
    int num_individuals;             /* Total number of individuals (nind) */
    int num_phenotyped_individuals;  /* Number with phenotype data (nt) */
    int num_training;                /* Number in training set */
    int num_testing;                 /* Number in test set */
    int nref;
    
    /* Phenotype data */
    double *phenotypes;              /* Phenotype values (why), size: nind */
    double *predicted_values;        /* Predicted genomic values (pred), size: nind */
    
    /* Genotype data (column-major: nloci x nt) */
    double *genotypes;               /* Standardized genotype matrix (X), size: nloci * nt */
    double *allele_frequencies;      /* Allele frequencies (freq), size: nloci */
    double *snp_correlations;        /* SNP correlations (xpx), size: nloci */
    double *included_loci;           /* Inclusion indicator (includedloci), size: nloci */
    
    /* Annotation data */
    int *categories;                 /* Annotation matrix (C), size: nloci * ncat [Flattened] */
    int *annotations_per_locus;      /* Number of annotations per SNP (nannot), size: nloci */
    int *current_category;           /* Current category assignment (a), size: nloci */
    int *current_distribution;       /* Current distribution assignment (vsnptrack), size: nloci */
    int *distribution_per_category;  /* Distribution per category (snptracker), size: nloci * ncat [Flattened] */
    double *effects_per_category;    /* Effects per category (gannot), size: nloci * ncat [Flattened] */
    
    /* Permutation and training data */
    int *trains;                     /* Training indicator (0=train, 1=test), size: nind */
    int *permvec;                    /* Permutation vector for SNP order, size: nloci */
    int *permannot;                  /* Permutation vector for categories, size: ncat */
    int *atemp;                      /* Temporary annotation buffer, size: ncat */
    
    /* Prediction metrics */
    double msep, bhat, ahat, corr;
} GenomicData;

/**
 * @struct MCMCState
 * @brief Current state of the MCMC sampler.
 */
typedef struct {
    /* Variance components */
    double mu;                       /* Intercept term */
    double variance_genetic;         /* Genetic variance (vara) */
    double variance_residual;        /* Residual variance (vare) */
    double scale;                    /* Scale parameter */
    double yhat;                     /* Mean phenotype */
    double vary;                     /* Phenotype variance */
    double nnind;                    /* Number of phenotyped individuals (as double) */
    
    /* Iteration tracking */
    int included;                    /* Number of included SNPs (non-zero effect) */
    int counter;                     /* Sample counter (post-burnin) */
    int rep;                         /* Current iteration number */
    
    /* Effect estimates */
    double *snp_effects;             /* SNP effect sizes (g), size: nloci */
    double *adjusted_phenotypes;     /* Adjusted phenotypes (yadj), size: nt */
    double *genomic_values;          /* Genomic values, size: ndist */
    
    /* Prior parameters */
    double *variance_scaling_factors;/* Variance scaling per distribution (gpin), size: ndist */
    double *dirichlet_priors;        /* Dirichlet prior parameters (delta), size: ndist */
    
    /* Distribution state */
    double *log_distribution_variances;                  /* Log of distribution variances (log_gp), size: ndist */
    double *residual_variance_over_distribution_variances; /* vare/gp ratio (vare_gp), size: ndist */
    double *p;                       /* Mixture proportions, size: ndist * ncat [Flattened] */
    double *log_p;                   /* Log mixture proportions, size: ndist * ncat [Flattened] */
    double *variance_per_distribution; /* Variance per distribution (varindist), size: ndist * ncat [Flattened] */
    int *snps_per_distribution;      /* SNP counts per distribution (snpindist), size: ndist * ncat [Flattened] */
    
    /* Scratch space for computations */
    double *dirichlet_scratch;       /* Dirichlet sampling scratch (dirx), size: ndist */
    double *category_probabilities;  /* Category probabilities, size: ncat */
    double *ytemp;                   /* Temporary phenotype vector, size: nloci */
    double *category_dirichlet_scratch; /* Category Dirichlet scratch (dira), size: ncat */
    double *log_likelihoods;         /* Log-likelihoods for distributions (s), size: ndist */
    double *selection_probs;         /* Selection probabilities (stemp), size: ndist */
    double *sstemp;                  /* Category selection scratch, size: ncat */
    double *ss;                      /* Category scores, size: ncat */
    double *z;                       /* Working vector, size: nt */
    
    /* Per-iteration scratch variables */
    double zz;                       /* Current SNP sum of squares (X'X) */
    double zz_vare;                  /* zz / vare ratio */
    double rhs;                      /* Right-hand side (X'y_adj) */
    double lhs;                      /* Left-hand side */
    double v1;                       /* Posterior variance denominator */
    double gk;                       /* Current SNP effect being sampled */
    double logdetV;                  /* Log determinant of V */
    double uhat;                     /* Posterior mean */
    double total_ssq;                /* Total sum of squares */
    double detV;                     /* Determinant of V */
    double maxs, maxtemp;            /* Maximum values for stability */
    double xhat, sk, skk, r, ssculm, clike;  /* Sampling scratch */
} MCMCState;

/**
 * @struct MCMCStorage
 * @brief Accumulated statistics for posterior summary.
 */
typedef struct {
    /* Per-SNP accumulators */
    double *sum_snp_effects;         /* Sum of SNP effects (gstore), size: nloci */
    double *varustore;               /* Variance accumulator, size: nloci */
    double *varistore;               /* Sum of squared effects, size: nloci */
    
    /* Hyperparameter accumulators: [0]=mu, [1]=nsnp, [2]=vara, [3]=vare */
    double *mu_vare_store;           /* Size: NUM_HYPERPARAMETER_STATS (4) */
    
    /* Per-distribution accumulators (ndist * ncat) [Flattened] */
    double *sum_snps_per_distribution;      /* Sum of SNP counts (snpstore) */
    double *sum_variance_per_distribution;  /* Sum of variances (varstore) */
    double *sum_mixture_proportions;        /* Sum of mixture proportions (pstore) */
    
    /* Per-SNP distribution counts */
    double *sum_distribution_counts; /* Size: nloci * ndist [Flattened] */
    double *sum_category_counts;     /* Size: nloci * ncat [Flattened] */
} MCMCStorage;

#endif /* BAYESRCO_H */
