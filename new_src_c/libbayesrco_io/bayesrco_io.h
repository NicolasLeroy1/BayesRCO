/**
 * @file bayesrco_io.h
 * @brief I/O and data management API for BayesRCO.
 * 
 * This header provides file I/O, data loading/saving, and PLINK format support.
 * It extends ModelConfig with file paths and handles for I/O operations.
 */

#ifndef BAYESRCO_IO_H
#define BAYESRCO_IO_H

#include "bayesrco_types.h"
#include "rng.h"
#include <stdio.h>

/* =========================================================================
 * I/O Configuration (extends ModelConfig with file I/O)
 * ========================================================================= */

/**
 * @struct IOConfig
 * @brief Configuration with file paths and handles for I/O operations.
 * 
 * Extends the pure ModelConfig from stats library with file I/O fields.
 */
typedef struct {
    /* Statistical configuration (from ModelConfig) */
    ModelConfig model;
    
    /* File paths */
    char geno_file_path[PATH_MAX_LENGTH];
    char pheno_file_path[PATH_MAX_LENGTH];
    char param_file_path[PATH_MAX_LENGTH];
    char model_file_path[PATH_MAX_LENGTH];
    char freq_file_path[PATH_MAX_LENGTH];
    char gv_file_path[PATH_MAX_LENGTH];
    char snp_output_file_path[PATH_MAX_LENGTH];
    char hyp_output_file_path[PATH_MAX_LENGTH];
    char cat_input_file_path[PATH_MAX_LENGTH];
    char cat_output_file_path[PATH_MAX_LENGTH];
    char beta_file_path[PATH_MAX_LENGTH];
    char log_file_path[PATH_MAX_LENGTH];
    
    /* File handles (opened during run, closed at cleanup) */
    FILE *fp_log;    /**< Log file for diagnostics */
    FILE *fp_hyp;    /**< Hyperparameters output file */
    FILE *fp_snp;    /**< SNP effects output file */
    FILE *fp_cat;    /**< Category output file */
    FILE *fp_beta;   /**< Beta effects output file */
    
    /* I/O flags */
    bool mcmc;       /**< MCMC mode (true) vs prediction mode (false) */
    bool cat;        /**< Use category input file */
    bool beta;       /**< Output beta effects */
    
    /* Extra options that don't belong in pure ModelConfig */
    int trait_column_index;  /**< 1-based index of phenotype column */
    uint64_t random_seed;    /**< Random seed */
} IOConfig;

/* =========================================================================
 * PLINK I/O Functions (from libplink)
 * ========================================================================= */

/**
 * Get dimensions from PLINK files.
 * 
 * @param prefix     PLINK file prefix (without extension)
 * @param num_ind    Output: number of individuals (from .fam)
 * @param num_snps   Output: number of SNPs (from .bim)
 * @return          SUCCESS or error code
 */
int plink_get_size(const char *prefix, int *num_ind, int *num_snps);

/**
 * Load phenotypes from PLINK .fam file.
 * 
 * @param prefix      PLINK file prefix
 * @param phenotypes  Output array (caller-allocated, size: num_ind)
 * @param num_ind     Number of individuals
 * @return           SUCCESS or error code
 */
int plink_load_phenotypes(const char *prefix, double *phenotypes, int num_ind);

/**
 * Load genotypes from PLINK .bed file.
 * 
 * @param prefix     PLINK file prefix
 * @param genotypes  Output array (caller-allocated, size: num_ind * num_snps)
 * @param num_ind    Number of individuals
 * @param num_snps   Number of SNPs
 * @return          SUCCESS or error code
 */
int plink_load_genotypes(const char *prefix, double *genotypes, int num_ind, int num_snps);

/* =========================================================================
 * Data Management Functions
 * ========================================================================= */

/**
 * Allocate GenomicData and related structures.
 * 
 * @param ioconfig   I/O configuration
 * @param gdata      Genomic data to allocate
 * @param mstate     MCMC state to allocate
 * @param mstore     MCMC storage to allocate
 * @return          SUCCESS or error code
 */
int io_allocate_data(IOConfig *ioconfig, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);

/**
 * Load parameters and model from files.
 */
int io_load_params(IOConfig *ioconfig, GenomicData *gdata, MCMCState *mstate);

/**
 * Load annotation categories from file.
 */
int io_load_categories(IOConfig *ioconfig, GenomicData *gdata);

/**
 * Write DGV (predicted genomic values) to file.
 */
int io_write_dgv(IOConfig *ioconfig, GenomicData *gdata);

/**
 * Write model output (parameters and effects) to files.
 */
int io_output_model(IOConfig *ioconfig, GenomicData *gdata, MCMCStorage *mstore);

/**
 * Write MCMC results to parameter and model files.
 * 
 * @param ioconfig   I/O configuration
 * @param results    Results from MCMC run
 * @param nloci      Number of SNPs
 * @param nind       Number of individuals (for phenotypic predictions)
 * @return          SUCCESS or error code
 */
int io_write_results(IOConfig *ioconfig, MCMCResults *results, int nloci, int nind);

/**
 * Get dimensions using IOConfig.
 */
int io_get_size(IOConfig *ioconfig, GenomicData *gdata);


/**
 * Write allele frequencies to file.
 */
int io_write_frequencies(IOConfig *ioconfig, GenomicData *gdata);

/**
 * Load allele frequencies from file.
 */
int io_load_frequencies(IOConfig *ioconfig, GenomicData *gdata);

/**
 * Write predicted genomic values (DGV) to file.
 */
int io_write_predictions(IOConfig *ioconfig, double *predicted_values, int nind, int *trains);

/**
 * Cleanup all allocated data and close files.
 */
void io_cleanup(IOConfig *ioconfig, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);

/**
 * Load phenotypes using IOConfig.
 * Uses ioconfig->pheno_file_path and ioconfig->trait_column_index.
 */
int io_load_phenotypes(IOConfig *ioconfig, GenomicData *gdata);

/**
 * Load genotypes using IOConfig.
 * Uses ioconfig->geno_file_path.
 */
int io_load_genotypes(IOConfig *ioconfig, GenomicData *gdata);

/* =========================================================================
 * Convenience Wrappers
 * ========================================================================= */

/**
 * Load and standardize genotypes using PLINK files.
 * 
 * High-level wrapper that:
 * 1. Loads genotypes from PLINK .bed
 * 2. Loads/saves allele frequencies
 * 3. Standardizes genotypes
 * 
 * @param ioconfig   I/O configuration
 * @param gdata      Genomic data (modified in place)
 * @return          SUCCESS or error code
 */
int io_load_and_standardize_plink(IOConfig *ioconfig, GenomicData *gdata);

#endif /* BAYESRCO_IO_H */
