/**
 * @file utils.c
 * @brief General utility functions for data processing and computation.
 */

#include "utils.h"
#include <stdlib.h>
#include <math.h>

/**
 * Permute an array of indices using Fisher-Yates shuffle.
 * 
 * @param rs  Random number generator state
 * @param n   Size of array
 * @param p   Array to permute (will contain 0..n-1 in random order)
 */
void permutate(prng_state *rs, int n, int *p) {
    int i, j, k, ipj, itemp, m;
    double u[PERMUTE_BATCH_SIZE];
    
    /* Initialize to identity permutation */
    for (i = 0; i < n; i++) {
        p[i] = i; 
    }
    
    /* Generate random numbers in batches for efficiency */
    for (i = 0; i < n; i += PERMUTE_BATCH_SIZE) {
        m = (n - i < PERMUTE_BATCH_SIZE) ? (n - i) : PERMUTE_BATCH_SIZE;
        for (int x = 0; x < PERMUTE_BATCH_SIZE; x++) {
            u[x] = rng_uniform(rs, 0.0, 1.0);
        }
        
        for (j = 0; j < m; j++) {
            ipj = i + j;
            k = (int)(u[j] * (n - ipj)) + ipj;
            if (k >= n) k = n - 1;

            itemp = p[ipj];
            p[ipj] = p[k];
            p[k] = itemp;
        }
    }
}

/**
 * Compute direct genomic values (predicted breeding values) for all individuals.
 * 
 * DGV = mu + sum(X * g) where g is the posterior mean of SNP effects.
 * 
 * @param gdata   Genomic data (genotypes)
 * @param mstate  MCMC state (intercept)
 * @param mstore  MCMC storage (posterior mean effects)
 */
void compute_dgv(GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    int nind = gdata->num_individuals;
    int nloci = gdata->num_loci;
    int nt = gdata->num_phenotyped_individuals;
    int tr = 0;
    
    /* Initialize predictions to missing */
    for (int i = 0; i < nind; i++) {
        gdata->predicted_values[i] = MISSING_VALUE;
    }

    /* Compute DGV for training individuals */
    for (int i = 0; i < nind; i++) {
        if (gdata->trains[i] == 0) {
            double dp = 0.0;
            for (int j = 0; j < nloci; j++) {
                /* Column-major access: X[j * nt + tr] */
                dp += gdata->genotypes[j * nt + tr] * mstore->sum_snp_effects[j];
            }
            gdata->predicted_values[i] = mstate->mu + dp;
            tr++;
        }
    }
}

/**
 * Compute residuals (adjusted phenotypes) for MCMC initialization.
 * 
 * residual = y - X*g - mu
 * 
 * @param gdata   Genomic data
 * @param mstate  MCMC state
 */
void compute_residuals(GenomicData *gdata, MCMCState *mstate) {
    int nind = gdata->num_individuals;
    int nloci = gdata->num_loci;
    int nt = gdata->num_phenotyped_individuals;
    int tr = 0;
    
    for (int i = 0; i < nind; i++) {
        if (gdata->trains[i] == 0) {
            double dp = 0.0;
            for (int j = 0; j < nloci; j++) {
                /* Column-major access */
                dp += gdata->genotypes[j * nt + tr] * mstate->snp_effects[j];
            }
            mstate->adjusted_phenotypes[tr] = gdata->phenotypes[i] - dp - mstate->mu;
            tr++;
        }
    }
}

/**
 * Center and standardize genotypes.
 * 
 * Transforms genotypes to have mean 0 and variance proportional to 2pq.
 * Missing values (>2) are imputed with the mean genotype.
 * 
 * @param config  Model configuration
 * @param gdata   Genomic data (modified in place)
 */
void xcenter(ModelConfig *config, GenomicData *gdata) {
    int nloci = gdata->num_loci;
    int nt = gdata->num_phenotyped_individuals;
    double q, qtest;
    FILE *fp;
    
    if (config->mcmc) {
        /* Compute allele frequencies and standardize */
        for (int j = 0; j < nloci; j++) {
            double sum = 0.0;
            int cnt = 0;
            
            /* First pass: compute mean (allele frequency) */
            for (int tr = 0; tr < nt; tr++) {
                double val = gdata->genotypes[j * nt + tr];
                if (val < GENOTYPE_MISSING_THRESHOLD) {
                    sum += val;
                    cnt++;
                }
            }
            
            q = (cnt > 0) ? sum / (2.0 * cnt) : 0.0;

            if (q == 1.0 || q == 0.0) {
                /* Fixed allele - set to zero */
                for (int tr = 0; tr < nt; tr++) {
                    gdata->genotypes[j * nt + tr] = 0.0;
                }
            } else {
                /* Standardize: (x - 2q) / sqrt(2pq) */
                double denom = sqrt(2.0 * q * (1.0 - q));
                for (int tr = 0; tr < nt; tr++) {
                    double val = gdata->genotypes[j * nt + tr];
                    if (val > 2.0) val = 2.0 * q;  /* Impute missing */
                    gdata->genotypes[j * nt + tr] = (val - 2.0 * q) / denom;
                }
            }
            gdata->allele_frequencies[j] = q;
        }
        
        /* Write frequencies to file */
        fp = fopen(config->freq_file_path, "w");
        if (fp) {
            for (int j = 0; j < nloci; j++) {
                fprintf(fp, "%10.6f\n", gdata->allele_frequencies[j]);
            }
            fclose(fp);
        }
        
    } else {
        /* Prediction mode: read frequencies from file */
        fp = fopen(config->freq_file_path, "r");
        if (fp) {
            for (int j = 0; j < nloci; j++) {
                if (fscanf(fp, "%lf", &gdata->allele_frequencies[j]) != 1) break;
            }
            fclose(fp);
        }
        
        /* Standardize using stored frequencies */
        for (int j = 0; j < nloci; j++) {
            q = gdata->allele_frequencies[j];
            if (q == 1.0 || q == 0.0) {
                for (int tr = 0; tr < nt; tr++) {
                    gdata->genotypes[j * nt + tr] = 0.0;
                }
            } else {
                /* Compute current mean for missing imputation */
                double sum = 0.0;
                int cnt = 0;
                for (int tr = 0; tr < nt; tr++) {
                    double val = gdata->genotypes[j * nt + tr];
                    if (val < GENOTYPE_MISSING_THRESHOLD) {
                        sum += val;
                        cnt++;
                    }
                }
                qtest = (cnt > 0) ? sum / (2.0 * cnt) : 0.0;
                
                double denom = sqrt(2.0 * q * (1.0 - q));
                for (int tr = 0; tr < nt; tr++) {
                    double val = gdata->genotypes[j * nt + tr];
                    if (val > 2.0) val = 2.0 * qtest;
                    gdata->genotypes[j * nt + tr] = (val - 2.0 * q) / denom;
                }
            }
        }
    }
}
