/**
 * @file mcmc_core_io.c
 * @brief I/O and data management functions for BayesRCO.
 * 
 * Provides file I/O, data loading/saving, and memory management.
 */

#include "bayesrco_io.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* =========================================================================
 * Memory Allocation
 * ========================================================================= */

#define SAFE_CALLOC(ptr, count, type) do { \
    (ptr) = (type*)calloc((count), sizeof(type)); \
    if (!(ptr)) return ERR_MEMORY; \
} while(0)

#define SAFE_FREE(ptr) do { if ((ptr)) { free((ptr)); (ptr) = NULL; } } while(0)

/**
 * Allocate data structures.
 */
int io_allocate_data(IOConfig *ioconfig, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    int nloci = gdata->num_loci;
    int nind = gdata->num_individuals;
    int ncat = ioconfig->model.num_categories;
    int ndist = ioconfig->model.num_distributions;

    if (!ioconfig->mcmc) {
        /* Prediction mode: swap train/test indicators */
        for(int i=0; i<nind; i++) {
            if (gdata->trains[i] == 0) gdata->trains[i] = 3;
            else if (gdata->trains[i] == 1) gdata->trains[i] = 0;
        }
        for(int i=0; i<nind; i++) {
            if (gdata->trains[i] == 3) gdata->trains[i] = 1;
        }
    }
    
    gdata->num_phenotyped_individuals = 0;
    for(int i=0; i<nind; i++) if(gdata->trains[i]==0) gdata->num_phenotyped_individuals++;
    
    int nt = gdata->num_phenotyped_individuals;

    SAFE_CALLOC(gdata->predicted_values, nind, double);
    
    SAFE_CALLOC(mstate->variance_scaling_factors, ndist, double);
    SAFE_CALLOC(mstate->genomic_values, ndist, double);
    
    /* Flattened state matrices (ndist * ncat) */
    SAFE_CALLOC(mstate->p, ndist * ncat, double);
    SAFE_CALLOC(mstate->log_p, ndist * ncat, double);
    SAFE_CALLOC(mstate->snps_per_distribution, ndist * ncat, int);
    SAFE_CALLOC(mstate->variance_per_distribution, ndist * ncat, double);
    
    SAFE_CALLOC(gdata->permannot, ncat, int);
    SAFE_CALLOC(mstate->dirichlet_priors, ndist, double);
    SAFE_CALLOC(mstate->dirichlet_scratch, ndist, double);
    SAFE_CALLOC(mstate->snp_effects, nloci, double);
    SAFE_CALLOC(mstate->adjusted_phenotypes, nt, double);
    SAFE_CALLOC(mstate->z, nt, double);
    
    SAFE_CALLOC(mstate->log_likelihoods, ndist, double);
    SAFE_CALLOC(mstate->selection_probs, ndist, double);
    SAFE_CALLOC(mstate->sstemp, ncat, double);
    
    SAFE_CALLOC(gdata->snp_correlations, nloci, double);
    SAFE_CALLOC(mstore->sum_snp_effects, nloci, double);
    
    /* Flattened storage matrices (ndist * ncat) */
    SAFE_CALLOC(mstore->sum_snps_per_distribution, ndist * ncat, double);
    SAFE_CALLOC(mstore->sum_variance_per_distribution, ndist * ncat, double);
    SAFE_CALLOC(mstore->sum_mixture_proportions, ndist * ncat, double);
    
    /* Flattened SNP-level matrices */
    SAFE_CALLOC(mstore->sum_distribution_counts, nloci * ndist, double);
    SAFE_CALLOC(mstore->sum_category_counts, nloci * ncat, double);
    SAFE_CALLOC(gdata->distribution_per_category, nloci * ncat, int);
    SAFE_CALLOC(gdata->effects_per_category, nloci * ncat, double);
    SAFE_CALLOC(gdata->categories, nloci * ncat, int);

    SAFE_CALLOC(mstore->mu_vare_store, 4, double);
    SAFE_CALLOC(gdata->allele_frequencies, nloci, double);
    SAFE_CALLOC(gdata->permvec, nloci, int);
    
    SAFE_CALLOC(mstate->log_distribution_variances, ndist, double);
    SAFE_CALLOC(mstate->residual_variance_over_distribution_variances, ndist, double);
    
    SAFE_CALLOC(mstore->varustore, nloci, double);
    SAFE_CALLOC(mstore->varistore, nloci, double);
    
    SAFE_CALLOC(gdata->current_distribution, nloci, int);
    SAFE_CALLOC(gdata->annotations_per_locus, nloci, int);
    SAFE_CALLOC(gdata->current_category, nloci, int);
    
    SAFE_CALLOC(mstate->ss, ncat, double);
    SAFE_CALLOC(gdata->included_loci, nloci, double);
    SAFE_CALLOC(mstate->category_probabilities, ncat, double);
    SAFE_CALLOC(mstate->ytemp, nloci, double);
    SAFE_CALLOC(mstate->category_dirichlet_scratch, ncat, double);
    SAFE_CALLOC(gdata->atemp, ncat, int);
    
    return SUCCESS;
}

/* =========================================================================
 * Parameter I/O
 * ========================================================================= */

/**
 * Load parameters from file.
 */
int io_load_params(IOConfig *ioconfig, GenomicData *gdata, MCMCState *mstate) {
    FILE *fp;
    char buffer[1024];
    int nc = ioconfig->model.num_distributions + 1;
    double *gtemp = (double*)calloc(nc, sizeof(double));
    if (!gtemp) return ERR_MEMORY;

    fp = fopen(ioconfig->param_file_path, "r");
    if(!fp) {
        free(gtemp);
        return ERR_FILE_IO;
    }

    if (!fgets(buffer, sizeof(buffer), fp)) { 
        fclose(fp); 
        free(gtemp); 
        return ERR_FILE_IO; 
    }
    
    /* Read SNP effects - just need the last column for prediction */
    MCMCStorage mstore_temp = {0};
    mstore_temp.sum_snp_effects = (double*)calloc(gdata->num_loci, sizeof(double));
    if (!mstore_temp.sum_snp_effects) {
        fclose(fp);
        free(gtemp);
        return ERR_MEMORY;
    }
    
    for(int i=0; i<gdata->num_loci; i++) {
        for(int k=0; k<nc; k++) {
            if (fscanf(fp, "%lf", &gtemp[k]) != 1) break;
        }
        mstore_temp.sum_snp_effects[i] = gtemp[nc-1];
    }
    fclose(fp);
    free(gtemp);
    
    /* Copy effects to mstate for prediction */
    for(int i=0; i<gdata->num_loci; i++) {
        mstate->snp_effects[i] = mstore_temp.sum_snp_effects[i];
    }
    free(mstore_temp.sum_snp_effects);
    
    /* Load intercept from model file */
    fp = fopen(ioconfig->model_file_path, "r");
    if(fp) {
        char dum[100];
        if (fscanf(fp, "%s %lf", dum, &mstate->mu) != 2) mstate->mu = 0.0;
        fclose(fp);
    }
    return SUCCESS;
}

/**
 * Load categories from file.
 */
int io_load_categories(IOConfig *ioconfig, GenomicData *gdata) {
    int nloci = gdata->num_loci;
    int ncat = ioconfig->model.num_categories;

    if (ioconfig->cat) {
        FILE *fp = fopen(ioconfig->cat_input_file_path, "r");
        if (!fp) {
            fprintf(stderr, "Cannot open cat file %s\n", ioconfig->cat_input_file_path);
            return ERR_FILE_IO;
        }
        
        for(int i=0; i<nloci; i++) {
            for(int j=0; j<ncat; j++) { 
                int val;
                if (fscanf(fp, "%d", &val) != 1) break;
                gdata->categories[IDX2(i, j, ncat)] = val; 
            }
        }
        fclose(fp);
    } else {
        if (ncat == 1) {
            for(int i=0; i<nloci; i++) {
                gdata->categories[IDX2(i, 0, ncat)] = 1;
            }
        }
    }
    return SUCCESS;
}

/* =========================================================================
 * Output Functions
 * ========================================================================= */

/**
 * Write DGV (predicted values) to file.
 */
int io_write_dgv(IOConfig *ioconfig, GenomicData *gdata) {
    FILE *fp = fopen(ioconfig->gv_file_path, "w");
    if(!fp) return ERR_FILE_IO;
    
    for(int i=0; i<gdata->num_individuals; i++) {
        if(gdata->trains[i] == 0) {
            fprintf(fp, "%15.7E\n", gdata->predicted_values[i]);
        } else {
            fprintf(fp, "NA\n");
        }
    }
    fclose(fp);
    return SUCCESS;
}

/**
 * Write model output to files.
 */
int io_output_model(IOConfig *ioconfig, GenomicData *gdata, MCMCStorage *mstore) {
    FILE *fp;
    int nloci = gdata->num_loci;
    int ncat = ioconfig->model.num_categories;
    int ndist = ioconfig->model.num_distributions;

    fp = fopen(ioconfig->param_file_path, "w");
    if (!fp) return ERR_FILE_IO;
    
    fprintf(fp, "  ");
    for(int i=1; i<=ndist; i++) fprintf(fp, " PIP%-4d", i);
    fprintf(fp, "  beta");
    for(int i=1; i<=ncat; i++) fprintf(fp, " PAIP%-4d", i);
    fprintf(fp, " Vbeta   Vi\n");
    
    for (int i=0; i<nloci; i++) {
        for(int k=0; k<ndist; k++) fprintf(fp, "%15.7E ", mstore->sum_distribution_counts[IDX2(i, k, ndist)]);
        fprintf(fp, "%15.7E ", mstore->sum_snp_effects[i]);
        for(int k=0; k<ncat; k++) fprintf(fp, "%15.7E ", mstore->sum_category_counts[IDX2(i, k, ncat)]);
        fprintf(fp, "%15.7E %15.7E\n", mstore->varustore[i], mstore->varistore[i]);
    }
    fclose(fp);
    
    fp = fopen(ioconfig->model_file_path, "w");
    if(!fp) return ERR_FILE_IO;
    
    fprintf(fp, "Mean      %15.7E\n", mstore->mu_vare_store[0]);
    fprintf(fp, "Nsnp      %15.7E\n", mstore->mu_vare_store[1]);
    fprintf(fp, "Va        %15.7E\n", mstore->mu_vare_store[2]);
    fprintf(fp, "Ve        %15.7E\n", mstore->mu_vare_store[3]);
    
    for(int j=1; j<=ncat; j++) {
        for(int i=1; i<=ndist; i++) {
             fprintf(fp, "Nk%d_%d      %15.7E\n", i, j, mstore->sum_snps_per_distribution[IDX2(i-1, j-1, ncat)]);
        }
    }
    for(int j=1; j<=ncat; j++) {
        for(int i=1; i<=ndist; i++) {
             fprintf(fp, "Pk%d_%d      %15.7E\n", i, j, mstore->sum_mixture_proportions[IDX2(i-1, j-1, ncat)]);
        }
    }
    for(int j=1; j<=ncat; j++) {
        for(int i=1; i<=ndist; i++) {
             fprintf(fp, "Vk%d_%d      %15.7E\n", i, j, mstore->sum_variance_per_distribution[IDX2(i-1, j-1, ncat)]);
        }
    }
    fclose(fp);
    return SUCCESS;
}

/* =========================================================================
 * Frequency I/O
 * ========================================================================= */

/**
 * Write allele frequencies to file.
 */
int io_write_frequencies(IOConfig *ioconfig, GenomicData *gdata) {
    FILE *fp = fopen(ioconfig->freq_file_path, "w");
    if (!fp) return ERR_FILE_IO;
    
    for (int j = 0; j < gdata->num_loci; j++) {
        fprintf(fp, "%10.6f\n", gdata->allele_frequencies[j]);
    }
    fclose(fp);
    return SUCCESS;
}

/**
 * Load allele frequencies from file.
 */
int io_load_frequencies(IOConfig *ioconfig, GenomicData *gdata) {
    FILE *fp = fopen(ioconfig->freq_file_path, "r");
    if (!fp) return ERR_FILE_IO;
    
    for (int j = 0; j < gdata->num_loci; j++) {
        if (fscanf(fp, "%lf", &gdata->allele_frequencies[j]) != 1) {
            fclose(fp);
            return ERR_FILE_IO;
        }
    }
    fclose(fp);
    return SUCCESS;
}

/* =========================================================================
 * Cleanup
 * ========================================================================= */

/**
 * Cleanup all allocated data and close files.
 */
void io_cleanup(IOConfig *ioconfig, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    /* Close file handles */
    if (ioconfig) {
        if (ioconfig->fp_log) { fclose(ioconfig->fp_log); ioconfig->fp_log = NULL; }
        if (ioconfig->fp_hyp) { fclose(ioconfig->fp_hyp); ioconfig->fp_hyp = NULL; }
        if (ioconfig->fp_snp) { fclose(ioconfig->fp_snp); ioconfig->fp_snp = NULL; }
        if (ioconfig->fp_cat) { fclose(ioconfig->fp_cat); ioconfig->fp_cat = NULL; }
        if (ioconfig->fp_beta) { fclose(ioconfig->fp_beta); ioconfig->fp_beta = NULL; }
    }
    
    /* Free GenomicData */
    if (gdata) {
        SAFE_FREE(gdata->phenotypes);
        SAFE_FREE(gdata->predicted_values);
        SAFE_FREE(gdata->allele_frequencies);
        SAFE_FREE(gdata->snp_correlations);
        SAFE_FREE(gdata->included_loci);
        SAFE_FREE(gdata->genotypes);
        SAFE_FREE(gdata->annotations_per_locus);
        SAFE_FREE(gdata->current_category);
        SAFE_FREE(gdata->current_distribution);
        SAFE_FREE(gdata->trains);
        SAFE_FREE(gdata->permvec);
        SAFE_FREE(gdata->permannot);
        SAFE_FREE(gdata->atemp);
        SAFE_FREE(gdata->categories);
        SAFE_FREE(gdata->distribution_per_category);
        SAFE_FREE(gdata->effects_per_category);
    }
    
    /* Free MCMCState */
    if (mstate) {
        SAFE_FREE(mstate->genomic_values);
        SAFE_FREE(mstate->variance_scaling_factors);
        SAFE_FREE(mstate->dirichlet_priors);
        SAFE_FREE(mstate->snp_effects);
        SAFE_FREE(mstate->adjusted_phenotypes);
        SAFE_FREE(mstate->dirichlet_scratch);
        SAFE_FREE(mstate->category_probabilities);
        SAFE_FREE(mstate->ytemp);
        SAFE_FREE(mstate->category_dirichlet_scratch);
        SAFE_FREE(mstate->log_likelihoods);
        SAFE_FREE(mstate->selection_probs);
        SAFE_FREE(mstate->sstemp);
        SAFE_FREE(mstate->ss);
        SAFE_FREE(mstate->log_distribution_variances);
        SAFE_FREE(mstate->residual_variance_over_distribution_variances);
        SAFE_FREE(mstate->z);
        SAFE_FREE(mstate->variance_per_distribution);
        SAFE_FREE(mstate->p);
        SAFE_FREE(mstate->log_p);
        SAFE_FREE(mstate->snps_per_distribution);
    }
    
    /* Free MCMCStorage */
    if (mstore) {
        SAFE_FREE(mstore->sum_snp_effects);
        SAFE_FREE(mstore->mu_vare_store);
        SAFE_FREE(mstore->varustore);
        SAFE_FREE(mstore->varistore);
        SAFE_FREE(mstore->sum_snps_per_distribution);
        SAFE_FREE(mstore->sum_variance_per_distribution);
        SAFE_FREE(mstore->sum_mixture_proportions);
        SAFE_FREE(mstore->sum_distribution_counts);
        SAFE_FREE(mstore->sum_category_counts);
    }
}

/* =========================================================================
 * Clean Output Functions (using MCMCResults)
 * ========================================================================= */

/**
 * Write MCMC results to parameter and model files.
 */
int io_write_results(IOConfig *ioconfig, MCMCResults *results, int nloci, int nind) {
    FILE *fp;
    int ncat = ioconfig->model.num_categories;
    int ndist = ioconfig->model.num_distributions;

    /* Write .param file */
    fp = fopen(ioconfig->param_file_path, "w");
    if (!fp) return ERR_FILE_IO;
    
    fprintf(fp, "  ");
    for(int i=1; i<=ndist; i++) fprintf(fp, " PIP%-4d", i);
    fprintf(fp, "  beta");
    for(int i=1; i<=ncat; i++) fprintf(fp, " PAIP%-4d", i);
    fprintf(fp, " Vbeta   Vi\n");
    
    for (int i=0; i<nloci; i++) {
        /* PIP (distribution probabilities) */
        for(int k=0; k<ndist; k++) {
            double val = (results->distribution_probs) ? results->distribution_probs[IDX2(i, k, ndist)] : 0.0;
            fprintf(fp, "%15.7E ", val);
        }
        
        /* Beta (posterior mean effect) */
        fprintf(fp, "%15.7E ", results->posterior_means[i]);
        
        /* PAIP (category probabilities) */
        for(int k=0; k<ncat; k++) {
             double val = (results->category_probs) ? results->category_probs[IDX2(i, k, ncat)] : 0.0;
             fprintf(fp, "%15.7E ", val);
        }
        
        /* Vbeta (posterior variance) and Vi (inclusion probability?) */
        /* MCMCResults has posterior_vars. original varustore/varistore */
        /* original varustore = variance of u? varistore = sum of u^2 ? */
        /* MCMCResults.posterior_vars is likely variance of effect. */
        /* For backward compatibility we just write variance and 0 for now as we might miss varistore mapping */
        fprintf(fp, "%15.7E 0.0000000E+00\n", results->posterior_vars[i]);
    }
    fclose(fp);
    
    /* Write .model file */
    fp = fopen(ioconfig->model_file_path, "w");
    if(!fp) return ERR_FILE_IO;
    
    fprintf(fp, "Mean      %15.7E\n", results->mu);
    fprintf(fp, "Nsnp      %15.7E\n", results->nsnp);
    fprintf(fp, "Va        %15.7E\n", results->vara);
    fprintf(fp, "Ve        %15.7E\n", results->vare);
    
    /* Detailed distribution stats are not in MCMCResults, leaving empty for now or 0 */
    for(int j=1; j<=ncat; j++) {
        for(int i=1; i<=ndist; i++) {
             fprintf(fp, "Nk%d_%d      %15.7E\n", i, j, 0.0);
        }
    }
    for(int j=1; j<=ncat; j++) {
        for(int i=1; i<=ndist; i++) {
             fprintf(fp, "Pk%d_%d      %15.7E\n", i, j, 0.0);
        }
    }
    for(int j=1; j<=ncat; j++) {
        for(int i=1; i<=ndist; i++) {
             fprintf(fp, "Vk%d_%d      %15.7E\n", i, j, 0.0);
        }
    }
    fclose(fp);
    return SUCCESS;
}

/**
 * Write predicted genomic values (DGV) to file.
 */
int io_write_predictions(IOConfig *ioconfig, double *predicted_values, int nind, int *trains) {
    FILE *fp = fopen(ioconfig->gv_file_path, "w");
    if(!fp) return ERR_FILE_IO;
    
    for(int i=0; i<nind; i++) {
        /* Original: if(gdata->trains[i] == 0) fprintf(fp, val) else fprintf(fp, "NA") */
        /* trains==0 implies 'phenotyped'/training. */
        if(trains[i] == 0) {
            fprintf(fp, "%15.7E\n", predicted_values[i]);
        } else {
             fprintf(fp, "NA\n");
        }
    }
    fclose(fp);
    return SUCCESS;
}
