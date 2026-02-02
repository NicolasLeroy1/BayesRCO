#include "mcmc_core_io.h"
#include <string.h>
#include <time.h>

void init_random_seed_custom(ModelConfig *config, prng_state *rs) {
    uint64_t seed[4];
    if (config->random_seed != 0) {
        for(int i=0; i<4; i++) seed[i] = abs(config->random_seed) + i;
    } else {
        time_t t = time(NULL);
        for(int i=0; i<4; i++) seed[i] = t + 37 * i;
    }
    rng_seed(rs, seed);
}

int allocate_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    int nloci = gdata->num_loci;
    int nind = gdata->num_individuals;
    int ncat = config->num_categories;
    int ndist = config->num_distributions;

    if (!config->mcmc) {
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

int load_param(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore, MCMCState *mstate) {
    FILE *fp;
    char buffer[1024];
    int nc = config->num_distributions + 1;
    double *gtemp = (double*)calloc(nc, sizeof(double));
    if (!gtemp) return ERR_MEMORY;

    fp = fopen(config->param_file_path, "r");
    if(!fp) {
        free(gtemp);
        return ERR_FILE_IO;
    }

    if (!fgets(buffer, sizeof(buffer), fp)) { 
        fclose(fp); 
        free(gtemp); 
        return ERR_FILE_IO; 
    }
    
    for(int i=0; i<gdata->num_loci; i++) {
        for(int k=0; k<nc; k++) {
            if (fscanf(fp, "%lf", &gtemp[k]) != 1) break;
        }
        mstore->sum_snp_effects[i] = gtemp[nc-1];
    }
    fclose(fp);
    free(gtemp);
    
    fp = fopen(config->model_file_path, "r");
    if(fp) {
        char dum[100];
        if (fscanf(fp, "%s %lf", dum, &mstate->mu) != 2) mstate->mu = 0.0;
        fclose(fp);
    }
    return SUCCESS;
}

int load_categories(ModelConfig *config, GenomicData *gdata) {
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;

    if (config->cat) {
        FILE *fp = fopen(config->cat_input_file_path, "r");
        if (!fp) {
            fprintf(stderr, "Cannot open cat file %s\n", config->cat_input_file_path);
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

void write_dgv(ModelConfig *config, GenomicData *gdata) {
    FILE *fp = fopen(config->gv_file_path, "w");
    if(!fp) return;
    
    for(int i=0; i<gdata->num_individuals; i++) {
        if(gdata->trains[i] == 0) {
            fprintf(fp, "%15.7E\n", gdata->predicted_values[i]);
        } else {
            fprintf(fp, "NA\n");
        }
    }
    fclose(fp);
}

void output_model(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore) {
    FILE *fp;
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;
    int ndist = config->num_distributions;

    fp = fopen(config->param_file_path, "w");
    if (!fp) return;
    
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
    
    fp = fopen(config->model_file_path, "w");
    if(!fp) return;
    
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
}

void output_beta(ModelConfig *config, MCMCState *mstate, GenomicData *gdata) {
    if (!config->fp_beta) return;
    for(int i=0; i<gdata->num_loci; i++) {
        fprintf(config->fp_beta, " %15.6E", mstate->snp_effects[i]*mstate->snp_effects[i]);
    }
    fprintf(config->fp_beta, "\n");
}

void cleanup_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    if (config->fp_log) { fclose(config->fp_log); config->fp_log = NULL; }
    if (config->fp_hyp) { fclose(config->fp_hyp); config->fp_hyp = NULL; }
    if (config->fp_snp) { fclose(config->fp_snp); config->fp_snp = NULL; }
    if (config->fp_cat) { fclose(config->fp_cat); config->fp_cat = NULL; }
    if (config->fp_beta) { fclose(config->fp_beta); config->fp_beta = NULL; }
    
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
