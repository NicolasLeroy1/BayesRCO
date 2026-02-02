#include "io.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

void get_size(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    char buffer[1024];
    
    gdata->num_individuals = 0;
    fp = fopen(config->phenotype_file_path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->num_individuals++;
        }
        fclose(fp);
    } else {
        printf("Error opening phen file: %s\n", config->phenotype_file_path);
        exit(1);
    }

    gdata->num_loci = 0;
    fp = fopen(config->bim_file_path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->num_loci++;
        }
        fclose(fp);
    } else {
        printf("Error opening bim file: %s\n", config->bim_file_path);
        exit(1);
    }
}

void load_phenos_plink(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    char line[4096];
    
    gdata->trains = (int*)calloc(gdata->num_individuals, sizeof(int));
    gdata->phenotypes = (double*)calloc(gdata->num_individuals, sizeof(double));
    
    fp = fopen(config->phenotype_file_path, "r");
    if (!fp) exit(1);
    
    for (int i = 0; i < gdata->num_individuals; i++) {
        if (!fgets(line, sizeof(line), fp)) break;
        
        char *ptr = line;
        int col = 0;
        char *token;
        char *saveptr;
        double val = MISSING_VALUE;
        int found = 0;
        
        // Manual tokenization to handle multiple spaces
        char *p = line;
        while (*p && isspace(*p)) p++; // skip leading
        
        int current_col = 0;
        int target_col = 5 + (config->trait_column_index - 1); // 0-based index of target column
        
        // Naive parsing assuming space delimited
        while (*p) {
            if (current_col == target_col) {
                // parse value
                if (strncmp(p, "NA", 2) == 0) {
                    val = MISSING_VALUE;
                } else {
                    char temp[100];
                    int k=0;
                    while (*p && !isspace(*p) && k<99) temp[k++] = *p++;
                    temp[k] = '\0';
                    if (strcmp(temp, "NA") == 0) {
                         val = MISSING_VALUE;
                    } else {
                         val = atof(temp);
                    }
                }
                found = 1;
                break;
            }
            
            // Skip word
            while (*p && !isspace(*p)) p++;
            // Skip space
            while (*p && isspace(*p)) p++;
            current_col++;
        }
        
        if (found && val != MISSING_VALUE) { 
             gdata->trains[i] = 0;
             gdata->phenotypes[i] = val;
        } else {
             gdata->trains[i] = 1;
             gdata->phenotypes[i] = MISSING_VALUE;
        }
    }
    fclose(fp);
}

void load_snp_binary(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    unsigned char b1, b2, b3;
    double igen[4] = {0.0, 3.0, 1.0, 2.0}; // HomRef, Missing, Het, HomAlt (PLINK binary format)
    
    // allocate genotypes in COLUMN-MAJOR order: X[locus * nt + individual]
    // This allows contiguous access when iterating over individuals for a given locus
    gdata->genotypes = (double*)calloc(gdata->num_loci * gdata->num_phenotyped_individuals, sizeof(double));
    
    fp = fopen(config->genotype_file_path, "rb");
    if (!fp) exit(1);
    
    // Header
    if (fread(&b1, 1, 1, fp) != 1) { fclose(fp); exit(1); }
    if (fread(&b2, 1, 1, fp) != 1) { fclose(fp); exit(1); }
    if (fread(&b3, 1, 1, fp) != 1) { fclose(fp); exit(1); }
    
    int nbytes = (gdata->num_individuals + 3) / 4;
    unsigned char *buffer = (unsigned char*)malloc(nbytes);
    
    for (int j = 0; j < gdata->num_loci; j++) {
        if (fread(buffer, 1, nbytes, fp) != nbytes) break;
        
        int tr = 0; // index in training individuals
        for (int i = 0; i < gdata->num_individuals; i++) {
            int byte_idx = i / 4;
            int bit_idx = (i % 4) * 2;
            unsigned char b = buffer[byte_idx];
            unsigned char code = (b >> bit_idx) & 3;
            double val = igen[code];
            
            if (gdata->trains[i] == 0) {
                // Column-major: X[j, tr] = X[j * nt + tr]
                gdata->genotypes[j * gdata->num_phenotyped_individuals + tr] = val;
                tr++;
            }
        }
    }
    
    free(buffer);
    fclose(fp);
}

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

void allocate_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    // Determine nt, ntrain etc.
    if (!config->mcmc) {
        for(int i=0; i<gdata->num_individuals; i++) {
            if (gdata->trains[i] == 0) gdata->trains[i] = 3;
            else if (gdata->trains[i] == 1) gdata->trains[i] = 0;
        }
        for(int i=0; i<gdata->num_individuals; i++) {
            if (gdata->trains[i] == 3) gdata->trains[i] = 1;
        }
    }
    
    gdata->num_phenotyped_individuals = 0;
    for(int i=0; i<gdata->num_individuals; i++) if(gdata->trains[i]==0) gdata->num_phenotyped_individuals++;
    
    // Allocations with NULL checking
    SAFE_CALLOC(gdata->predicted_values, gdata->num_individuals, double);
    
    SAFE_CALLOC(mstate->variance_scaling_factors, config->num_distributions, double);
    SAFE_CALLOC(mstate->genomic_values, config->num_distributions, double);
    SAFE_MALLOC(mstate->p, config->num_distributions, double*);
    SAFE_MALLOC(mstate->log_p, config->num_distributions, double*);
    SAFE_MALLOC(mstate->snps_per_distribution, config->num_distributions, int*);
    SAFE_MALLOC(mstate->variance_per_distribution, config->num_distributions, double*);
    
    for(int i=0; i<config->num_distributions; i++) {
        SAFE_CALLOC(mstate->p[i], config->num_categories, double);
        SAFE_CALLOC(mstate->log_p[i], config->num_categories, double);
        SAFE_CALLOC(mstate->snps_per_distribution[i], config->num_categories, int);
        SAFE_CALLOC(mstate->variance_per_distribution[i], config->num_categories, double);
    }

    SAFE_CALLOC(gdata->permannot, config->num_categories, int);
    SAFE_CALLOC(mstate->dirichlet_priors, config->num_distributions, double);
    SAFE_CALLOC(mstate->dirichlet_scratch, config->num_distributions, double);
    SAFE_CALLOC(mstate->snp_effects, gdata->num_loci, double);
    SAFE_CALLOC(mstate->adjusted_phenotypes, gdata->num_phenotyped_individuals, double);
    SAFE_CALLOC(mstate->z, gdata->num_phenotyped_individuals, double);
    
    SAFE_CALLOC(mstate->log_likelihoods, config->num_distributions, double);
    SAFE_CALLOC(mstate->selection_probs, config->num_distributions, double);
    SAFE_CALLOC(mstate->sstemp, config->num_categories, double);
    
    SAFE_CALLOC(gdata->snp_correlations, gdata->num_loci, double);
    
    SAFE_CALLOC(mstore->sum_snp_effects, gdata->num_loci, double);
    
    // Rank 2 stores
    SAFE_MALLOC(mstore->sum_snps_per_distribution, config->num_distributions, double*);
    SAFE_MALLOC(mstore->sum_variance_per_distribution, config->num_distributions, double*);
    SAFE_MALLOC(mstore->sum_mixture_proportions, config->num_distributions, double*);
    for(int i=0; i<config->num_distributions; i++) {
        SAFE_CALLOC(mstore->sum_snps_per_distribution[i], config->num_categories, double);
        SAFE_CALLOC(mstore->sum_variance_per_distribution[i], config->num_categories, double);
        SAFE_CALLOC(mstore->sum_mixture_proportions[i], config->num_categories, double);
    }
    
    SAFE_MALLOC(mstore->sum_distribution_counts, gdata->num_loci, double*);
    for(int i=0; i<gdata->num_loci; i++) {
        SAFE_CALLOC(mstore->sum_distribution_counts[i], config->num_distributions, double);
    }
    
    SAFE_CALLOC(mstore->mu_vare_store, 4, double);
    SAFE_CALLOC(gdata->allele_frequencies, gdata->num_loci, double);
    SAFE_CALLOC(gdata->permvec, gdata->num_loci, int);
    
    SAFE_MALLOC(gdata->distribution_per_category, gdata->num_loci, int*);
    for(int i=0; i<gdata->num_loci; i++) {
        SAFE_CALLOC(gdata->distribution_per_category[i], config->num_categories, int);
    }
    
    SAFE_CALLOC(mstate->log_distribution_variances, config->num_distributions, double);
    SAFE_CALLOC(mstate->residual_variance_over_distribution_variances, config->num_distributions, double);
    
    SAFE_CALLOC(mstore->varustore, gdata->num_loci, double);
    SAFE_CALLOC(mstore->varistore, gdata->num_loci, double);
    
    SAFE_CALLOC(gdata->current_distribution, gdata->num_loci, int);
    SAFE_CALLOC(gdata->annotations_per_locus, gdata->num_loci, int);
    SAFE_CALLOC(gdata->current_category, gdata->num_loci, int);
    
    SAFE_CALLOC(mstate->ss, config->num_categories, double);
    SAFE_CALLOC(gdata->included_loci, gdata->num_loci, double);
    SAFE_CALLOC(mstate->category_probabilities, config->num_categories, double);
    SAFE_CALLOC(mstate->ytemp, gdata->num_loci, double);
    SAFE_CALLOC(mstate->category_dirichlet_scratch, config->num_categories, double);
    SAFE_CALLOC(gdata->atemp, config->num_categories, int);
    
    SAFE_MALLOC(mstore->sum_category_counts, gdata->num_loci, double*);
    SAFE_MALLOC(gdata->effects_per_category, gdata->num_loci, double*);
    for(int i=0; i<gdata->num_loci; i++) {
        SAFE_CALLOC(mstore->sum_category_counts[i], config->num_categories, double);
        SAFE_CALLOC(gdata->effects_per_category[i], config->num_categories, double);
    }
}

void load_param(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore, MCMCState *mstate) {
    FILE *fp;
    char buffer[1024];
    int nc = config->num_distributions + 1;
    double *gtemp = (double*)calloc(nc, sizeof(double));
    
    fp = fopen(config->param_file_path, "r");
    if(!fp) return;
    if (!fgets(buffer, sizeof(buffer), fp)) { fclose(fp); free(gtemp); return; } // Header
    
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
}

void load_categories(ModelConfig *config, GenomicData *gdata) {
    gdata->categories = (int**)malloc(gdata->num_loci * sizeof(int*));
    for(int i=0; i<gdata->num_loci; i++) {
        gdata->categories[i] = (int*)calloc(config->num_categories, sizeof(int));
    }

    if (config->cat) {
        FILE *fp = fopen(config->cat_input_file_path, "r");
        if (!fp) {
            printf("Cannot open cat file %s\n", config->cat_input_file_path);
            exit(1);
        }
        
        for(int i=0; i<gdata->num_loci; i++) {
            for(int j=0; j<config->num_categories; j++) { 
                int val;
                if (fscanf(fp, "%d", &val) != 1) break;
                gdata->categories[i][j] = val; 
            }
        }
        fclose(fp);
    } else {
        // Default: 1 category, all SNPs in it.
        if (config->num_categories == 1) {
            for(int i=0; i<gdata->num_loci; i++) {
                gdata->categories[i][0] = 1;
            }
        } else {
             printf("Warning: ncat > 1 but no cat file provided. All categories set to 0.\n");
        }
    }
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
    
    fp = fopen(config->param_file_path, "w");
    if (!fp) exit(1);
    
    // Header
    fprintf(fp, "  ");
    for(int i=1; i<=config->num_distributions; i++) fprintf(fp, " PIP%-4d", i);
    fprintf(fp, "  beta");
    for(int i=1; i<=config->num_categories; i++) fprintf(fp, " PAIP%-4d", i);
    fprintf(fp, " Vbeta   Vi\n");
    
    for (int i=0; i<gdata->num_loci; i++) {
        for(int k=0; k<config->num_distributions; k++) fprintf(fp, "%15.7E ", mstore->sum_distribution_counts[i][k]);
        fprintf(fp, "%15.7E ", mstore->sum_snp_effects[i]);
        for(int k=0; k<config->num_categories; k++) fprintf(fp, "%15.7E ", mstore->sum_category_counts[i][k]);
        fprintf(fp, "%15.7E %15.7E\n", mstore->varustore[i], mstore->varistore[i]);
    }
    fclose(fp);
    
    fp = fopen(config->model_file_path, "w");
    if(!fp) exit(1);
    
    fprintf(fp, "Mean      %15.7E\n", mstore->mu_vare_store[0]);
    fprintf(fp, "Nsnp      %15.7E\n", mstore->mu_vare_store[1]);
    fprintf(fp, "Va        %15.7E\n", mstore->mu_vare_store[2]);
    fprintf(fp, "Ve        %15.7E\n", mstore->mu_vare_store[3]);
    
    for(int j=1; j<=config->num_categories; j++) {
        for(int i=1; i<=config->num_distributions; i++) {
             fprintf(fp, "Nk%d_%d      %15.7E\n", i, j, mstore->sum_snps_per_distribution[i-1][j-1]);
        }
    }
    for(int j=1; j<=config->num_categories; j++) {
        for(int i=1; i<=config->num_distributions; i++) {
             fprintf(fp, "Pk%d_%d      %15.7E\n", i, j, mstore->sum_mixture_proportions[i-1][j-1]);
        }
    }
    for(int j=1; j<=config->num_categories; j++) {
        for(int i=1; i<=config->num_distributions; i++) {
             fprintf(fp, "Vk%d_%d      %15.7E\n", i, j, mstore->sum_variance_per_distribution[i-1][j-1]);
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
    int i;
    
    // Close any open file handles
    if (config->fp_log) { fclose(config->fp_log); config->fp_log = NULL; }
    if (config->fp_hyp) { fclose(config->fp_hyp); config->fp_hyp = NULL; }
    if (config->fp_snp) { fclose(config->fp_snp); config->fp_snp = NULL; }
    if (config->fp_cat) { fclose(config->fp_cat); config->fp_cat = NULL; }
    if (config->fp_beta) { fclose(config->fp_beta); config->fp_beta = NULL; }
    
    // Free GenomicData arrays
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
    
    // Free 2D arrays in GenomicData
    if (gdata->categories) {
        for (i = 0; i < gdata->num_loci; i++) {
            SAFE_FREE(gdata->categories[i]);
        }
        SAFE_FREE(gdata->categories);
    }
    
    if (gdata->distribution_per_category) {
        for (i = 0; i < gdata->num_loci; i++) {
            SAFE_FREE(gdata->distribution_per_category[i]);
        }
        SAFE_FREE(gdata->distribution_per_category);
    }
    
    if (gdata->effects_per_category) {
        for (i = 0; i < gdata->num_loci; i++) {
            SAFE_FREE(gdata->effects_per_category[i]);
        }
        SAFE_FREE(gdata->effects_per_category);
    }
    
    // Free MCMCState arrays
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
    
    // Free 2D arrays in MCMCState
    if (mstate->variance_per_distribution) {
        for (i = 0; i < config->num_distributions; i++) {
            SAFE_FREE(mstate->variance_per_distribution[i]);
        }
        SAFE_FREE(mstate->variance_per_distribution);
    }
    
    if (mstate->p) {
        for (i = 0; i < config->num_distributions; i++) {
            SAFE_FREE(mstate->p[i]);
        }
        SAFE_FREE(mstate->p);
    }
    
    if (mstate->log_p) {
        for (i = 0; i < config->num_distributions; i++) {
            SAFE_FREE(mstate->log_p[i]);
        }
        SAFE_FREE(mstate->log_p);
    }
    
    if (mstate->snps_per_distribution) {
        for (i = 0; i < config->num_distributions; i++) {
            SAFE_FREE(mstate->snps_per_distribution[i]);
        }
        SAFE_FREE(mstate->snps_per_distribution);
    }
    
    // Free MCMCStorage arrays
    SAFE_FREE(mstore->sum_snp_effects);
    SAFE_FREE(mstore->mu_vare_store);
    SAFE_FREE(mstore->varustore);
    SAFE_FREE(mstore->varistore);
    
    // Free 2D arrays in MCMCStorage
    if (mstore->sum_snps_per_distribution) {
        for (i = 0; i < config->num_distributions; i++) {
            SAFE_FREE(mstore->sum_snps_per_distribution[i]);
        }
        SAFE_FREE(mstore->sum_snps_per_distribution);
    }
    
    if (mstore->sum_variance_per_distribution) {
        for (i = 0; i < config->num_distributions; i++) {
            SAFE_FREE(mstore->sum_variance_per_distribution[i]);
        }
        SAFE_FREE(mstore->sum_variance_per_distribution);
    }
    
    if (mstore->sum_mixture_proportions) {
        for (i = 0; i < config->num_distributions; i++) {
            SAFE_FREE(mstore->sum_mixture_proportions[i]);
        }
        SAFE_FREE(mstore->sum_mixture_proportions);
    }
    
    if (mstore->sum_distribution_counts) {
        for (i = 0; i < gdata->num_loci; i++) {
            SAFE_FREE(mstore->sum_distribution_counts[i]);
        }
        SAFE_FREE(mstore->sum_distribution_counts);
    }
    
    if (mstore->sum_category_counts) {
        for (i = 0; i < gdata->num_loci; i++) {
            SAFE_FREE(mstore->sum_category_counts[i]);
        }
        SAFE_FREE(mstore->sum_category_counts);
    }
}
