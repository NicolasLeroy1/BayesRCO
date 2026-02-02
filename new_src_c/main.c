/**
 * @file main.c
 * @brief Entry point for the BayesRCO C implementation.
 * 
 * Handles command-line parsing, configuration initialization, data loading,
 * and orchestrates the MCMC run.
 */

#include "bayesRCO.h"
#include "rng.h"
#include "io.h"
#include "utils.h"
#include "mcmc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

/* =========================================================================
 * Additional Constants
 * ========================================================================= */

#define DEFAULT_VARIANCE_GENETIC 0.01
#define DEFAULT_VARIANCE_RESIDUAL 0.01
#define DEFAULT_DIRICHLET_PRIOR 1.0
#define DEFAULT_NUM_CATEGORIES 1
#define DEFAULT_DFVARA -2.0
#define DEFAULT_DFVARE -2.0

/* Default variance scaling factors for 4 distributions */
static const double GPIN_DEFAULTS[DEFAULT_NUM_DISTRIBUTIONS] = {0.0, 0.0001, 0.001, 0.01};

/* =========================================================================
 * Helper Functions
 * ========================================================================= */

/**
 * Check if a file exists and is accessible.
 */
static bool file_exists(const char *filename) {
    return access(filename, F_OK) != -1;
}

/**
 * Safe string copy with bounds checking.
 */
static void safe_strcpy(char *dest, const char *src, size_t dest_size) {
    if (dest_size == 0) return;
    strncpy(dest, src, dest_size - 1);
    dest[dest_size - 1] = '\0';
}

/**
 * Parse a comma-separated list of doubles.
 * Returns the count of parsed values, or -1 on error.
 */
static int parse_double_list(const char *arg, double **out_array) {
    if (!arg || !out_array) return -1;
    
    /* Make a copy since strtok modifies the string */
    char *arg_copy = strdup(arg);
    if (!arg_copy) return -1;
    
    int count = 0;
    double *array = NULL;
    
    char *token = strtok(arg_copy, ",");
    while (token) {
        count++;
        double *tmp = realloc(array, count * sizeof(double));
        if (!tmp) {
            free(array);
            free(arg_copy);
            return -1;
        }
        array = tmp;
        array[count - 1] = atof(token);
        token = strtok(NULL, ",");
    }
    
    free(arg_copy);
    *out_array = array;
    return count;
}

/**
 * Print usage information and exit.
 */
static void print_usage(const char *program_name) {
    fprintf(stderr, "Usage: %s [options]\n\n", program_name);
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -bfile PREFIX    Input file prefix (.bed/.bim/.fam)\n");
    fprintf(stderr, "  -out PREFIX      Output file prefix\n\n");
    fprintf(stderr, "MCMC Options:\n");
    fprintf(stderr, "  -numit N         Number of iterations (default: %d)\n", DEFAULT_NUM_ITERATIONS);
    fprintf(stderr, "  -burnin N        Burn-in iterations (default: %d)\n", DEFAULT_BURNIN);
    fprintf(stderr, "  -thin N          Thinning interval (default: %d)\n", DEFAULT_THINNING);
    fprintf(stderr, "  -ndist N         Number of distributions (default: %d)\n", DEFAULT_NUM_DISTRIBUTIONS);
    fprintf(stderr, "  -seed N          Random seed\n\n");
    fprintf(stderr, "Model Options:\n");
    fprintf(stderr, "  -vara VAL        Prior genetic variance (default: %.2f)\n", DEFAULT_VARIANCE_GENETIC);
    fprintf(stderr, "  -vare VAL        Prior residual variance (default: %.2f)\n", DEFAULT_VARIANCE_RESIDUAL);
    fprintf(stderr, "  -gpin V1,V2,...  Variance scaling factors\n");
    fprintf(stderr, "  -delta V1,...    Dirichlet prior (default: %.1f)\n", DEFAULT_DIRICHLET_PRIOR);
    fprintf(stderr, "  -additive        Use additive model (default: mixture)\n");
    fprintf(stderr, "  -bayesCpi        Enable BayesCpi (estimate pi)\n\n");
    fprintf(stderr, "Annotation Options:\n");
    fprintf(stderr, "  -cat             Enable annotation categories\n");
    fprintf(stderr, "  -ncat N          Number of categories (default: %d)\n", DEFAULT_NUM_CATEGORIES);
    fprintf(stderr, "  -catfile FILE    Category file path\n\n");
    fprintf(stderr, "Output Options:\n");
    fprintf(stderr, "  -snpout          Output SNP-level results\n");
    fprintf(stderr, "  -beta            Output beta values\n");
    fprintf(stderr, "  -permute         Permute marker order\n\n");
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -predict         Prediction mode (instead of MCMC)\n");
    fprintf(stderr, "  -h               Show this help message\n");
}

/**
 * Initialize configuration with default values.
 */
static void init_default_config(ModelConfig *config, MCMCState *mstate) {
    memset(config, 0, sizeof(ModelConfig));
    memset(mstate, 0, sizeof(MCMCState));
    
    /* Algorithm parameters */
    config->trait_column_index = 1;
    config->marker_set_size = 0;
    config->marker_replicates = DEFAULT_MARKER_REPLICATES;
    config->num_iterations = DEFAULT_NUM_ITERATIONS;
    config->burnin_iterations = DEFAULT_BURNIN;
    config->thinning_interval = DEFAULT_THINNING;
    config->num_distributions = DEFAULT_NUM_DISTRIBUTIONS;
    config->num_categories = DEFAULT_NUM_CATEGORIES;
    config->random_seed = 0;
    
    /* Hyperparameters */
    config->dfvara = DEFAULT_DFVARA;
    config->dfvare = DEFAULT_DFVARE;
    mstate->variance_genetic = DEFAULT_VARIANCE_GENETIC;
    mstate->variance_residual = DEFAULT_VARIANCE_RESIDUAL;
    
    /* Flags */
    config->mcmc = true;
    config->snpout = false;
    config->permute = false;
    config->cat = false;
    config->beta = false;
    config->mixture = true;
    config->nobayesCpi = true;
    config->VCE = false;
}

/**
 * Derive output file paths from prefix.
 */
static void derive_file_paths(ModelConfig *config) {
    snprintf(config->genotype_file_path, PATH_MAX_LENGTH, "%s.bed", config->input_prefix);
    snprintf(config->phenotype_file_path, PATH_MAX_LENGTH, "%s.fam", config->input_prefix);
    snprintf(config->bim_file_path, PATH_MAX_LENGTH, "%s.bim", config->input_prefix);
    
    if (strlen(config->output_prefix) > 0) {
        snprintf(config->log_file_path, PATH_MAX_LENGTH, "%s.log", config->output_prefix);
        snprintf(config->freq_file_path, PATH_MAX_LENGTH, "%s.frq", config->output_prefix);
        snprintf(config->gv_file_path, PATH_MAX_LENGTH, "%s.gv", config->output_prefix);
        snprintf(config->hyp_file_path, PATH_MAX_LENGTH, "%s.hyp", config->output_prefix);
        snprintf(config->snp_file_path, PATH_MAX_LENGTH, "%s.snp", config->output_prefix);
        snprintf(config->cat_file_path, PATH_MAX_LENGTH, "%s.catit", config->output_prefix);
        snprintf(config->beta_file_path, PATH_MAX_LENGTH, "%s.beta", config->output_prefix);
        snprintf(config->model_file_path, PATH_MAX_LENGTH, "%s.model", config->output_prefix);
        snprintf(config->param_file_path, PATH_MAX_LENGTH, "%s.param", config->output_prefix);
    }
}

/**
 * Initialize variance component estimation parameters.
 */
static void init_variance_components(ModelConfig *config, GenomicData *gdata, MCMCState *mstate) {
    if (config->dfvara < DEFAULT_DFVARA) {
        config->VCE = false;
        
        /* Compute mean phenotype for training individuals */
        double sum_y = 0.0;
        int count = 0;
        for (int i = 0; i < gdata->num_individuals; i++) {
            if (gdata->trains[i] == 0) {
                sum_y += gdata->phenotypes[i];
                count++;
            }
        }
        mstate->yhat = sum_y / mstate->nnind;
        
        /* Compute phenotype variance */
        double sum_sq = 0.0;
        for (int i = 0; i < gdata->num_individuals; i++) {
            if (gdata->trains[i] == 0) {
                double diff = gdata->phenotypes[i] - mstate->yhat;
                sum_sq += diff * diff;
            }
        }
        mstate->vary = sum_sq / (mstate->nnind - 1.0);
        mstate->variance_genetic = mstate->variance_genetic * mstate->vary;
    } else {
        config->VCE = true;
        config->vara_ap = mstate->variance_genetic;
        config->vare_ap = mstate->variance_residual;
        if (config->dfvara == DEFAULT_DFVARA) config->vara_ap = 0.0;
        if (config->dfvare == DEFAULT_DFVARE) config->vare_ap = 0.0;
    }
}

/**
 * Open log file and write header.
 */
static bool init_logging(ModelConfig *config) {
    config->fp_log = fopen(config->log_file_path, "w");
    if (!config->fp_log) {
        fprintf(stderr, "Error: Cannot open log file: %s\n", config->log_file_path);
        return false;
    }
    
    fprintf(config->fp_log, "Program BayesRCO C Port\n");
    
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    char datetime[64];
    strftime(datetime, sizeof(datetime), "%Y-%m-%d %H:%M:%S", t);
    fprintf(config->fp_log, "Run started at %s\n", datetime);
    fprintf(config->fp_log, "Prefix for input files           : %s\n", config->input_prefix);
    fprintf(config->fp_log, "Prefix for output files          : %s\n", config->output_prefix);
    
    return true;
}

/* =========================================================================
 * Main Entry Point
 * ========================================================================= */

int main(int argc, char **argv) {
    ModelConfig config;
    GenomicData gdata;
    MCMCState mstate;
    MCMCStorage mstore;
    prng_state rs;
    
    /* Initialize with defaults */
    init_default_config(&config, &mstate);
    memset(&gdata, 0, sizeof(GenomicData));
    memset(&mstore, 0, sizeof(MCMCStorage));
    
    /* Temporary storage for array arguments */
    int gpin_count = 0;
    double *gpin_args = NULL;
    int delta_count = 0;
    double *delta_args = NULL;
    
    /* Print usage if no arguments */
    if (argc == 1) {
        print_usage(argv[0]);
        return 0;
    }
    
    /* Parse command-line arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-bfile") == 0) {
            if (i + 1 < argc) safe_strcpy(config.input_prefix, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-out") == 0) {
            if (i + 1 < argc) safe_strcpy(config.output_prefix, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-n") == 0) {
            if (i + 1 < argc) config.trait_column_index = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-vara") == 0) {
            if (i + 1 < argc) mstate.variance_genetic = atof(argv[++i]);
        } else if (strcmp(argv[i], "-vare") == 0) {
            if (i + 1 < argc) mstate.variance_residual = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dfvara") == 0) {
            if (i + 1 < argc) config.dfvara = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dfvare") == 0) {
            if (i + 1 < argc) config.dfvare = atof(argv[++i]);
        } else if (strcmp(argv[i], "-delta") == 0) {
            if (i + 1 < argc) {
                delta_count = parse_double_list(argv[++i], &delta_args);
                if (delta_count < 0) {
                    fprintf(stderr, "Error: Failed to parse -delta argument\n");
                    return 1;
                }
            }
        } else if (strcmp(argv[i], "-msize") == 0) {
            if (i + 1 < argc) config.marker_set_size = atoi(argv[++i]);
            config.permute = true;
        } else if (strcmp(argv[i], "-mrep") == 0) {
            if (i + 1 < argc) config.marker_replicates = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-numit") == 0) {
            if (i + 1 < argc) config.num_iterations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-burnin") == 0) {
            if (i + 1 < argc) config.burnin_iterations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-thin") == 0) {
            if (i + 1 < argc) config.thinning_interval = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-ndist") == 0) {
            if (i + 1 < argc) config.num_distributions = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-gpin") == 0) {
            if (i + 1 < argc) {
                gpin_count = parse_double_list(argv[++i], &gpin_args);
                if (gpin_count < 0) {
                    fprintf(stderr, "Error: Failed to parse -gpin argument\n");
                    free(delta_args);
                    return 1;
                }
            }
        } else if (strcmp(argv[i], "-seed") == 0) {
            if (i + 1 < argc) config.random_seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-predict") == 0) {
            config.mcmc = false;
        } else if (strcmp(argv[i], "-snpout") == 0) {
            config.snpout = true;
        } else if (strcmp(argv[i], "-permute") == 0) {
            config.permute = true;
        } else if (strcmp(argv[i], "-model") == 0) {
            if (i + 1 < argc) safe_strcpy(config.model_file_path, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-freq") == 0) {
            if (i + 1 < argc) safe_strcpy(config.freq_file_path, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-param") == 0) {
            if (i + 1 < argc) safe_strcpy(config.param_file_path, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-cat") == 0) {
            config.cat = true;
        } else if (strcmp(argv[i], "-beta") == 0) {
            config.beta = true;
        } else if (strcmp(argv[i], "-ncat") == 0) {
            if (i + 1 < argc) config.num_categories = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-catfile") == 0) {
            if (i + 1 < argc) safe_strcpy(config.cat_input_file_path, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-additive") == 0) {
            config.mixture = false;
        } else if (strcmp(argv[i], "-bayesCpi") == 0) {
            config.nobayesCpi = false;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            free(gpin_args);
            free(delta_args);
            return 0;
        } else {
            fprintf(stderr, "Warning: Unknown option '%s'\n", argv[i]);
        }
    }
    
    /* Validate required inputs */
    if (strlen(config.input_prefix) == 0 && config.mcmc) {
        fprintf(stderr, "Error: -bfile is required for MCMC run\n");
        free(gpin_args);
        free(delta_args);
        return 1;
    }
    
    /* Derive file paths */
    derive_file_paths(&config);
    
    /* Check input files exist */
    if (config.mcmc) {
        if (!file_exists(config.genotype_file_path)) {
            fprintf(stderr, "Error: File not found: %s\n", config.genotype_file_path);
            free(gpin_args);
            free(delta_args);
            return 1;
        }
    }
    
    /* Initialize logging */
    if (!init_logging(&config)) {
        free(gpin_args);
        free(delta_args);
        return 1;
    }
    
    /* Load data */
    get_size(&config, &gdata);
    load_phenos_plink(&config, &gdata);
    allocate_data(&config, &gdata, &mstate, &mstore);
    
    /* Set variance scaling factors */
    if (gpin_count > 0) {
        if (gpin_count != config.num_distributions) {
            fprintf(stderr, "Error: -gpin count (%d) must match -ndist (%d)\n", 
                    gpin_count, config.num_distributions);
            cleanup_data(&config, &gdata, &mstate, &mstore);
            free(gpin_args);
            free(delta_args);
            return 1;
        }
        for (int i = 0; i < config.num_distributions; i++) {
            mstate.variance_scaling_factors[i] = gpin_args[i];
        }
    } else {
        /* Use defaults, extending with zeros if needed */
        for (int i = 0; i < config.num_distributions; i++) {
            if (i < DEFAULT_NUM_DISTRIBUTIONS) {
                mstate.variance_scaling_factors[i] = GPIN_DEFAULTS[i];
            } else {
                mstate.variance_scaling_factors[i] = 0.0;
            }
        }
    }
    free(gpin_args);
    gpin_args = NULL;
    
    /* Set Dirichlet priors */
    if (delta_count > 0) {
        if (delta_count == 1) {
            for (int i = 0; i < config.num_distributions; i++) {
                mstate.dirichlet_priors[i] = delta_args[0];
            }
        } else if (delta_count == config.num_distributions) {
            for (int i = 0; i < config.num_distributions; i++) {
                mstate.dirichlet_priors[i] = delta_args[i];
            }
        } else {
            fprintf(stderr, "Error: -delta count must be 1 or match -ndist\n");
            cleanup_data(&config, &gdata, &mstate, &mstore);
            free(delta_args);
            return 1;
        }
    } else {
        for (int i = 0; i < config.num_distributions; i++) {
            mstate.dirichlet_priors[i] = DEFAULT_DIRICHLET_PRIOR;
        }
    }
    free(delta_args);
    delta_args = NULL;
    
    /* Load annotation categories */
    load_categories(&config, &gdata);
    
    /* Log data summary */
    if (config.mcmc && config.fp_log) {
        fprintf(config.fp_log, "Phenotype column               = %8d\n", config.trait_column_index);
        fprintf(config.fp_log, "No. of loci                    = %8d\n", gdata.num_loci);
        fprintf(config.fp_log, "No. of individuals             = %8d\n", gdata.num_individuals);
        fprintf(config.fp_log, "No. of training individuals    = %8d\n", gdata.num_phenotyped_individuals);
        fflush(config.fp_log);
    }
    
    /* Load genotypes and initialize */
    load_snp_binary(&config, &gdata);
    xcenter(&config, &gdata);
    init_random_seed_custom(&config, &rs);
    
    if (config.mcmc) {
        mstate.nnind = (double)gdata.num_phenotyped_individuals;
        
        /* Open output files */
        if (config.snpout) {
            config.fp_snp = fopen(config.snp_file_path, "w");
            if (!config.fp_snp) {
                fprintf(stderr, "Warning: Cannot open SNP output file: %s\n", config.snp_file_path);
            }
        }
        if (config.cat) {
            config.fp_cat = fopen(config.cat_file_path, "w");
            if (!config.fp_cat) {
                fprintf(stderr, "Warning: Cannot open category output file: %s\n", config.cat_file_path);
            }
        }
        if (config.beta) {
            config.fp_beta = fopen(config.beta_file_path, "w");
            if (!config.fp_beta) {
                fprintf(stderr, "Warning: Cannot open beta output file: %s\n", config.beta_file_path);
            }
        }
        
        config.fp_hyp = fopen(config.hyp_file_path, "w");
        if (config.fp_hyp) {
            fprintf(config.fp_hyp, " Replicate       Nsnp              Va              Ve ");
            fprintf(config.fp_hyp, "\n");
        } else {
            fprintf(stderr, "Warning: Cannot open hyperparameter file: %s\n", config.hyp_file_path);
        }
        
        /* Initialize variance components */
        init_variance_components(&config, &gdata, &mstate);
        
        /* Run MCMC */
        run_mcmc(&config, &gdata, &mstate, &mstore, &rs);
        
    } else {
        /* Prediction mode */
        fprintf(stderr, "Prediction mode not fully implemented in C port yet.\n");
    }
    
    /* Cleanup all allocated memory and close file handles */
    cleanup_data(&config, &gdata, &mstate, &mstore);
    
    return 0;
}
