/**
 * @file main.c
 * @brief CLI entry point for BayesRCO.
 * 
 * Uses decoupled libraries:
 * - libbayesrco_stats: Pure statistical computation
 * - libbayesrco_io: I/O and data management
 */

#include "../libbayesrco_stats/bayesrco_stats.h"
#include "../libbayesrco_io/bayesrco_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#define DEFAULT_VARIANCE_GENETIC 0.01
#define DEFAULT_VARIANCE_RESIDUAL 0.01
#define DEFAULT_DIRICHLET_PRIOR 1.0
#define DEFAULT_NUM_CATEGORIES 1
#define DEFAULT_DFVARA -2.0
#define DEFAULT_DFVARE -2.0

static const double GPIN_DEFAULTS[DEFAULT_NUM_DISTRIBUTIONS] = {0.0, 0.0001, 0.001, 0.01};

static bool file_exists(const char *filename) {
    return access(filename, F_OK) != -1;
}

static void safe_strcpy(char *dest, const char *src, size_t dest_size) {
    if (dest_size == 0) return;
    strncpy(dest, src, dest_size - 1);
    dest[dest_size - 1] = '\0';
}

static int parse_double_list(const char *arg, double **out_array) {
    if (!arg || !out_array) return -1;
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

static void init_default_config(IOConfig *config) {
    memset(config, 0, sizeof(IOConfig));
    
    /* IO Flags & Options */
    config->mcmc = true;
    config->cat = false;
    config->beta = false;
    config->trait_column_index = 1;
    config->random_seed = 0;
    
    /* Stats Model Config */
    config->model.marker_set_size = 0;
    config->model.marker_replicates = DEFAULT_MARKER_REPLICATES;
    config->model.num_iterations = DEFAULT_NUM_ITERATIONS;
    config->model.burnin_iterations = DEFAULT_BURNIN;
    config->model.thinning_interval = DEFAULT_THINNING;
    config->model.num_distributions = DEFAULT_NUM_DISTRIBUTIONS;
    config->model.num_categories = DEFAULT_NUM_CATEGORIES;
    config->model.dfvara = DEFAULT_DFVARA;
    config->model.dfvare = DEFAULT_DFVARE;
    config->model.vara_ap = DEFAULT_VARIANCE_GENETIC; 
    config->model.vare_ap = DEFAULT_VARIANCE_RESIDUAL;
    config->model.mixture = true;
    config->model.nobayesCpi = true;
    config->model.VCE = false;
}

static void derive_file_paths(IOConfig *config, const char *input_prefix, const char *output_prefix) {
    if (strlen(input_prefix) > 0) {
        snprintf(config->geno_file_path, PATH_MAX_LENGTH, "%s.bed", input_prefix);
        snprintf(config->pheno_file_path, PATH_MAX_LENGTH, "%s.fam", input_prefix);
    }
    
    if (strlen(output_prefix) > 0) {
        snprintf(config->log_file_path, PATH_MAX_LENGTH, "%s.log", output_prefix);
        snprintf(config->freq_file_path, PATH_MAX_LENGTH, "%s.frq", output_prefix);
        snprintf(config->gv_file_path, PATH_MAX_LENGTH, "%s.gv", output_prefix);
        snprintf(config->hyp_output_file_path, PATH_MAX_LENGTH, "%s.hyp", output_prefix);
        snprintf(config->snp_output_file_path, PATH_MAX_LENGTH, "%s.snp", output_prefix);
        snprintf(config->cat_output_file_path, PATH_MAX_LENGTH, "%s.catit", output_prefix);
        snprintf(config->beta_file_path, PATH_MAX_LENGTH, "%s.beta", output_prefix);
        snprintf(config->model_file_path, PATH_MAX_LENGTH, "%s.model", output_prefix);
        snprintf(config->param_file_path, PATH_MAX_LENGTH, "%s.param", output_prefix);
    }
}

static bool init_logging(IOConfig *config, const char *input_prefix, const char *output_prefix) {
    config->fp_log = fopen(config->log_file_path, "w");
    if (!config->fp_log) {
        fprintf(stderr, "Error: Cannot open log file: %s\n", config->log_file_path);
        return false;
    }
    fprintf(config->fp_log, "Program BayesRCO C Port (Decoupled)\n");
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    char datetime[64];
    strftime(datetime, sizeof(datetime), "%Y-%m-%d %H:%M:%S", t);
    fprintf(config->fp_log, "Run started at %s\n", datetime);
    fprintf(config->fp_log, "Prefix for input files           : %s\n", input_prefix);
    fprintf(config->fp_log, "Prefix for output files          : %s\n", output_prefix);
    return true;
}

int main(int argc, char **argv) {
    IOConfig config;
    GenomicData gdata;
    
    char input_prefix[PATH_MAX_LENGTH] = "";
    char output_prefix[PATH_MAX_LENGTH] = "";
    
    init_default_config(&config);
    memset(&gdata, 0, sizeof(GenomicData));
    
    int gpin_count = 0;
    double *gpin_args = NULL;
    int delta_count = 0;
    double *delta_args = NULL;
    double init_vara = DEFAULT_VARIANCE_GENETIC;
    double init_vare = DEFAULT_VARIANCE_RESIDUAL;
    
    if (argc == 1) {
        print_usage(argv[0]);
        return 0;
    }
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-bfile") == 0) {
            if (i + 1 < argc) safe_strcpy(input_prefix, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-out") == 0) {
            if (i + 1 < argc) safe_strcpy(output_prefix, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-n") == 0) {
            if (i + 1 < argc) config.trait_column_index = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-vara") == 0) {
            if (i + 1 < argc) init_vara = atof(argv[++i]);
        } else if (strcmp(argv[i], "-vare") == 0) {
            if (i + 1 < argc) init_vare = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dfvara") == 0) {
            if (i + 1 < argc) config.model.dfvara = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dfvare") == 0) {
            if (i + 1 < argc) config.model.dfvare = atof(argv[++i]);
        } else if (strcmp(argv[i], "-delta") == 0) {
            if (i + 1 < argc) {
                delta_count = parse_double_list(argv[++i], &delta_args);
                if (delta_count < 0) return 1;
            }
        } else if (strcmp(argv[i], "-msize") == 0) {
            if (i + 1 < argc) config.model.marker_set_size = atoi(argv[++i]);
            config.model.permute = true;
        } else if (strcmp(argv[i], "-mrep") == 0) {
            if (i + 1 < argc) config.model.marker_replicates = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-numit") == 0) {
            if (i + 1 < argc) config.model.num_iterations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-burnin") == 0) {
            if (i + 1 < argc) config.model.burnin_iterations = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-thin") == 0) {
            if (i + 1 < argc) config.model.thinning_interval = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-ndist") == 0) {
            if (i + 1 < argc) config.model.num_distributions = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-gpin") == 0) {
            if (i + 1 < argc) {
                gpin_count = parse_double_list(argv[++i], &gpin_args);
                if (gpin_count < 0) {
                    free(delta_args);
                    return 1;
                }
            }
        } else if (strcmp(argv[i], "-seed") == 0) {
             if (i + 1 < argc) config.random_seed = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "-predict") == 0) {
            config.mcmc = false;
        } else if (strcmp(argv[i], "-snpout") == 0) {
            /* Handled by IO config flags implicit logic for now */
        } else if (strcmp(argv[i], "-permute") == 0) {
            config.model.permute = true;
        } else if (strcmp(argv[i], "-cat") == 0) {
            config.cat = true;
        } else if (strcmp(argv[i], "-beta") == 0) {
            config.beta = true;
        } else if (strcmp(argv[i], "-ncat") == 0) {
            if (i + 1 < argc) config.model.num_categories = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-catfile") == 0) {
            if (i + 1 < argc) safe_strcpy(config.cat_input_file_path, argv[++i], PATH_MAX_LENGTH);
        } else if (strcmp(argv[i], "-additive") == 0) {
            config.model.mixture = false;
        } else if (strcmp(argv[i], "-bayesCpi") == 0) {
            config.model.nobayesCpi = false;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            free(gpin_args); free(delta_args);
            return 0;
        }
    }
    
    if (strlen(input_prefix) == 0 && config.mcmc) return 1;
    derive_file_paths(&config, input_prefix, output_prefix);
    
    if (config.mcmc && !file_exists(config.geno_file_path)) {
        fprintf(stderr, "Genotype file not found: %s\n", config.geno_file_path);
        return 1;
    }
    
    if (!init_logging(&config, input_prefix, output_prefix)) return 1;
    
    fprintf(config.fp_log, "Reading PLINK files...\n");
    if (plink_get_size(input_prefix, &gdata.num_individuals, &gdata.num_loci) != SUCCESS) {
        fprintf(stderr, "Error: Failed to determine dataset dimensions\n");
        return 1;
    }
    
    if (io_load_phenotypes(&config, &gdata) != SUCCESS) {
        fprintf(stderr, "Error: Failed to load phenotype data\n");
        return 1;
    }
    
    gdata.num_phenotyped_individuals = 0;
    for(int i=0; i<gdata.num_individuals; i++) {
        if(gdata.trains[i] == 0) gdata.num_phenotyped_individuals++;
    }
    
    if (io_load_genotypes(&config, &gdata) != SUCCESS) {
        fprintf(stderr, "Error: Failed to load genotype data\n");
        return 1;
    }
    
    gdata.categories = (int*)calloc(gdata.num_loci * config.model.num_categories, sizeof(int));
    if (io_load_categories(&config, &gdata) != SUCCESS) {
         fprintf(stderr, "Error: Failed to load category file\n");
         return 1;
    }
    
    gdata.allele_frequencies = (double*)calloc(gdata.num_loci, sizeof(double));
    extern void standardize_genotypes(GenomicData *gdata, bool compute_freqs);
    
    if (config.mcmc) {
        standardize_genotypes(&gdata, true);
        io_write_frequencies(&config, &gdata);
    } else {
        io_load_frequencies(&config, &gdata);
        standardize_genotypes(&gdata, false);
    }
    
    MCMCParams params;
    memset(&params, 0, sizeof(MCMCParams));
    params.num_iterations = config.model.num_iterations;
    params.burnin_iterations = config.model.burnin_iterations;
    params.thinning_interval = config.model.thinning_interval;
    params.num_distributions = config.model.num_distributions;
    params.num_categories = config.model.num_categories;
    params.permute = config.model.permute;
    params.VCE = config.model.VCE;
    params.dfvara = config.model.dfvara;
    params.dfvare = config.model.dfvare;
    params.vara_ap = config.model.vara_ap;
    params.vare_ap = config.model.vare_ap;
    params.variance_genetic = init_vara;
    params.variance_residual = init_vare;
    
    params.variance_scaling_factors = (double*)calloc(params.num_distributions, sizeof(double));
    params.dirichlet_priors = (double*)calloc(params.num_distributions, sizeof(double));
    
    if (gpin_count > 0) {
        for(int i=0; i<params.num_distributions; i++) params.variance_scaling_factors[i] = gpin_args[i];
    } else {
         for(int i=0; i<params.num_distributions; i++) {
            if (i < DEFAULT_NUM_DISTRIBUTIONS) params.variance_scaling_factors[i] = GPIN_DEFAULTS[i];
            else params.variance_scaling_factors[i] = 0.0;
         }
    }
    free(gpin_args);
    
    if (delta_count > 0) {
        for(int i=0; i<params.num_distributions; i++) 
            params.dirichlet_priors[i] = (delta_count == 1) ? delta_args[0] : delta_args[i];
    } else {
        for(int i=0; i<params.num_distributions; i++) params.dirichlet_priors[i] = DEFAULT_DIRICHLET_PRIOR;
    }
    free(delta_args);
    
    MCMCResults results;
    memset(&results, 0, sizeof(MCMCResults));

    if (allocate_results(&results, gdata.num_loci, params.num_distributions, params.num_categories, gdata.num_individuals) != SUCCESS) {
        fprintf(stderr, "Error: Failed to allocate results structure\n");
        return 1;
    }
    
    fprintf(config.fp_log, "Starting MCMC...\n");
    
    int ret = SUCCESS;
    if (config.model.mixture) {
        ret = run_bayesrco_mixture(
            &params, 
            gdata.genotypes, 
            gdata.phenotypes, 
            gdata.categories, 
            gdata.num_loci, 
            gdata.num_phenotyped_individuals, 
            config.model.num_categories,
            gdata.num_individuals, 
            gdata.trains,
            config.random_seed,
            &results
        );
    } else {
         if (!config.model.nobayesCpi) {
             ret = run_bayesrco_bayesCpi(
                &params, gdata.genotypes, gdata.phenotypes, gdata.categories,
                gdata.num_loci, gdata.num_phenotyped_individuals, config.model.num_categories,
                gdata.num_individuals, gdata.trains, config.random_seed, &results
             );
         } else {
             ret = run_bayesrco_additive(
                &params, gdata.genotypes, gdata.phenotypes, gdata.categories,
                gdata.num_loci, gdata.num_phenotyped_individuals, config.model.num_categories,
                gdata.num_individuals, gdata.trains, config.random_seed, &results
             );
         }
    }
    
    if (ret == SUCCESS) {
        fprintf(config.fp_log, "MCMC completed successfully.\n");
        io_write_results(&config, &results, gdata.num_loci, gdata.num_individuals);
        io_write_predictions(&config, results.predicted_values, gdata.num_individuals, gdata.trains);
    } else {
         fprintf(stderr, "MCMC execution failed with code %d\n", ret);
    }
    
    io_cleanup(&config, &gdata, NULL, NULL);
    free(params.variance_scaling_factors);
    free(params.dirichlet_priors);
    free_results(&results);
    
    return 0;
}
