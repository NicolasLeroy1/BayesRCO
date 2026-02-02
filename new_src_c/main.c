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

// Helper to check file existence
bool file_exists(const char *filename) {
    if (access(filename, F_OK) != -1) return true;
    return false;
}


int main(int argc, char **argv) {
    // 1. Initialize Default Configuration
    ModelConfig config;
    GenomicData gdata;
    MCMCState mstate;
    MCMCStorage mstore;
    prng_state rs;

    memset(&config, 0, sizeof(ModelConfig));
    memset(&gdata, 0, sizeof(GenomicData));
    memset(&mstate, 0, sizeof(MCMCState));
    memset(&mstore, 0, sizeof(MCMCStorage));

    // Defaults from mod_cmd.f90
    config.trait_pos = 1;
    mstate.vara = 0.01;
    mstate.vare = 0.01;
    config.dfvara = -2.0;
    config.dfvare = -2.0;
    // delta default 1.0, but handled as array usually
    // mod_data says real(dp) :: delta (allocatable)
    // We will allocate delta later, but need initial value.
    double delta_default = 1.0;
    
    config.msize = 0;
    config.mrep = 5000;
    config.numit = 50000;
    config.burnin = 20000;
    config.thin = 10;
    config.ndist = 4;
    // gpin default: 0.0, 0.0001, 0.001, 0.01 (Wait, 4 values)
    double gpin_defaults[] = {0.0, 0.0001, 0.001, 0.01};
    int gpin_defaults_len = 4;
    
    config.seed1 = 0;
    config.mcmc = true; // -predict defaults false -> mcmc true
    config.snpout = false;
    config.permute = false;
    config.cat = false;
    config.beta = false;
    config.ncat = 1;
    config.mixture = true; // -additive defaults false -> mixture true
    config.nobayesCpi = true; // -bayesCpi defaults false -> nobayesCpi true
    
    // Command line parsing
    if (argc == 1) {
        printf("Usage: bayesRCO_c [options]\n");
        exit(0);
    }
    
    // We need temporary storage for array args
    int gpin_count = 0;
    double *gpin_args = NULL;
    int delta_count = 0;
    double *delta_args = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-bfile") == 0) {
            if (i+1 < argc) strcpy(config.inprefix, argv[++i]);
        } else if (strcmp(argv[i], "-out") == 0) {
            if (i+1 < argc) strcpy(config.outprefix, argv[++i]);
        } else if (strcmp(argv[i], "-n") == 0) {
            if (i+1 < argc) config.trait_pos = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-vara") == 0) {
            if (i+1 < argc) mstate.vara = atof(argv[++i]);
        } else if (strcmp(argv[i], "-vare") == 0) {
            if (i+1 < argc) mstate.vare = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dfvara") == 0) {
            if (i+1 < argc) config.dfvara = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dfvare") == 0) {
            if (i+1 < argc) config.dfvare = atof(argv[++i]);
        } else if (strcmp(argv[i], "-delta") == 0) {
            if (i+1 < argc) {
                // Parse comma separated list
                char *token = strtok(argv[++i], ",");
                while (token) {
                    delta_count++;
                    delta_args = (double*)realloc(delta_args, delta_count * sizeof(double));
                    delta_args[delta_count-1] = atof(token);
                    token = strtok(NULL, ",");
                }
            }
        } else if (strcmp(argv[i], "-msize") == 0) {
            if (i+1 < argc) config.msize = atoi(argv[++i]);
            config.permute = true;
        } else if (strcmp(argv[i], "-mrep") == 0) {
            if (i+1 < argc) config.mrep = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-numit") == 0) {
            if (i+1 < argc) config.numit = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-burnin") == 0) {
            if (i+1 < argc) config.burnin = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-thin") == 0) {
            if (i+1 < argc) config.thin = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-ndist") == 0) {
            if (i+1 < argc) config.ndist = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-gpin") == 0) {
            if (i+1 < argc) {
                char *token = strtok(argv[++i], ",");
                while (token) {
                    gpin_count++;
                    gpin_args = (double*)realloc(gpin_args, gpin_count * sizeof(double));
                    gpin_args[gpin_count-1] = atof(token);
                    token = strtok(NULL, ",");
                }
            }
        } else if (strcmp(argv[i], "-seed") == 0) {
            if (i+1 < argc) config.seed1 = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-predict") == 0) {
            config.mcmc = false;
        } else if (strcmp(argv[i], "-snpout") == 0) {
            config.snpout = true;
        } else if (strcmp(argv[i], "-permute") == 0) {
            config.permute = true;
        } else if (strcmp(argv[i], "-model") == 0) {
             if (i+1 < argc) strcpy(config.modfil, argv[++i]);
        } else if (strcmp(argv[i], "-freq") == 0) {
             if (i+1 < argc) strcpy(config.freqfil, argv[++i]);
        } else if (strcmp(argv[i], "-param") == 0) {
             if (i+1 < argc) strcpy(config.paramfil, argv[++i]);
        } else if (strcmp(argv[i], "-cat") == 0) {
            config.cat = true;
        } else if (strcmp(argv[i], "-beta") == 0) {
            config.beta = true;
        } else if (strcmp(argv[i], "-ncat") == 0) {
             if (i+1 < argc) config.ncat = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-catfile") == 0) {
             if (i+1 < argc) strcpy(config.catRC, argv[++i]);
        } else if (strcmp(argv[i], "-additive") == 0) {
            config.mixture = false;
        } else if (strcmp(argv[i], "-bayesCpi") == 0) {
            config.nobayesCpi = false;
        } else if (strcmp(argv[i], "-h") == 0) {
            printf("Usage: bayesRCO_c [options] ... (see mod_cmd.f90 for full list)\n");
            exit(0);
        }
    }
    
    // Validate inputs
    if (strlen(config.inprefix) == 0 && config.mcmc) {
        printf("Error: -bfile is required for MCMC run\n");
        exit(1);
    }
    
    // Derived file paths
    sprintf(config.genfil, "%s.bed", config.inprefix);
    sprintf(config.phenfil, "%s.fam", config.inprefix);
    sprintf(config.bimfil, "%s.bim", config.inprefix);
    if (strlen(config.outprefix) > 0) {
        sprintf(config.logfil, "%s.log", config.outprefix);
        sprintf(config.freqfil, "%s.frq", config.outprefix);
        sprintf(config.mbvfil, "%s.gv", config.outprefix);
        sprintf(config.hypfil, "%s.hyp", config.outprefix);
        sprintf(config.locfil, "%s.snp", config.outprefix);
        sprintf(config.catfil, "%s.catit", config.outprefix);
        sprintf(config.betafil, "%s.beta", config.outprefix);
        sprintf(config.modfil, "%s.model", config.outprefix);
        sprintf(config.paramfil, "%s.param", config.outprefix);
    }

    if (config.mcmc) {
        if (!file_exists(config.genfil)) { printf("File not found: %s\n", config.genfil); exit(1); }
    } else {
        // Checking predict files?
    }
    
    // Init Logging
    if (config.mcmc) {
        config.fp_log = fopen(config.logfil, "w");
        if (config.fp_log) {
            fprintf(config.fp_log, "Program BayesRCO C Port\n");
            // Add date/time logic
            time_t now = time(NULL);
            struct tm *t = localtime(&now);
            char cdate[20], ctime_str[20];
            strftime(cdate, sizeof(cdate), "%Y%m%d", t);
            strftime(ctime_str, sizeof(ctime_str), "%H%M%S", t);
            fprintf(config.fp_log, "Run started at %s %s\n", cdate, ctime_str);
            fprintf(config.fp_log, "Prefix for input files           : %s\n", config.inprefix);
            fprintf(config.fp_log, "Prefix for output files          : %s\n", config.outprefix);
            // ... more logging matching Fortran
        }
    } else {
        // Predict logging
        config.fp_log = fopen(config.logfil, "w"); 
        fprintf(config.fp_log, "Program BayesR C Port\n");
    }
    
    // Load Data Phase
    get_size(&config, &gdata);
    load_phenos_plink(&config, &gdata);
    
    // Allocations handled in allocate_data for arrays dependent on gpin?
    // Note: gpin needs to be allocated in mstate before usage? 
    // allocate_data allocates mstate->gpin.
    // So we should temporarily store gpin args, and fill after allocation.
    
    allocate_data(&config, &gdata, &mstate, &mstore);
    
    // Parse Priors Logic (filling allocated arrays)
    if (gpin_count > 0) {
        if (gpin_count != config.ndist) {
             printf("Error: -gpin count mismatch with -ndist\n");
             exit(1);
        }
        for(int i=0; i<config.ndist; i++) mstate.gpin[i] = gpin_args[i];
    } else {
        if (config.ndist != 4) {
             // If ndist changed but no gpin provided, we might have issue if defaults assume 4
             // For now use defaults up to 4 or 0 if more?
             for(int i=0; i<config.ndist; i++) {
                 if(i < 4) mstate.gpin[i] = gpin_defaults[i];
                 else mstate.gpin[i] = 0.0;
             }
        } else {
             for(int i=0; i<4; i++) mstate.gpin[i] = gpin_defaults[i];
        }
    }
    if (gpin_args) free(gpin_args);
    
    // Delta
    if (delta_count > 0) {
        if (delta_count == 1) {
            for(int i=0; i<config.ndist; i++) mstate.delta[i] = delta_args[0];
        } else if (delta_count == config.ndist) {
            for(int i=0; i<config.ndist; i++) mstate.delta[i] = delta_args[i];
        } else {
             printf("Error: -delta count must be 1 or ndist\n");
             exit(1);
        }
    } else {
        for(int i=0; i<config.ndist; i++) mstate.delta[i] = delta_default;
    }
    if (delta_args) free(delta_args);
    
    load_categories(&config, &gdata);

    if (config.mcmc) {
        // Logging continued
        if (config.fp_log) {
            fprintf(config.fp_log, "Phenotype column               = %8d\n", config.trait_pos);
            fprintf(config.fp_log, "No. of loci                    = %8d\n", gdata.nloci);
            fprintf(config.fp_log, "No. of individuals             = %8d\n", gdata.nind);
            fprintf(config.fp_log, "No. of training individuals    = %8d\n", gdata.nt);
            // ...
            fflush(config.fp_log);
        }
    }
    
    load_snp_binary(&config, &gdata);
    xcenter(&config, &gdata);
    init_random_seed_custom(&config, &rs);

    if (config.mcmc) {
        mstate.nnind = (double)gdata.nt;
        
        // Open output files
        if (config.snpout) config.fp_loc = fopen(config.locfil, "w");
        if (config.cat) config.fp_cat = fopen(config.catfil, "w");
        if (config.beta) config.fp_beta = fopen(config.betafil, "w");
        
        config.fp_hyp = fopen(config.hypfil, "w");
        if (config.fp_hyp) {
             fprintf(config.fp_hyp, " Replicate       Nsnp              Va              Ve "); // Header matching Fortran
             // ...
             fprintf(config.fp_hyp, "\n");
        }
        
        // VCE Logic
        if (config.dfvara < -2.0) {
            config.VCE = false;
            // Replicate Fortran logic: calc yhat, vary, update vara
            double sum_y = 0.0;
            int cnt = 0;
            for(int i=0; i<gdata.nind; i++) {
                if (gdata.trains[i] == 0) {
                    sum_y += gdata.why[i];
                    cnt++;
                }
            }
            mstate.yhat = sum_y / mstate.nnind;
            
            double sum_sq = 0.0;
            for(int i=0; i<gdata.nind; i++) {
                if (gdata.trains[i] == 0) {
                    double diff = gdata.why[i] - mstate.yhat;
                    sum_sq += diff * diff;
                }
            }
            mstate.vary = sum_sq / (mstate.nnind - 1.0);
            mstate.vara = mstate.vara * mstate.vary;
        } else {
            config.VCE = true;
            config.vara_ap = mstate.vara;
            config.vare_ap = mstate.vare;
            if (config.dfvara == -2.0) config.vara_ap = 0.0;
            if (config.dfvare == -2.0) config.vare_ap = 0.0;
        }

        run_mcmc(&config, &gdata, &mstate, &mstore, &rs);
        
    } else {
        // Predict mode
        // load_param
        // compute_dgv
        // write_dgv
        printf("Prediction mode not fully implemented in C port yet.\n");
    }

    // Cleanup
    if (config.fp_log) fclose(config.fp_log);
    // Free memory...

    return 0;
}
