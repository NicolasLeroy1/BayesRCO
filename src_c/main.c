#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bayes_structs.h"
#include "stats.h"
#include "data_io.h"
#include "mcmc.h"

void parse_arguments(int argc, char **argv, BayesConfig *config) {
    // Defaults
    config->numit = 1000;
    config->burnin = 100;
    config->thin = 1;
    config->seed = 0;
    config->trait_pos = 1;
    config->ndist = 4;
    config->ncat = 1;
    config->vara = 0.01;
    config->vare = 0.01;
    config->dfvara = -2.0;
    config->dfvare = -2.0;
    
    strcpy(config->inprefix, "input");
    strcpy(config->outprefix, "output");
    config->catRC[0] = '\0';

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-bfile") == 0) strcpy(config->inprefix, argv[++i]);
        else if (strcmp(argv[i], "-out") == 0) strcpy(config->outprefix, argv[++i]);
        else if (strcmp(argv[i], "-numit") == 0) config->numit = atoi(argv[++i]);
        else if (strcmp(argv[i], "-burnin") == 0) config->burnin = atoi(argv[++i]);
        else if (strcmp(argv[i], "-thin") == 0) config->thin = atoi(argv[++i]);
        else if (strcmp(argv[i], "-seed") == 0) config->seed = atoi(argv[++i]);
        else if (strcmp(argv[i], "-vara") == 0) config->vara = atof(argv[++i]);
        else if (strcmp(argv[i], "-vare") == 0) config->vare = atof(argv[++i]);
        else if (strcmp(argv[i], "-catfile") == 0) strcpy(config->catRC, argv[++i]);
    }

    config->gpin = calloc(config->ndist, sizeof(double));
    config->delta = calloc(config->ndist, sizeof(double));
    if (config->ndist >= 4) {
        config->gpin[0] = 0.0;
        config->gpin[1] = 0.0001;
        config->gpin[2] = 0.001;
        config->gpin[3] = 0.01;
        
        for(int i=0; i<config->ndist; i++) config->delta[i] = 1.0;
    }

    if (config->catRC[0] != '\0') {
        // Simple column counter for first line
        FILE *fcat = fopen(config->catRC, "r");
        if (fcat) {
            char line[1024];
            if (fgets(line, sizeof(line), fcat)) {
                int cols = 0;
                char *p = strtok(line, " \t\n");
                while(p) { cols++; p = strtok(NULL, " \t\n"); }
                config->ncat = cols;
            }
            fclose(fcat);
        }
    }
}

int main(int argc, char **argv) {
    printf("BayesRCO Refactored C Version\n");

    BayesConfig config;
    parse_arguments(argc, argv, &config);

    char phenfil[512], bimfil[512], genfil[512];
    sprintf(phenfil, "%s.fam", config.inprefix);
    sprintf(bimfil, "%s.bim", config.inprefix);
    sprintf(genfil, "%s.bed", config.inprefix);

    BayesData data;
    if (get_data_size(phenfil, bimfil, &data.nind, &data.nloci) != 0) return 1;

    data.why = malloc(data.nind * sizeof(double));
    data.trains = malloc(data.nind * sizeof(int));
    if (load_phenotypes(phenfil, data.nind, config.trait_pos, data.why, data.trains, &data.nt) != 0) return 1;

    data.X = malloc((size_t)data.nloci * data.nt * sizeof(double));
    if (load_genotypes(genfil, data.nind, data.nloci, data.trains, data.nt, data.X) != 0) return 1;

    data.C = malloc(data.nloci * config.ncat * sizeof(int));
    load_categories(config.catRC, data.nloci, config.ncat, data.C);

    data.freq = malloc(data.nloci * sizeof(double));
    compute_frequencies_and_center(data.nt, data.nloci, data.X, data.freq);

    // Save frequencies for verification (to match verify.sh expected output)
    char freqfil[512];
    sprintf(freqfil, "%s.frq", config.outprefix);
    FILE *fp_freq = fopen(freqfil, "w");
    if (fp_freq) {
        for (int j = 0; j < data.nloci; j++) fprintf(fp_freq, "%10.6f\n", data.freq[j]);
        fclose(fp_freq);
    }

    init_rng(config.seed);

    BayesModel model;
    init_model(&config, &data, &model);

    run_mcmc(&config, &data, &model);

    // Cleanup
    free(config.gpin);
    free(config.delta);
    free_bayes_data(&data);
    free_bayes_model(&model);

    printf("Done.\n");
    return 0;
}
