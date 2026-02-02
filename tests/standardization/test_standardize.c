#include "../../new_src_c/bayesRCO.h"
#include "../../new_src_c/io.h"
#include "../../new_src_c/utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    ModelConfig config;
    GenomicData gdata;

    memset(&config, 0, sizeof(ModelConfig));
    memset(&gdata, 0, sizeof(GenomicData));
    
    strcpy(config.inprefix, "test_data");
    strcpy(config.phenfil, "test_data.fam");
    strcpy(config.bimfil, "test_data.bim");
    strcpy(config.genfil, "test_data.bed");
    strcpy(config.freqfil, "test_data.frq");
    config.trait_pos = 1;
    config.mcmc = true;

    get_size(&config, &gdata);
    load_phenos_plink(&config, &gdata);
    
    // Manual allocation for test
    gdata.nt = 0;
    for(int i=0; i<gdata.nind; i++) if(gdata.trains[i]==0) gdata.nt++;
    gdata.X = (double*)calloc(gdata.nt * gdata.nloci, sizeof(double));
    gdata.freqstore = (double*)calloc(gdata.nloci, sizeof(double));
    
    load_snp_binary(&config, &gdata);

    // Run xcenter
    xcenter(&config, &gdata);

    printf("Allele Frequencies (freqstore):\n");
    for (int j = 0; j < gdata.nloci; j++) {
        printf("%20.16f\n", gdata.freqstore[j]);
    }

    printf("Standardized Genotypes (X):\n");
    for (int i = 0; i < gdata.nt; i++) {
        for (int j = 0; j < gdata.nloci; j++) {
            printf("%20.16f\n", gdata.X[i * gdata.nloci + j]);
        }
    }

    return 0;
}
