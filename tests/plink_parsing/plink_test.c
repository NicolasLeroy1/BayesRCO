#include "bayesrco_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    IOConfig ioconfig;
    GenomicData gdata;

    memset(&ioconfig, 0, sizeof(IOConfig));
    memset(&gdata, 0, sizeof(GenomicData));
    
    strcpy(ioconfig.geno_file_path, "test_data.bed");
    strcpy(ioconfig.pheno_file_path, "test_data.fam");
    ioconfig.trait_column_index = 1;
    ioconfig.mcmc = true;
    io_get_size(&ioconfig, &gdata);
    io_load_phenotypes(&ioconfig, &gdata);
    io_load_genotypes(&ioconfig, &gdata);

    printf("Phenotypes (why):\n");
    for (int i = 0; i < gdata.num_individuals; i++) {
        if (gdata.phenotypes[i] == MISSING_VALUE) {
            printf(" NA\n");
        } else {
            printf("%20.6f\n", gdata.phenotypes[i]);
        }
    }

    printf("Genotypes (X):\n");
    for (int i = 0; i < gdata.num_phenotyped_individuals; i++) {
        for (int j = 0; j < gdata.num_loci; j++) {
            printf("%5.1f\n", gdata.genotypes[j * gdata.num_phenotyped_individuals + i]);
        }
    }

    io_cleanup(&ioconfig, &gdata, NULL, NULL);

    return 0;
}
