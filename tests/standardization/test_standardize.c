#include "bayesrco_io.h"
#include "utils.h"
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

    // Load dimensions
    if (io_get_size(&ioconfig, &gdata) != SUCCESS) {
        fprintf(stderr, "Error in io_get_size\n");
        return 1;
    }
    
    // Load phenotypes (implicitly allocates trains)
    if (io_load_phenotypes(&ioconfig, &gdata) != SUCCESS) {
        fprintf(stderr, "Error in io_load_phenotypes\n");
        return 1;
    }
    
    // Allocate logic for X and frequencies
    int nt = gdata.num_phenotyped_individuals;
    gdata.genotypes = (double*)calloc(gdata.num_loci * nt, sizeof(double));
    gdata.allele_frequencies = (double*)calloc(gdata.num_loci, sizeof(double));
    
    // Load genotypes
    if (io_load_genotypes(&ioconfig, &gdata) != SUCCESS) {
        fprintf(stderr, "Error in io_load_genotypes\n");
        return 1;
    }

    // Run standardization
    standardize_genotypes(&gdata, true);

    printf("Allele Frequencies (freqstore):\n");
    for (int j = 0; j < gdata.num_loci; j++) {
        printf("%25.16E\n", gdata.allele_frequencies[j]);
    }

    printf("Standardized Genotypes (X):\n");
    for (int i = 0; i < nt; i++) {
        for (int j = 0; j < gdata.num_loci; j++) {
            // C genotypes are column-major (nloci x nt), so X[j * nt + i]
            printf("%25.16E\n", gdata.genotypes[j * nt + i]);
        }
    }

    // Cleanup
    io_cleanup(&ioconfig, &gdata, NULL, NULL);

    return 0;
}
