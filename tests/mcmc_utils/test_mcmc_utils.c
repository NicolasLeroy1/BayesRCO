#include "bayesrco_io.h"
#include "mcmc_utils.h"
#include "rng.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

int main() {
    IOConfig ioconfig;
    GenomicData gdata;
    MCMCState mstate;
    MCMCStorage mstore;
    prng_state rs;
    int i, j;

    memset(&ioconfig, 0, sizeof(IOConfig));
    memset(&gdata, 0, sizeof(GenomicData));
    memset(&mstate, 0, sizeof(MCMCState));
    memset(&mstore, 0, sizeof(MCMCStorage));

    // Seed RNG for reproducibility (matching Fortran seed pattern)
    uint32_t p[8];
    for (i = 0; i < 8; i++) p[i] = 12345 + i;
    uint64_t s_packed[4];
    s_packed[0] = ((uint64_t)p[6] << 32) | (uint64_t)p[7];
    s_packed[1] = ((uint64_t)p[4] << 32) | (uint64_t)p[5];
    s_packed[2] = ((uint64_t)p[2] << 32) | (uint64_t)p[3];
    s_packed[3] = ((uint64_t)p[0] << 32) | (uint64_t)p[1];
    rng_seed(&rs, s_packed);

    // Setup minimal config
    ioconfig.model.num_categories = 2;
    ioconfig.model.num_distributions = 4;
    ioconfig.model.VCE = true;
    ioconfig.model.dfvara = 4.0;
    ioconfig.model.dfvare = 4.0;
    ioconfig.model.vara_ap = 0.01;
    ioconfig.model.vare_ap = 0.01;
    ioconfig.mcmc = true;

    // Setup minimal gdata dimensions
    gdata.num_loci = 10;
    gdata.num_individuals = 5;
    
    // Allocate trains manually so io_allocate_data can compute num_phenotyped_individuals
    gdata.trains = (int*)calloc(gdata.num_individuals, sizeof(int));
    for (i = 0; i < gdata.num_individuals; i++) gdata.trains[i] = 0;

    // Allocate others via io_allocate_data
    io_allocate_data(&ioconfig, &gdata, &mstate, &mstore);

    // Initialize test values
    double why_vals[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    gdata.phenotypes = (double*)calloc(gdata.num_individuals, sizeof(double));
    for (i = 0; i < gdata.num_individuals; i++) gdata.phenotypes[i] = why_vals[i];
    
    mstate.nnind = 5.0;
    mstate.variance_genetic = 1.0;
    mstate.variance_residual = 1.0;
    mstate.mu = 0.0;
    
    double vscaling[] = {0.0, 0.0001, 0.001, 0.01};
    for (i = 0; i < ioconfig.model.num_distributions; i++) mstate.variance_scaling_factors[i] = vscaling[i];
    
    mstate.dirichlet_priors = (double*)calloc(ioconfig.model.num_distributions, sizeof(double));
    for (i = 0; i < ioconfig.model.num_distributions; i++) mstate.dirichlet_priors[i] = 1.0;

    // Fill genotypes with simple test pattern (nloci x nt, column-major)
    int nt = gdata.num_phenotyped_individuals;
    gdata.genotypes = (double*)calloc(gdata.num_loci * nt, sizeof(double));
    for (j = 0; j < gdata.num_loci; j++) {
        for (i = 0; i < nt; i++) {
            // Fortran: gdata%X(i, j) = dble(i + j) / 10.0d0
            gdata.genotypes[j * nt + i] = (double)(i + 1 + j + 1) / 10.0;
        }
    }

    // Fill categories with test pattern (alternating annotations, size: nloci * ncat)
    int ncat = ioconfig.model.num_categories;
    for (j = 0; j < gdata.num_loci; j++) {
        for (i = 0; i < ncat; i++) gdata.categories[j * ncat + i] = 0;
        gdata.categories[j * ncat + (j % ncat)] = 1;
    }

    // Test mcmc_init_common
    mcmc_init_common(&(ioconfig.model), &gdata, &mstore);

    printf("xpx:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%25.16E\n", gdata.snp_correlations[j]);
    }

    printf("nannot:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%4d\n", gdata.annotations_per_locus[j]);
    }

    // Test mcmc_start_values_common
    mcmc_start_values_common(&(ioconfig.model), &gdata, &mstate, &rs);

    printf("yhat:\n");
    printf("%25.16E\n", mstate.yhat);

    printf("vary:\n");
    printf("%25.16E\n", mstate.vary);

    printf("g_init:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%25.16E\n", mstate.snp_effects[j]);
    }

    printf("p_init:\n");
    for (j = 0; j < ncat; j++) {
        for (i = 0; i < ioconfig.model.num_distributions; i++) {
            printf("%25.16E\n", mstate.p[i * ncat + j]);
        }
    }

    // Test mcmc_iteration_pre_common
    mstate.rep = 1;
    mcmc_iteration_pre_common(&(ioconfig.model), &mstate, &rs);

    printf("mu_after_pre:\n");
    printf("%25.16E\n", mstate.mu);

    printf("vare_after_pre:\n");
    printf("%25.16E\n", mstate.variance_residual);

    printf("log_gp:\n");
    for (i = 1; i < ioconfig.model.num_distributions; i++) {
        printf("%25.16E\n", mstate.log_distribution_variances[i]);
    }

    printf("vare_gp:\n");
    for (i = 1; i < ioconfig.model.num_distributions; i++) {
        printf("%25.16E\n", mstate.residual_variance_over_distribution_variances[i]);
    }

    // Cleanup
    io_cleanup(&ioconfig, &gdata, &mstate, &mstore);

    return 0;
}
