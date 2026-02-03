#include "bayesrco_io.h"
#include "mcmc_utils.h"
#include "mcmc_additive.h"
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

    // Seed RNG for reproducibility
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
    ioconfig.model.marker_set_size = 0;
    ioconfig.mcmc = true;

    // Setup minimal gdata dimensions
    gdata.num_loci = 10;
    gdata.num_individuals = 5;
    
    // Allocate trains manually so io_allocate_data can compute num_phenotyped_individuals
    gdata.trains = (int*)calloc(gdata.num_individuals, sizeof(int));
    for (i = 0; i < gdata.num_individuals; i++) gdata.trains[i] = 0;

    // Allocate arrays
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

    // Fill categories with test pattern - add some overlapping annotations
    int ncat = ioconfig.model.num_categories;
    for (j = 0; j < gdata.num_loci; j++) {
        for (i = 0; i < ncat; i++) gdata.categories[j * ncat + i] = 0;
        gdata.categories[j * ncat + (j % ncat)] = 1;
        // Add overlap: locus 2 is in both categories
        if (j == 2) gdata.categories[j * ncat + 1] = 1;
    }

    // Initialize permvec and permannot
    for (j = 0; j < gdata.num_loci; j++) gdata.permvec[j] = j;
    for (j = 0; j < ncat; j++) gdata.permannot[j] = j;

    // Compute nannot and xpx
    for (j = 0; j < gdata.num_loci; j++) {
        int sum_c = 0;
        for (i = 0; i < ncat; i++) sum_c += gdata.categories[j * ncat + i];
        gdata.annotations_per_locus[j] = sum_c;
        
        double sum_xx = 0.0;
        for (i = 0; i < nt; i++) {
            double val = gdata.genotypes[j * nt + i];
            sum_xx += val * val;
        }
        gdata.snp_correlations[j] = sum_xx;
    }

    // Initialize p and log_p
    int ndist = ioconfig.model.num_distributions;
    for (j = 0; j < ncat; j++) {
        mstate.p[IDX2(0, j, ncat)] = 0.5;
        mstate.p[IDX2(1, j, ncat)] = 0.25;
        mstate.p[IDX2(2, j, ncat)] = 0.15;
        mstate.p[IDX2(3, j, ncat)] = 0.10;
        
        for (i = 0; i < ndist; i++) {
            mstate.log_p[IDX2(i, j, ncat)] = log(mstate.p[IDX2(i, j, ncat)]);
        }
    }

    // Initialize gp and vare_gp
    for (i = 0; i < ndist; i++) mstate.distribution_variances[i] = mstate.variance_scaling_factors[i] * mstate.variance_genetic;
    for (i = 1; i < ndist; i++) {
        mstate.log_distribution_variances[i] = log(mstate.distribution_variances[i]);
        mstate.residual_variance_over_distribution_variances[i] = mstate.variance_residual / mstate.distribution_variances[i];
    }

    // Initialize g, gannot, vsnptrack, snptracker
    for (j = 0; j < gdata.num_loci; j++) mstate.snp_effects[j] = 0.1;
    for (j = 0; j < gdata.num_loci; j++) {
        for (i = 0; i < ncat; i++) {
            gdata.effects_per_category[j * ncat + i] = 0.0;
            gdata.distribution_per_category[j * ncat + i] = 2; // BayesC like
        }
        gdata.current_distribution[j] = 2;
    }

    // Initialize yadj
    for (i = 0; i < nt; i++) {
        mstate.adjusted_phenotypes[i] = gdata.phenotypes[i] - 3.0;
    }

    // Run one iteration of additive kernel
    mstate.rep = 1;
    mcmc_additive_kernel(&(ioconfig.model), &gdata, &mstate, &mstore, &rs);

    printf("g_after_kernel:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%25.16E\n", mstate.snp_effects[j]);
    }

    printf("gannot:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        for (i = 0; i < ncat; i++) {
            printf("%25.16E\n", gdata.effects_per_category[j * ncat + i]);
        }
    }

    printf("snpindist_after_kernel:\n");
    for (j = 0; j < ncat; j++) {
        for (i = 0; i < ndist; i++) {
            printf("%6d\n", mstate.snps_per_distribution[IDX2(i, j, ncat)]);
        }
    }

    printf("included_after_kernel:\n");
    printf("%6d\n", mstate.included);

    printf("vsnptrack:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%4d\n", gdata.current_distribution[j]);
    }

    // Cleanup
    io_cleanup(&ioconfig, &gdata, &mstate, &mstore);

    return 0;
}
