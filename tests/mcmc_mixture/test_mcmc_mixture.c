#include "../../new_src_c/bayesRCO.h"
#include "../../new_src_c/mcmc_utils.h"
#include "../../new_src_c/mcmc_mixture.h"
#include "../../new_src_c/rng.h"
#include "../../new_src_c/utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

int main() {
    ModelConfig config;
    GenomicData gdata;
    MCMCState mstate;
    MCMCStorage mstore;
    prng_state rs;
    int i, j;

    memset(&config, 0, sizeof(ModelConfig));
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
    config.num_categories = 2;
    config.num_distributions = 4;
    config.VCE = true;
    config.dfvara = 4.0;
    config.dfvare = 4.0;
    config.vara_ap = 0.01;
    config.vare_ap = 0.01;
    config.marker_set_size = 0;

    // Setup minimal gdata
    gdata.num_loci = 10;
    gdata.num_individuals = 5;
    gdata.num_phenotyped_individuals = 5;

    // Allocate arrays
    gdata.genotypes = (double*)calloc(gdata.num_phenotyped_individuals * gdata.num_loci, sizeof(double));
    gdata.snp_correlations = (double*)calloc(gdata.num_loci, sizeof(double));
    gdata.annotations_per_locus = (int*)calloc(gdata.num_loci, sizeof(int));
    gdata.categories = (int**)calloc(gdata.num_loci, sizeof(int*));
    for (i = 0; i < gdata.num_loci; i++) gdata.categories[i] = (int*)calloc(config.num_categories, sizeof(int));
    gdata.permvec = (int*)calloc(gdata.num_loci, sizeof(int));
    gdata.phenotypes = (double*)calloc(gdata.num_individuals, sizeof(double));
    gdata.trains = (int*)calloc(gdata.num_individuals, sizeof(int));
    gdata.distribution_per_category = (int**)calloc(gdata.num_loci, sizeof(int*));
    for (i = 0; i < gdata.num_loci; i++) gdata.distribution_per_category[i] = (int*)calloc(config.num_categories, sizeof(int));
    gdata.current_distribution = (int*)calloc(gdata.num_loci, sizeof(int));
    gdata.current_category = (int*)calloc(gdata.num_loci, sizeof(int));
    gdata.atemp = (int*)calloc(config.num_categories, sizeof(int));

    mstate.snp_effects = (double*)calloc(gdata.num_loci, sizeof(double));
    mstate.adjusted_phenotypes = (double*)calloc(gdata.num_phenotyped_individuals, sizeof(double));
    mstate.genomic_values = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.variance_scaling_factors = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.p = (double**)calloc(config.num_distributions, sizeof(double*));
    for (i = 0; i < config.num_distributions; i++) mstate.p[i] = (double*)calloc(config.num_categories, sizeof(double));
    mstate.log_p = (double**)calloc(config.num_distributions, sizeof(double*));
    for (i = 0; i < config.num_distributions; i++) mstate.log_p[i] = (double*)calloc(config.num_categories, sizeof(double));
    mstate.log_distribution_variances = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.residual_variance_over_distribution_variances = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.dirichlet_priors = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.dirichlet_scratch = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.category_dirichlet_scratch = (double*)calloc(config.num_categories, sizeof(double));
    mstate.category_probabilities = (double*)calloc(config.num_categories, sizeof(double));
    mstate.snps_per_distribution = (int**)calloc(config.num_distributions, sizeof(int*));
    for (i = 0; i < config.num_distributions; i++) mstate.snps_per_distribution[i] = (int*)calloc(config.num_categories, sizeof(int));
    mstate.variance_per_distribution = (double**)calloc(config.num_distributions, sizeof(double*));
    for (i = 0; i < config.num_distributions; i++) mstate.variance_per_distribution[i] = (double*)calloc(config.num_categories, sizeof(double));
    mstate.ytemp = (double*)calloc(gdata.num_phenotyped_individuals, sizeof(double));
    mstate.log_likelihoods = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.selection_probs = (double*)calloc(config.num_distributions, sizeof(double));
    mstate.ss = (double*)calloc(config.num_categories, sizeof(double));
    mstate.sstemp = (double*)calloc(config.num_categories, sizeof(double));

    mstore.sum_snp_effects = (double*)calloc(gdata.num_loci, sizeof(double));
    mstore.sum_mixture_proportions = (double**)calloc(config.num_distributions, sizeof(double*));
    for (i = 0; i < config.num_distributions; i++) mstore.sum_mixture_proportions[i] = (double*)calloc(config.num_categories, sizeof(double));
    mstore.sum_snps_per_distribution = (double**)calloc(config.num_distributions, sizeof(double*));
    for (i = 0; i < config.num_distributions; i++) mstore.sum_snps_per_distribution[i] = (double*)calloc(config.num_categories, sizeof(double));
    mstore.sum_variance_per_distribution = (double**)calloc(config.num_distributions, sizeof(double*));
    for (i = 0; i < config.num_distributions; i++) mstore.sum_variance_per_distribution[i] = (double*)calloc(config.num_categories, sizeof(double));
    mstore.varistore = (double*)calloc(gdata.num_loci, sizeof(double));
    mstore.varustore = (double*)calloc(gdata.num_loci, sizeof(double));
    mstore.sum_distribution_counts = (double**)calloc(gdata.num_loci, sizeof(double*));
    for (i = 0; i < gdata.num_loci; i++) mstore.sum_distribution_counts[i] = (double*)calloc(config.num_distributions, sizeof(double));
    mstore.sum_category_counts = (double**)calloc(gdata.num_loci, sizeof(double*));
    for (i = 0; i < gdata.num_loci; i++) mstore.sum_category_counts[i] = (double*)calloc(config.num_categories, sizeof(double));

    // Initialize test values
    for (i = 0; i < gdata.num_individuals; i++) gdata.trains[i] = 0;
    double why_vals[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    for (i = 0; i < gdata.num_individuals; i++) gdata.phenotypes[i] = why_vals[i];
    mstate.nnind = 5.0;
    mstate.variance_genetic = 1.0;
    mstate.variance_residual = 1.0;
    mstate.mu = 0.0;
    double gpin_vals[] = {0.0, 0.0001, 0.001, 0.01};
    for (i = 0; i < config.num_distributions; i++) mstate.variance_scaling_factors[i] = gpin_vals[i];
    double delta_vals[] = {1.0, 1.0, 1.0, 1.0};
    for (i = 0; i < config.num_distributions; i++) mstate.dirichlet_priors[i] = delta_vals[i];

    // Fill X with simple test pattern (row-major in C)
    for (i = 0; i < gdata.num_phenotyped_individuals; i++) {
        for (j = 0; j < gdata.num_loci; j++) {
            gdata.genotypes[i * gdata.num_loci + j] = (double)(i + 1 + j + 1) / 10.0;
        }
    }

    // Fill C with test pattern - some loci have multiple annotations
    for (j = 0; j < gdata.num_loci; j++) {
        gdata.categories[j][0] = 1;  
        if (j % 2 == 1) gdata.categories[j][1] = 1;  
    }

    // Initialize permvec
    for (j = 0; j < gdata.num_loci; j++) {
        gdata.permvec[j] = j;
    }

    // Compute nannot
    for (j = 0; j < gdata.num_loci; j++) {
        gdata.annotations_per_locus[j] = 0;
        for (i = 0; i < config.num_categories; i++) gdata.annotations_per_locus[j] += gdata.categories[j][i];
    }

    // Compute xpx
    for (j = 0; j < gdata.num_loci; j++) {
        gdata.snp_correlations[j] = dot_product_col(gdata.genotypes, j, gdata.num_phenotyped_individuals, gdata.num_loci, &gdata.genotypes[j]);
        // Manual calculation for column dot product (sanity check logic preserved)
        double sum = 0.0;
        for (i = 0; i < gdata.num_phenotyped_individuals; i++) {
            double val = gdata.genotypes[i * gdata.num_loci + j];
            sum += val * val;
        }
        gdata.snp_correlations[j] = sum;
    }

    // Initialize p
    for (j = 0; j < config.num_categories; j++) {
        mstate.p[0][j] = 0.5;
        mstate.p[1][j] = 0.25;
        mstate.p[2][j] = 0.15;
        mstate.p[3][j] = 0.10;
    }

    // Initialize log_p
    for (j = 0; j < config.num_categories; j++) {
        for (i = 0; i < config.num_distributions; i++) {
            mstate.log_p[i][j] = log(mstate.p[i][j]);
        }
    }

    // Initialize gp and vare_gp
    for (i = 0; i < config.num_distributions; i++) mstate.genomic_values[i] = mstate.variance_scaling_factors[i] * mstate.variance_genetic;
    for (i = 1; i < config.num_distributions; i++) {
        mstate.log_distribution_variances[i] = log(mstate.genomic_values[i]);
        mstate.residual_variance_over_distribution_variances[i] = mstate.variance_residual / mstate.genomic_values[i];
    }

    // Initialize g
    for (j = 0; j < gdata.num_loci; j++) mstate.snp_effects[j] = 0.1;

    // Initialize yadj (residuals)
    for (i = 0; i < gdata.num_phenotyped_individuals; i++) {
        mstate.adjusted_phenotypes[i] = gdata.phenotypes[i] - 3.0;  // y - mean
    }

    // Run mixture init
    mcmc_mixture_init(&config, &gdata, &mstate);

    printf("vsnptrack_after_init:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%4d\n", gdata.current_distribution[j]);
    }

    printf("a_after_init:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%4d\n", gdata.current_category[j]);
    }

    // Run one iteration of mixture kernel
    mstate.rep = 1;
    mcmc_mixture_kernel(&config, &gdata, &mstate, &mstore, &rs);

    printf("g_after_kernel:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        printf("%20.16f\n", mstate.snp_effects[j]);
    }

    printf("snpindist_after_kernel:\n");
    for (j = 0; j < config.num_categories; j++) {
        for (i = 0; i < config.num_distributions; i++) {
            printf("%6d\n", mstate.snps_per_distribution[i][j]);
        }
    }

    printf("included_after_kernel:\n");
    printf("%6d\n", mstate.included);

    printf("snptracker:\n");
    for (j = 0; j < gdata.num_loci; j++) {
        for (i = 0; i < config.num_categories; i++) {
            printf("%4d\n", gdata.distribution_per_category[j][i]);
        }
    }

    return 0;
}
