#include "../../new_src_c/bayesRCO.h"
#include "../../new_src_c/mcmc_utils.h"
#include "../../new_src_c/rng.h"
#include "../../new_src_c/utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

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

    // Seed RNG for reproducibility (matching Fortran seed pattern)
    uint32_t p[8];
    for (i = 0; i < 8; i++) p[i] = 12345 + i;
    uint64_t s_packed[4];
    s_packed[0] = ((uint64_t)p[6] << 32) | (uint64_t)p[7];
    s_packed[1] = ((uint64_t)p[4] << 32) | (uint64_t)p[5];
    s_packed[2] = ((uint64_t)p[2] << 32) | (uint64_t)p[3];
    s_packed[3] = ((uint64_t)p[0] << 32) | (uint64_t)p[1];
    manual_seed(&rs, s_packed);

    // Setup minimal config
    config.ncat = 2;
    config.ndist = 4;
    config.VCE = true;
    config.dfvara = 4.0;
    config.dfvare = 4.0;
    config.vara_ap = 0.01;
    config.vare_ap = 0.01;

    // Setup minimal gdata
    gdata.nloci = 10;
    gdata.nind = 5;
    gdata.nt = 5;

    // Allocate arrays
    gdata.X = (double*)calloc(gdata.nt * gdata.nloci, sizeof(double));
    gdata.xpx = (double*)calloc(gdata.nloci, sizeof(double));
    gdata.nannot = (int*)calloc(gdata.nloci, sizeof(int));
    gdata.C = (int**)calloc(gdata.nloci, sizeof(int*));
    for (i = 0; i < gdata.nloci; i++) gdata.C[i] = (int*)calloc(config.ncat, sizeof(int));
    gdata.permvec = (int*)calloc(gdata.nloci, sizeof(int));
    gdata.why = (double*)calloc(gdata.nind, sizeof(double));
    gdata.trains = (int*)calloc(gdata.nind, sizeof(int));
    gdata.snptracker = (int**)calloc(gdata.nloci, sizeof(int*));
    for (i = 0; i < gdata.nloci; i++) gdata.snptracker[i] = (int*)calloc(config.ncat, sizeof(int));

    mstate.g = (double*)calloc(gdata.nloci, sizeof(double));
    mstate.yadj = (double*)calloc(gdata.nt, sizeof(double));
    mstate.gp = (double*)calloc(config.ndist, sizeof(double));
    mstate.gpin = (double*)calloc(config.ndist, sizeof(double));
    mstate.p = (double**)calloc(config.ndist, sizeof(double*));
    for (i = 0; i < config.ndist; i++) mstate.p[i] = (double*)calloc(config.ncat, sizeof(double));
    mstate.log_p = (double**)calloc(config.ndist, sizeof(double*));
    for (i = 0; i < config.ndist; i++) mstate.log_p[i] = (double*)calloc(config.ncat, sizeof(double));
    mstate.log_gp = (double*)calloc(config.ndist, sizeof(double));
    mstate.vare_gp = (double*)calloc(config.ndist, sizeof(double));
    mstate.delta = (double*)calloc(config.ndist, sizeof(double));
    mstate.dirx = (double*)calloc(config.ndist, sizeof(double));
    mstate.snpindist = (int**)calloc(config.ndist, sizeof(int*));
    for (i = 0; i < config.ndist; i++) mstate.snpindist[i] = (int*)calloc(config.ncat, sizeof(int));
    mstate.varindist = (double**)calloc(config.ndist, sizeof(double*));
    for (i = 0; i < config.ndist; i++) mstate.varindist[i] = (double*)calloc(config.ncat, sizeof(double));
    mstate.ytemp = (double*)calloc(gdata.nloci, sizeof(double));

    mstore.gstore = (double*)calloc(gdata.nloci, sizeof(double));
    mstore.pstore = (double**)calloc(config.ndist, sizeof(double*));
    for (i = 0; i < config.ndist; i++) mstore.pstore[i] = (double*)calloc(config.ncat, sizeof(double));
    mstore.snpstore = (double**)calloc(config.ndist, sizeof(double*));
    for (i = 0; i < config.ndist; i++) mstore.snpstore[i] = (double*)calloc(config.ncat, sizeof(double));
    mstore.varstore = (double**)calloc(config.ndist, sizeof(double*));
    for (i = 0; i < config.ndist; i++) mstore.varstore[i] = (double*)calloc(config.ncat, sizeof(double));
    mstore.varistore = (double*)calloc(gdata.nloci, sizeof(double));
    mstore.varustore = (double*)calloc(gdata.nloci, sizeof(double));

    // Initialize test values
    for (i = 0; i < gdata.nind; i++) gdata.trains[i] = 0;
    double why_vals[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    for (i = 0; i < gdata.nind; i++) gdata.why[i] = why_vals[i];
    mstate.nnind = 5.0;
    mstate.vara = 1.0;
    mstate.vare = 1.0;
    mstate.mu = 0.0;
    double gpin_vals[] = {0.0, 0.0001, 0.001, 0.01};
    for (i = 0; i < config.ndist; i++) mstate.gpin[i] = gpin_vals[i];
    double delta_vals[] = {1.0, 1.0, 1.0, 1.0};
    for (i = 0; i < config.ndist; i++) mstate.delta[i] = delta_vals[i];

    // Fill X with simple test pattern (row-major in C)
    for (i = 0; i < gdata.nt; i++) {
        for (j = 0; j < gdata.nloci; j++) {
            gdata.X[i * gdata.nloci + j] = (double)(i + 1 + j + 1) / 10.0;
        }
    }

    // Fill C with test pattern (alternating annotations)
    for (j = 0; j < gdata.nloci; j++) {
        for (i = 0; i < config.ncat; i++) gdata.C[j][i] = 0;
        gdata.C[j][j % config.ncat] = 1;
    }

    // Test mcmc_init_common
    mcmc_init_common(&config, &gdata, &mstore);

    printf("xpx:\n");
    for (j = 0; j < gdata.nloci; j++) {
        printf("%20.16f\n", gdata.xpx[j]);
    }

    printf("nannot:\n");
    for (j = 0; j < gdata.nloci; j++) {
        printf("%4d\n", gdata.nannot[j]);
    }

    // Test mcmc_start_values_common
    mcmc_start_values_common(&config, &gdata, &mstate, &rs);

    printf("yhat:\n");
    printf("%20.16f\n", mstate.yhat);

    printf("vary:\n");
    printf("%20.16f\n", mstate.vary);

    printf("g_init:\n");
    for (j = 0; j < gdata.nloci; j++) {
        printf("%20.16f\n", mstate.g[j]);
    }

    printf("p_init:\n");
    for (j = 0; j < config.ncat; j++) {
        for (i = 0; i < config.ndist; i++) {
            printf("%20.16f\n", mstate.p[i][j]);
        }
    }

    // Test mcmc_iteration_pre_common
    mstate.rep = 1;
    mcmc_iteration_pre_common(&config, &mstate, &rs);

    printf("mu_after_pre:\n");
    printf("%20.16f\n", mstate.mu);

    printf("vare_after_pre:\n");
    printf("%20.16f\n", mstate.vare);

    printf("log_gp:\n");
    for (i = 1; i < config.ndist; i++) {
        printf("%20.16f\n", mstate.log_gp[i]);
    }

    printf("vare_gp:\n");
    for (i = 1; i < config.ndist; i++) {
        printf("%20.16f\n", mstate.vare_gp[i]);
    }

    return 0;
}
