#include "../../new_src_c/bayesRCO.h"
#include "../../new_src_c/mcmc_utils.h"
#include "../../new_src_c/mcmc_bayesCpi.h"
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
    manual_seed(&rs, s_packed);

    // Setup minimal config
    config.ncat = 2;
    config.ndist = 4;
    config.VCE = true;
    config.dfvara = 4.0;
    config.dfvare = 4.0;
    config.vara_ap = 0.01;
    config.vare_ap = 0.01;
    config.msize = 0;

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
    gdata.vsnptrack = (int*)calloc(gdata.nloci, sizeof(int));
    gdata.gannot = (double**)calloc(gdata.nloci, sizeof(double*));
    for (i = 0; i < gdata.nloci; i++) gdata.gannot[i] = (double*)calloc(config.ncat, sizeof(double));
    gdata.includedloci = (double*)calloc(gdata.nloci, sizeof(double));

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
    mstate.s = (double*)calloc(config.ndist, sizeof(double));
    mstate.stemp = (double*)calloc(config.ndist, sizeof(double));

    mstore.gstore = (double*)calloc(gdata.nloci, sizeof(double));
    mstore.pstore = (double**)calloc(config.ndist, sizeof(double*));
    for (i = 0; i < config.ndist; i++) mstore.pstore[i] = (double*)calloc(config.ncat, sizeof(double));
    mstore.annotstore = (double**)calloc(gdata.nloci, sizeof(double*));
    for (i = 0; i < gdata.nloci; i++) mstore.annotstore[i] = (double*)calloc(config.ncat, sizeof(double));

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

    // Fill C with test pattern
    for (j = 0; j < gdata.nloci; j++) {
        for (i = 0; i < config.ncat; i++) gdata.C[j][i] = 0;
        gdata.C[j][j % config.ncat] = 1;
    }

    // Initialize permvec
    for (j = 0; j < gdata.nloci; j++) gdata.permvec[j] = j;

    // Compute nannot
    for (j = 0; j < gdata.nloci; j++) {
        gdata.nannot[j] = 0;
        for (i = 0; i < config.ncat; i++) gdata.nannot[j] += gdata.C[j][i];
    }

    // Compute xpx
    for (j = 0; j < gdata.nloci; j++) {
        double sum = 0.0;
        for (i = 0; i < gdata.nt; i++) {
            double val = gdata.X[i * gdata.nloci + j];
            sum += val * val;
        }
        gdata.xpx[j] = sum;
    }

    // Initialize p
    for (j = 0; j < config.ncat; j++) {
        mstate.p[0][j] = 0.5;
        mstate.p[1][j] = 0.25;
        mstate.p[2][j] = 0.15;
        mstate.p[3][j] = 0.10;
    }

    // Initialize log_p
    for (j = 0; j < config.ncat; j++) {
        for (i = 0; i < config.ndist; i++) {
            mstate.log_p[i][j] = log(mstate.p[i][j]);
        }
    }

    // Initialize gp and vare_gp
    for (i = 0; i < config.ndist; i++) mstate.gp[i] = mstate.gpin[i] * mstate.vara;
    for (i = 1; i < config.ndist; i++) {
        mstate.log_gp[i] = log(mstate.gp[i]);
        mstate.vare_gp[i] = mstate.vare / mstate.gp[i];
    }

    // Initialize g and gannot
    for (j = 0; j < gdata.nloci; j++) mstate.g[j] = 0.1;
    for (j = 0; j < gdata.nloci; j++) {
        for (i = 0; i < config.ncat; i++) gdata.gannot[j][i] = 0.0;
    }

    // Initialize yadj
    for (i = 0; i < gdata.nt; i++) {
        mstate.yadj[i] = gdata.why[i] - 3.0;
    }

    // Run init
    mcmc_bayesCpi_init(&config, &gdata, &mstate, &mstore);

    printf("gannot_after_init:\n");
    for (j = 0; j < gdata.nloci; j++) {
        for (i = 0; i < config.ncat; i++) {
            printf("%20.16f\n", gdata.gannot[j][i]);
        }
    }

    printf("snptracker_after_init:\n");
    for (j = 0; j < gdata.nloci; j++) {
        for (i = 0; i < config.ncat; i++) {
            printf("%4d\n", gdata.snptracker[j][i]);
        }
    }

    // Run one iteration of bayesCpi kernel
    mstate.rep = 1;
    mcmc_bayesCpi_kernel(&config, &gdata, &mstate, &mstore, &rs);

    printf("g_after_kernel:\n");
    for (j = 0; j < gdata.nloci; j++) {
        printf("%20.16f\n", mstate.g[j]);
    }

    printf("included_after_kernel:\n");
    printf("%6d\n", mstate.included);

    return 0;
}
