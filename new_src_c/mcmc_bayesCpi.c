#include "mcmc_bayesCpi.h"
#include "mcmc_utils.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

void mcmc_bayesCpi_init(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    // mstore->annotstore = dble(gdata->C(:, 1:ncat))
    for (int i = 0; i < gdata->nloci; i++) {
        for (int j = 0; j < config->ncat; j++) {
            mstore->annotstore[i][j] = (double)gdata->C[i][j];
        }
    }
    
    // gdata->gannot = 0.0
    for (int i = 0; i < gdata->nloci; i++) {
        for (int j = 0; j < config->ncat; j++) {
            gdata->gannot[i][j] = 0.0;
        }
    }
    
    // where (C(:, j) == 1) snptracker(:, j) = 2
    for (int j = 0; j < config->ncat; j++) {
        for (int i = 0; i < gdata->nloci; i++) {
            if (gdata->C[i][j] == 1) {
                gdata->snptracker[i][j] = 2;
            }
        }
    }
}

void mcmc_bayesCpi_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int i, j, k, kk, snploc, l;
    double skk, sk, clike, ssculm, r;
    bool overflow;

    // For each category
    for (j = 0; j < config->ncat; j++) {
        // Initialize log_p for this category
        for (i = 0; i < config->ndist; i++) {
            mstate->log_p[i][j] = log(mstate->p[i][j]);
        }

        // For each locus (in permuted order)
        for (k = 0; k < gdata->nloci; k++) {
            snploc = gdata->permvec[k];
            
            // Only process if this locus is in this category
            if (gdata->C[snploc][j] == 1) {
                mstate->zz = gdata->xpx[snploc];
                mstate->zz_vare = mstate->zz / mstate->vare;
                mstate->gk = gdata->gannot[snploc][j];
                
                // If previously included, add back effect to yadj
                if (gdata->snptracker[snploc][j] > 1) {
                    add_col_scalar(mstate->yadj, gdata->X, snploc, gdata->nt, gdata->nloci, mstate->gk);
                }
                
                // Compute rhs
                mstate->rhs = dot_product_col(gdata->X, snploc, gdata->nt, gdata->nloci, mstate->yadj);
                
                // Compute selection probabilities
                mstate->s[0] = mstate->log_p[0][j];
                for (kk = 1; kk < config->ndist; kk++) {
                    mstate->logdetV = log(mstate->gp[kk] * mstate->zz_vare + 1.0);
                    mstate->uhat = mstate->rhs / (mstate->zz + mstate->vare_gp[kk]);
                    mstate->s[kk] = -0.5 * (mstate->logdetV - 
                        (mstate->rhs * mstate->uhat / mstate->vare)) + mstate->log_p[kk][j];
                }
                
                // Stabilize and convert to probabilities
                for (kk = 0; kk < config->ndist; kk++) {
                    skk = mstate->s[kk];
                    sk = 0.0;
                    overflow = false;
                    for (l = 0; l < config->ndist; l++) {
                        if (l == kk) continue;
                        clike = mstate->s[l] - skk;
                        if (clike > LOG_UPPER_LIMIT) {
                            overflow = true;
                            break;
                        }
                        if (clike < -LOG_UPPER_LIMIT) continue;
                        sk += exp(clike);
                    }
                    if (overflow) {
                        mstate->stemp[kk] = 0.0;
                    } else {
                        mstate->stemp[kk] = 1.0 / (1.0 + sk);
                    }
                }
                
                // Sample distribution
                ssculm = 0.0;
                r = rand_uniform(rs, 0.0, 1.0);
                config->indistflag = 0;  // 0-based (dist 1 in Fortran = 0 in C)
                for (kk = 0; kk < config->ndist; kk++) {
                    ssculm += mstate->stemp[kk];
                    if (r < ssculm) {
                        config->indistflag = kk;
                        break;
                    }
                }
                
                // Store tracker (1-based for compatibility)
                gdata->snptracker[snploc][j] = config->indistflag + 1;
                gdata->vsnptrack[snploc] = config->indistflag + 1;
                
                // Sample effect
                if (config->indistflag == 0) {
                    mstate->gk = 0.0;
                } else {
                    mstate->v1 = mstate->zz + mstate->vare / mstate->gp[config->indistflag];
                    mstate->gk = rand_normal(rs, mstate->rhs / mstate->v1, sqrt(mstate->vare / mstate->v1));
                    add_col_scalar(mstate->yadj, gdata->X, snploc, gdata->nt, gdata->nloci, -mstate->gk);
                    mstate->included++;
                }
                
                gdata->gannot[snploc][j] = mstate->gk;
                
                // Early exit if msize reached
                if (config->msize > 0 && mstate->rep > config->mrep) {
                    if (mstate->included >= config->msize) break;
                }
            }
        }
    }
    
    // Sum loci effects: g = sum(gannot, dim=2)
    for (i = 0; i < gdata->nloci; i++) {
        double sum = 0.0;
        for (j = 0; j < config->ncat; j++) {
            sum += gdata->gannot[i][j];
        }
        mstate->g[i] = sum;
    }
    
    // Statistics: snpindist and varindist
    for (j = 0; j < config->ncat; j++) {
        for (i = 0; i < config->ndist; i++) {
            int cnt = 0;
            double sum_g2 = 0.0;
            for (k = 0; k < gdata->nloci; k++) {
                if (gdata->snptracker[k][j] == i + 1) {  // 1-based comparison
                    cnt++;
                    sum_g2 += gdata->gannot[k][j] * gdata->gannot[k][j];
                }
            }
            mstate->snpindist[i][j] = cnt;
            mstate->varindist[i][j] = sum_g2;
        }
    }
    
    // Count included loci (any snptracker > 1 means included)
    for (i = 0; i < gdata->nloci; i++) {
        gdata->includedloci[i] = 0.0;
    }
    for (j = 0; j < config->ncat; j++) {
        for (i = 0; i < gdata->nloci; i++) {
            if (gdata->snptracker[i][j] > 1) {
                gdata->includedloci[i] = 1.0;
            }
        }
    }
    double included_sum = 0.0;
    for (i = 0; i < gdata->nloci; i++) {
        included_sum += gdata->includedloci[i];
    }
    mstate->included = (int)included_sum;
}
