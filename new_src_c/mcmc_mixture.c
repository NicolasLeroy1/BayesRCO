#include "mcmc_mixture.h"
#include "mcmc_utils.h"
#include <math.h>
#include <stdbool.h>

void mcmc_mixture_init(ModelConfig *config, GenomicData *gdata, MCMCState *mstate) {
    // gdata->vsnptrack = 2
    for(int k=0; k<gdata->nloci; k++) gdata->vsnptrack[k] = 2; // ? Fortran: gdata%vsnptrack = 2
    
    for(int k=0; k<gdata->nloci; k++) {
        if (gdata->nannot[k] == 1) {
            for(int j=0; j<config->ncat; j++) {
                if (gdata->C[k][j] == 1) gdata->a[k] = j; // 0-based
            }
        }
    }
}

void mcmc_mixture_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int i, j, k, kk, snploc, l;
    double skk, sk, clike, ssculm, r;
    bool overflow;

    // gdata%snptracker = 0
    for(k=0; k<gdata->nloci; k++) {
        for(j=0; j<config->ncat; j++) gdata->snptracker[k][j] = 0;
    }

    // Choose annotation
    for (k = 0; k < gdata->nloci; k++) {
        snploc = gdata->permvec[k];
        if (gdata->nannot[snploc] > 1) {
             mstate->gk = mstate->g[snploc];
             // z => gdata%X(:, snploc)
             mstate->zz = gdata->xpx[snploc];
             mstate->zz_vare = mstate->zz / mstate->vare;
             
             // gdata%atemp = 0
             for(j=0; j<config->ncat; j++) gdata->atemp[j] = 0;
             
             if (mstate->rep != 1) {
                 gdata->atemp[gdata->a[snploc]] = 1;
             }
             
             // mstate%dira = dble(gdata%C...) + dble(gdata%atemp)
             for(j=0; j<config->ncat; j++) mstate->dira[j] = (double)gdata->C[snploc][j] + (double)gdata->atemp[j];
             
             // pia = rdirichlet2
             rdirichlet2(rs, config->ncat, mstate->dira, mstate->pia);
             
             if (gdata->vsnptrack[snploc] > 1) {
                 // ytemp = yadj + z * gk
                 for(int row=0; row<gdata->nt; row++) mstate->ytemp[row] = mstate->yadj[row] + gdata->X[row * gdata->nloci + snploc] * mstate->gk;
             } else {
                 for(int row=0; row<gdata->nt; row++) mstate->ytemp[row] = mstate->yadj[row];
             }
             
             mstate->rhs = dot_product_col(gdata->X, snploc, gdata->nt, gdata->nloci, mstate->ytemp);
             // mstate->lhs not used?
             
             mstate->maxs = 0.0;
             for (i = 1; i < config->ndist; i++) {
                 mstate->uhat = mstate->rhs / (mstate->zz + mstate->vare_gp[i]);
                 mstate->maxtemp = 0.5 * mstate->uhat * mstate->rhs / mstate->vare;
                 if (mstate->maxtemp > mstate->maxs) mstate->maxs = mstate->maxtemp;
             }
             
             for (j = 0; j < config->ncat; j++) {
                 if (gdata->C[snploc][j] == 1) {
                     mstate->ss[j] = mstate->p[0][j] * exp(-mstate->maxs); // p(1,j)
                     for(kk=1; kk<config->ndist; kk++) { // p(kk+1, j)
                         mstate->detV = mstate->gp[kk] * mstate->zz_vare + 1.0;
                         mstate->uhat = mstate->rhs / (mstate->zz + mstate->vare_gp[kk]);
                         mstate->ss[j] += mstate->p[kk][j] * pow(mstate->detV, -0.5) * exp(0.5 * mstate->uhat * mstate->rhs / mstate->vare - mstate->maxs);
                     }
                     mstate->ss[j] = log(mstate->pia[j]) + log(mstate->ss[j]);
                 }
             }
             
             // Stabilize
             for(kk=0; kk<config->ncat; kk++) {
                 if (gdata->C[snploc][kk] == 1) {
                     skk = mstate->ss[kk];
                     sk = 0.0;
                     overflow = false;
                     for (l=0; l<config->ncat; l++) {
                         if (gdata->C[snploc][l] == 1) {
                             if (l==kk) continue;
                             clike = mstate->ss[l] - skk;
                             if (clike > LOG_UPPER_LIMIT) { overflow = true; break; }
                             if (clike < -LOG_UPPER_LIMIT) continue;
                             sk += exp(clike);
                         }
                     }
                     if (overflow) mstate->sstemp[kk] = 0.0;
                     else mstate->sstemp[kk] = 1.0 / (1.0 + sk);
                 } else {
                     mstate->sstemp[kk] = 0.0;
                 }
             }
             
             // Sample annotation
             ssculm = 0.0;
             r = rand_uniform(rs, 0.0, 1.0);
             config->annotflag = 0;
             for (kk=0; kk<config->ncat; kk++) {
                 if (gdata->C[snploc][kk] == 1) {
                     ssculm += mstate->sstemp[kk];
                     if (r < ssculm) {
                         config->annotflag = kk;
                         break;
                     }
                 }
             }
             gdata->a[snploc] = config->annotflag;
        }
    }
    
    // Sample effect
    for (k = 0; k < gdata->nloci; k++) {
        snploc = gdata->permvec[k];
        j = gdata->a[snploc]; // annotation index
        
        mstate->gk = mstate->g[snploc];
        mstate->zz = gdata->xpx[snploc];
        mstate->zz_vare = mstate->zz / mstate->vare;
        
        if (gdata->vsnptrack[snploc] > 1) {
             add_col_scalar(mstate->yadj, gdata->X, snploc, gdata->nt, gdata->nloci, mstate->gk);
        }
        
        mstate->rhs = dot_product_col(gdata->X, snploc, gdata->nt, gdata->nloci, mstate->yadj);
        
        mstate->s[0] = mstate->log_p[0][j];
        for(kk=1; kk<config->ndist; kk++) {
            mstate->logdetV = log(mstate->gp[kk] * mstate->zz_vare + 1.0);
            mstate->uhat = mstate->rhs / (mstate->zz + mstate->vare_gp[kk]);
            mstate->s[kk] = -0.5 * (mstate->logdetV - (mstate->rhs * mstate->uhat / mstate->vare)) + mstate->log_p[kk][j];
        }
        
        // Stabilize
        for(kk=0; kk<config->ndist; kk++) {
            skk = mstate->s[kk];
            sk = 0.0;
            overflow = false;
            for(l=0; l<config->ndist; l++) {
                if (l==kk) continue;
                clike = mstate->s[l] - skk;
                if (clike > LOG_UPPER_LIMIT) { overflow=true; break; }
                if (clike < -LOG_UPPER_LIMIT) continue;
                sk += exp(clike);
            }
            if (overflow) mstate->stemp[kk] = 0.0;
            else mstate->stemp[kk] = 1.0 / (1.0 + sk);
        }
        
        // Sample dist
        ssculm = 0.0;
        r = rand_uniform(rs, 0.0, 1.0);
        config->indistflag = 0;
        for(kk=0; kk<config->ndist; kk++) {
            ssculm += mstate->stemp[kk];
            if (r < ssculm) {
                config->indistflag = kk;
                break;
            }
        }
        
        gdata->snptracker[snploc][j] = config->indistflag + 1; // Store 1-based?
        gdata->vsnptrack[snploc] = config->indistflag + 1; // 1-based? Fortran: vsnptrack = indistflag (where indistflag is 1-based)
        
        if (config->indistflag == 0) { // Dist 1 in Fortran
            mstate->gk = 0.0;
        } else {
            mstate->v1 = mstate->zz + mstate->vare / mstate->gp[config->indistflag];
            mstate->gk = rand_normal(rs, mstate->rhs / mstate->v1, sqrt(mstate->vare / mstate->v1));
            add_col_scalar(mstate->yadj, gdata->X, snploc, gdata->nt, gdata->nloci, -mstate->gk);
            mstate->included++;
        }
        mstate->g[snploc] = mstate->gk;
        // msize omitted
    }
    
    // Stats update mixture
    for(j=0; j<config->ncat; j++) {
        for(i=0; i<config->ndist; i++) {
             int cnt = 0;
             double sum_g2 = 0.0;
             for(k=0; k<gdata->nloci; k++) {
                 if (gdata->snptracker[k][j] == i+1) {
                     cnt++;
                     // sum(g2) but mask gdata%gannot? 
                     // In mixture: sum(mstate%g * mstate%g, mask=...)
                     // mstate->g is updated.
                     // IMPORTANT: in mixture, g is updated per snp.
                     // In Fortran additive kernel: g is sum(gannot).
                     // In mixture kernel: mstate%g(snploc) = mstate%gk. 
                     // So usage is correct.
                     sum_g2 += mstate->g[k] * mstate->g[k];
                 }
             }
             mstate->snpindist[i][j] = cnt;
             mstate->varindist[i][j] = sum_g2;
        }
    }
    
    // Included
    double included_sum = 0.0;
     for(i=0; i<config->ndist; i++) included_sum += mstate->snpindist[i][0]; // Wait, Fortran included = nloci - sum(snpindist(1,:)). 
     // snpindist(1,:) is "zero effect" count (dist 1).
     // Total loci = sum all dists. 
     // So included = nloci - count_of_dist_1.
     // In C: dist 0 is index 0.
     int count_zero = 0;
     for(j=0; j<config->ncat; j++) count_zero += mstate->snpindist[0][j];
     mstate->included = gdata->nloci - count_zero;
}
