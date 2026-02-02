#include "mcmc_additive.h"
#include "mcmc_utils.h"
#include <math.h>
#include <stdbool.h>

void mcmc_additive_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int i, j, k, kk, snploc, l, kcat;
    double skk, sk, clike, ssculm, r;
    bool overflow;
    
    for (k = 0; k < gdata->nloci; k++) {
        snploc = gdata->permvec[k];
        mstate->gk = mstate->g[snploc];
        
        mstate->zz = gdata->xpx[snploc];
        mstate->zz_vare = mstate->zz / mstate->vare;
        
        if (gdata->vsnptrack[snploc] > 1) { // >1 means not index 0 in Fortran (0=no effect, 1=no effect/base? 
            // Fortran: 1..ndist. 1 is usually the "zero effect" distribution or spike.
            // BayesR usually has dist 1 as zero effect.
            // If vsnptrack > 1, means it has an effect.
            // Un-correct yadj
            add_col_scalar(mstate->yadj, gdata->X, snploc, gdata->nt, gdata->nloci, mstate->gk);
        }
        
        mstate->rhs = dot_product_col(gdata->X, snploc, gdata->nt, gdata->nloci, mstate->yadj);
        
        for (kcat = 0; kcat < config->ncat; kcat++) {
            j = gdata->permannot[kcat] - 1; // 0-based
            // Update log_p (already done in update_hypers, but additive kernel refreshes it?)
            // Fortran does: mstate%log_p(1, j) = dlog(mstate%p(1, j)) ...
            // This is just calculating log locally? Okay.
            for(i=0; i<config->ndist; i++) mstate->log_p[i][j] = log(mstate->p[i][j]);
            
             if (gdata->C[snploc][j] == 1) { // C is (nloci x ncat)
                 // lhs
                 // mstate->lhs = mstate->zz / mstate->vare; // Not used below?
                 
                 mstate->s[0] = mstate->log_p[0][j]; // s(1) in Fortran
                 for (kk = 1; kk < config->ndist; kk++) {
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
                         if (clike > LOG_UPPER_LIMIT) { overflow = true; break; }
                         if (clike < -LOG_UPPER_LIMIT) continue;
                         sk += exp(clike);
                     }
                     if (overflow) mstate->stemp[kk] = 0.0;
                     else mstate->stemp[kk] = 1.0 / (1.0 + sk);
                 }
                 
                 // Sample
                 ssculm = 0.0;
                 r = rand_uniform(rs, 0.0, 1.0); // Fortran random_number(r) is [0,1)
                 config->indistflag = 0; // 0-based index for dist
                 for(kk=0; kk<config->ndist; kk++) {
                     ssculm += mstate->stemp[kk];
                     if (r < ssculm) {
                         config->indistflag = kk;
                         break;
                     }
                 }
                 
                 // Need to map back to 1-based for storage if needed, or stick to 0-based.
                 // gdata->snptracker is int. Using 1-based to match Fortran might be safer for comparisons.
                 // Let's store 1-based in snptracker to avoid confusion with 0 initialization.
                 int dist_idx_1based = config->indistflag + 1;
                 gdata->snptracker[snploc][j] = dist_idx_1based;
                 
                 if (dist_idx_1based == 1) {
                     mstate->gk = 0.0;
                 } else {
                     mstate->v1 = mstate->zz + mstate->vare / mstate->gp[config->indistflag];
                     mstate->gk = rand_normal(rs, mstate->rhs / mstate->v1, sqrt(mstate->vare / mstate->v1));
                     add_col_scalar(mstate->yadj, gdata->X, snploc, gdata->nt, gdata->nloci, -mstate->gk);
                 }
                 gdata->gannot[snploc][j] = mstate->gk;
                 
                 // msize logic omitted
             }
        }
    }
    
    // Sum effects
    for(i=0; i<gdata->nloci; i++) {
        double sum_g = 0.0;
        int max_track = 0;
        for(j=0; j<config->ncat; j++) {
            sum_g += gdata->gannot[i][j];
            if (gdata->snptracker[i][j] > max_track) max_track = gdata->snptracker[i][j];
        }
        mstate->g[i] = sum_g;
        gdata->vsnptrack[i] = max_track;
    }
    
    // Update stats
    for(j=0; j<config->ncat; j++) {
        for(i=0; i<config->ndist; i++) { // i is 0-based (0..ndist-1) -> corresponds to 1..ndist
             int count = 0;
             double sum_sq = 0.0;
             for(int l=0; l<gdata->nloci; l++) {
                 if(gdata->snptracker[l][j] == i+1) {
                     count++;
                     sum_sq += gdata->gannot[l][j] * gdata->gannot[l][j];
                 }
             }
             mstate->snpindist[i][j] = count;
             mstate->varindist[i][j] = sum_sq;
        }
    }
    
    // Included
    double included_sum = 0.0;
    for(i=0; i<gdata->nloci; i++) {
        if (gdata->vsnptrack[i] > 1) {
            gdata->includedloci[i] = 1.0;
            included_sum += 1.0;
        } else {
            gdata->includedloci[i] = 0.0;
        }
    }
    mstate->included = (int)included_sum;
}
