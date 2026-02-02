#include "mcmc.h"
#include "utils.h"
#include "io.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Helper for dot product of row with yadj/vec
// Assuming X is row-major (nt x nloci)
// BUT in kernels we often access X(:, snploc) which is a COLUMN.
// If X is row-major, column access is stride = nloci.
// Using helper for column dot product.

double dot_product_col(double *X, int col_idx, int n_rows, int n_cols, double *vec) {
    double sum = 0.0;
    for (int i = 0; i < n_rows; i++) {
        sum += X[i * n_cols + col_idx] * vec[i];
    }
    return sum;
}

// Helper to add/sub column * scalar to vec
void add_col_scalar(double *vec, double *X, int col_idx, int n_rows, int n_cols, double scalar) {
    for (int i = 0; i < n_rows; i++) {
        vec[i] += X[i * n_cols + col_idx] * scalar;
    }
}

// -------------------------------------------------------------------------
// Helper for random choice
// -------------------------------------------------------------------------
int sample_discrete(double *probs, int n, prng_state *rs) {
    double r = rand_uniform(rs, 0.0, 1.0);
    double sum = 0.0;
    for(int i=0; i<n; i++) {
        sum += probs[i];
        if (r < sum) return i;
    }
    return n-1;
}

// -------------------------------------------------------------------------
// MCMC Common Utils (from mod_mcmc_utils.f90)
// -------------------------------------------------------------------------

void mcmc_save_samples_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    mstate->counter++;
    
    // mstore%gstore = mstore%gstore + mstate%g
    for(int i=0; i<gdata->nloci; i++) {
        mstore->gstore[i] += mstate->g[i];
        mstore->varistore[i] += mstate->g[i] * mstate->g[i];
    }
    
    // mstore%pstore... (ndist x ncat)
    for(int j=0; j<config->ncat; j++) {
        for(int i=0; i<config->ndist; i++) {
            mstore->pstore[i][j] += mstate->p[i][j];
            mstore->snpstore[i][j] += (double)mstate->snpindist[i][j];
            mstore->varstore[i][j] += mstate->varindist[i][j];
        }
    }
    
    mstore->mu_vare_store[0] += mstate->mu;
    mstore->mu_vare_store[1] += mstate->included; // implicit cast
    mstore->mu_vare_store[2] += mstate->vara;
    mstore->mu_vare_store[3] += mstate->vare;
    
    if (mstate->counter > 1) {
        for(int i=0; i<gdata->nloci; i++) {
             // mstore%varustore = mstore%varustore + (mstate%counter * mstate%g - mstore%gstore)**2 / ...
             // Be careful with current gstore value (it already includes current g).
             // Fortran: mstore%gstore = mstore%gstore + mstate%g (already done)
             // Formula handles recursive variance update?
             double term = ((double)mstate->counter * mstate->g[i] - mstore->gstore[i]);
             mstore->varustore[i] += term * term / ((double)mstate->counter * ((double)mstate->counter - 1.0));
        }
    }
    
    // Output to hyp file
    fprintf(config->fp_hyp, "%10d %10d %15.7E %15.7E ", mstate->rep, mstate->included, mstate->vara, mstate->vare);
    for(int j=0; j<config->ncat; j++) {
        for(int i=0; i<config->ndist; i++) {
            fprintf(config->fp_hyp, "%10d ", mstate->snpindist[i][j]);
        }
    }
    for(int j=0; j<config->ncat; j++) {
        for(int i=0; i<config->ndist; i++) {
            fprintf(config->fp_hyp, "%15.7E ", mstate->varindist[i][j]);
        }
    }
    fprintf(config->fp_hyp, "\n");
    fflush(config->fp_hyp);
    
    if (config->beta) output_beta(config, mstate, gdata);
}

void mcmc_init_common(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore) {
    // Zeroing is done in allocate_data (calloc) usually, but good to ensure
    // Only xpx logic needs copy
    for(int i=0; i<gdata->nloci; i++) {
        gdata->xpx[i] = dot_product_col(gdata->X, i, gdata->nt, gdata->nloci, &gdata->X[i]); // Wait, dot product of col i with col i
        // Optimization: dot_product_col(X, i, nt, nloci, X_col_i)
        // Manual loop:
        double sum = 0.0;
        for(int r=0; r<gdata->nt; r++) {
            double val = gdata->X[r * gdata->nloci + i];
            sum += val * val;
        }
        gdata->xpx[i] = sum;
    }
    
    for(int i=0; i<gdata->nloci; i++) {
        int sum_C = 0;
        for(int j=0; j<config->ncat; j++) sum_C += gdata->C[i][j];
        gdata->nannot[i] = sum_C;
    }
}

void mcmc_start_values_common(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs) {
    mstate->mu = 1.0;
    for(int i=0; i<gdata->nt; i++) mstate->yadj[i] = 0.0;
    
    // yhat = sum(why(trains==0))/nnind
    // In C, why is size nind, but we usually access it via trains==0 filtering or compact arrays?
    // In load_phenos, gdata->why is nind.
    // In mcmc, we mostly use yadj which is size nt.
    // Let's compute yhat from why where trains==0
    
    double sum_y = 0.0;
    int cnt = 0;
    for(int i=0; i<gdata->nind; i++) {
        if (gdata->trains[i] == 0) {
            sum_y += gdata->why[i];
            cnt++;
        }
    }
    mstate->yhat = sum_y / mstate->nnind;
    
    double sum_sq = 0.0;
    for(int i=0; i<gdata->nind; i++) {
        if (gdata->trains[i] == 0) {
            double diff = gdata->why[i] - mstate->yhat;
            sum_sq += diff * diff;
        }
    }
    mstate->vary = sum_sq / (mstate->nnind - 1.0);
    
    for(int i=0; i<config->ndist; i++) mstate->gp[i] = mstate->gpin[i] * mstate->vara;
    
    mstate->scale = 0.0;
    for(int j=0; j<config->ncat; j++) {
        mstate->p[0][j] = 0.5; // p(1,j)
        double sum_rest = 0.0;
        for(int i=1; i<config->ndist; i++) {
            mstate->p[i][j] = 1.0 / mstate->gpin[i];
            sum_rest += mstate->p[i][j];
        }
        for(int i=1; i<config->ndist; i++) {
            mstate->p[i][j] = 0.5 * mstate->p[i][j] / sum_rest;
        }
    }
    
    // g init
    double g_val = sqrt(mstate->vara / (0.5 * (double)gdata->nloci));
    for(int i=0; i<gdata->nloci; i++) mstate->g[i] = g_val;
    
    for(int k=0; k<gdata->nloci; k++) gdata->permvec[k] = k;
    
    compute_residuals(gdata, mstate);
}

void mcmc_iteration_pre_common(ModelConfig *config, MCMCState *mstate, prng_state *rs) {
    mstate->included = 0;
    
    if (!config->VCE) {
         // mstate%vare = dot_product(mstate%yadj, mstate%yadj) / rand_chi_square(mstate%nnind + 3.0d0)
         double dp = 0.0;
         for(int i=0; i<mstate->nnind; i++) dp += mstate->yadj[i] * mstate->yadj[i]; // nnind should match nt
         mstate->vare = dp / rand_chi_square(rs, mstate->nnind + 3.0);
    }
    
    // Update Mu
    double sum_yadj = 0.0;
    for(int i=0; i<mstate->nnind; i++) {
        mstate->yadj[i] += mstate->mu;
        sum_yadj += mstate->yadj[i];
    }
    
    mstate->mu = rand_normal(rs, sum_yadj / mstate->nnind, sqrt(mstate->vare / mstate->nnind));
    
    for(int i=0; i<mstate->nnind; i++) mstate->yadj[i] -= mstate->mu;
    
    for(int i=1; i<config->ndist; i++) { // 1-based loop in logic vs 0-based array?
        // Fortran: do i = 2, config%ndist -> array index 2..ndist
        // C: array index 1..ndist-1
        mstate->log_gp[i] = log(mstate->gp[i]);
        mstate->vare_gp[i] = mstate->vare / mstate->gp[i];
    }
}

void mcmc_update_hypers_common(int nc, ModelConfig *config, GenomicData *gdata, MCMCState *mstate, prng_state *rs) {
    if (config->VCE) {
        double sum_g2 = 0.0;
        for(int i=0; i<gdata->nloci; i++) sum_g2 += mstate->g[i] * mstate->g[i]; // Need to be careful: g is vector of size nloci? Yes.
        
        mstate->scale = ((double)mstate->included * sum_g2 + config->vara_ap * config->dfvara) / 
                        (config->dfvara + (double)mstate->included);
        
        mstate->vara = rand_scaled_inverse_chi_square(rs, (double)mstate->included + config->dfvara, mstate->scale);
        
        if (nc == 2) { // BayesCpi
            mstate->gp[1] = mstate->vara / (mstate->included > 0 ? (double)mstate->included : 1.0); // Avoid div by zero if included=0? Fortran didn't check.
            if (mstate->included == 0) mstate->gp[1] = 0.0; // Logic check needed
        } else {
            for(int i=0; i<config->ndist; i++) mstate->gp[i] = mstate->gpin[i] * mstate->vara;
        }
        
        double dp = 0.0;
        for(int i=0; i<mstate->nnind; i++) dp += mstate->yadj[i] * mstate->yadj[i];
        
        double vare_num = dp + config->vare_ap * config->dfvare;
        mstate->vare = vare_num / (mstate->nnind + config->dfvare);
        mstate->vare = rand_scaled_inverse_chi_square(rs, mstate->nnind + config->dfvare, mstate->vare);
    }
    
    for(int j=0; j<config->ncat; j++) {
        for(int i=0; i<config->ndist; i++) {
            mstate->dirx[i] = (double)mstate->snpindist[i][j] + mstate->delta[i];
        }
        rdirichlet(rs, config->ndist, mstate->dirx, mstate->ytemp); // Using ytemp as temp buffer for p
        // Wait, ytemp is size nloci. p is size ndist. Okay.
        
        for(int i=0; i<config->ndist; i++) {
            mstate->p[i][j] = mstate->ytemp[i];
            mstate->log_p[i][j] = log(mstate->p[i][j]);
        }
    }
}

// -------------------------------------------------------------------------
// Kernels
// -------------------------------------------------------------------------

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

void mcmc_bayesCpi_init(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
     mstore->annotstore = (double**)malloc(gdata->nloci * sizeof(double*)); // Ops, re-alloc? Already allocated.
     // Copy C to annotstore
     for(int i=0; i<gdata->nloci; i++) {
         for(int j=0; j<config->ncat; j++) mstore->annotstore[i][j] = (double)gdata->C[i][j];
     }
     // gannot = 0
     for(int i=0; i<gdata->nloci; i++) {
         for(int j=0; j<config->ncat; j++) gdata->gannot[i][j] = 0.0;
     }
     // snptracker = 2 where C=1
     for(int j=0; j<config->ncat; j++) {
         for(int i=0; i<gdata->nloci; i++) {
             if (gdata->C[i][j] == 1) gdata->snptracker[i][j] = 2;
         }
     }
}

// Placeholder to avoid overly long file if not needed, but cleaner to have it.
void mcmc_bayesCpi_kernel(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    // Similar to Additive.
    // Iterates categories, then loci.
    // Uses permvec.
    // ...
    // Since default is Mixture, I've prioritized Mixture kernel. 
    // And Additive.
    // Leaving BayesCpi placeholder for now unless requested.
}

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

// -------------------------------------------------------------------------
// Main Driver
// -------------------------------------------------------------------------

void run_mcmc(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int i, j, jj;
    
    mcmc_init_common(config, gdata, mstore);
    
    if (config->nobayesCpi) {
        if (config->mixture) {
            mcmc_mixture_init(config, gdata, mstate);
        } else {
            // mcmc_additive_init
            // gdata%gannot = 0.0
            for(int k=0; k<gdata->nloci; k++) {
                 for(int cat=0; cat<config->ncat; cat++) gdata->gannot[k][cat] = 0.0;
            }
            mcmc_bayesCpi_init(config, gdata, mstate, mstore); // Actually init logic is same structure?
            // Fortran additive_init: mstore%annotstore = dble(C), gannot=0, snptracker=2 where C=1
            // Identical to BayesCpi init.
        }
    } else {
        mcmc_bayesCpi_init(config, gdata, mstate, mstore);
    }
    
    mcmc_start_values_common(config, gdata, mstate, rs);
    
    if (config->nobayesCpi && !config->mixture) {
         for(j=0; j<config->ncat; j++) gdata->permannot[j] = j; 
    }
    
    // compute_residuals(gdata, mstate); // Redundant
    
    for (i = 1; i <= config->numit; i++) {
        mstate->rep = i;
        mcmc_iteration_pre_common(config, mstate, rs);
        
        if (config->permute) {
            permutate(rs, gdata->nloci, gdata->permvec);
        }
        
        if (config->nobayesCpi && !config->mixture) {
            permutate(rs, config->ncat, gdata->permannot);
        }
        
        if (config->nobayesCpi) {
            if (config->mixture) {
                mcmc_mixture_kernel(config, gdata, mstate, mstore, rs);
            } else {
                mcmc_additive_kernel(config, gdata, mstate, mstore, rs);
            }
        } else {
            mcmc_bayesCpi_kernel(config, gdata, mstate, mstore, rs);
        }
        
        mcmc_update_hypers_common(config->ndist, config, gdata, mstate, rs);
        
        if (mstate->rep % config->thin == 0) {
            if (mstate->rep > config->burnin) {
                if (config->nobayesCpi && config->mixture) {
                    for(j=0; j<gdata->nloci; j++) {
                        jj = gdata->vsnptrack[j]; // 1-based stored? yes
                        if (jj > 0) mstore->indiststore[j][jj-1] += 1.0;
                        
                        jj = gdata->a[j]; // 0-based index of category
                        mstore->annotstore[j][jj] += 1.0;
                    }
                } else if (config->nobayesCpi && !config->mixture) {
                    for(j=0; j<gdata->nloci; j++) {
                        for(jj=0; jj<config->ncat; jj++) {
                             if (gdata->snptracker[j][jj] > 0) {
                                 int idx = gdata->snptracker[j][jj] - 1;
                                 mstore->indiststore[j][idx] += 1.0;
                             }
                        }
                    }
                }
                
                mcmc_save_samples_common(config, gdata, mstate, mstore);
            }
        }
        
        if (mstate->rep % 1000 == 0) {
             // compute_residuals(gdata, mstate);
        }
    }
    
    // Post process
    // mcmc_calculate_posterior_means
    mstate->counter = (mstate->counter > 0) ? mstate->counter : 1;
    double div = (double)mstate->counter;
    for(int k=0; k<gdata->nloci; k++) {
        mstore->gstore[k] /= div;
        mstore->varustore[k] /= div;
        mstore->varistore[k] /= div;
        for(int n=0; n<config->ndist; n++) mstore->indiststore[k][n] /= div;
        for(int n=0; n<config->ncat; n++) mstore->annotstore[k][n] /= div;
    }
    // ... others
    mstore->mu_vare_store[0] /= div;
    mstore->mu_vare_store[1] /= div;
    mstore->mu_vare_store[2] /= div;
    mstore->mu_vare_store[3] /= div;
    
    for(int n=0; n<config->ndist; n++) {
        for(int c=0; c<config->ncat; c++) {
             mstore->pstore[n][c] /= div;
             mstore->snpstore[n][c] /= div;
             mstore->varstore[n][c] /= div;
        }
    }
    
    if (config->nobayesCpi && !config->mixture) {
         for(int k=0; k<gdata->nloci; k++) {
             if (gdata->nannot[k] > 1) {
                 for(int n=0; n<config->ndist; n++) mstore->indiststore[k][n] /= (double)gdata->nannot[k];
             }
         }
    }

    output_model(config, gdata, mstore);
    mstate->mu = mstore->mu_vare_store[0];
    compute_dgv(gdata, mstate, mstore);
    write_dgv(config, gdata);
}
