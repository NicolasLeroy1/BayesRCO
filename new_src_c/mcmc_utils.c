#include "mcmc_utils.h"
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
