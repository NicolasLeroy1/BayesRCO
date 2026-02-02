#include "mcmc.h"
#include "mcmc_utils.h"
#include "mcmc_mixture.h"
#include "mcmc_additive.h"
#include "mcmc_bayesCpi.h"
#include "utils.h"
#include "io.h"

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
