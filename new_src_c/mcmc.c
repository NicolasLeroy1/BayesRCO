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
            for(int k=0; k<gdata->num_loci; k++) {
                 for(int cat=0; cat<config->num_categories; cat++) gdata->effects_per_category[k][cat] = 0.0;
            }
            mcmc_bayesCpi_init(config, gdata, mstate, mstore); 
        }
    } else {
        mcmc_bayesCpi_init(config, gdata, mstate, mstore);
    }
    
    mcmc_start_values_common(config, gdata, mstate, rs);
    
    if (config->nobayesCpi && !config->mixture) {
         for(j=0; j<config->num_categories; j++) gdata->permannot[j] = j; 
    }
    
    // compute_residuals(gdata, mstate); // Redundant
    
    for (i = 1; i <= config->num_iterations; i++) {
        mstate->rep = i;
        mcmc_iteration_pre_common(config, mstate, rs);
        
        if (config->permute) {
            permutate(rs, gdata->num_loci, gdata->permvec);
        }
        
        if (config->nobayesCpi && !config->mixture) {
            permutate(rs, config->num_categories, gdata->permannot);
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
        
        mcmc_update_hypers_common(config->num_distributions, config, gdata, mstate, rs);
        
        if (mstate->rep % config->thinning_interval == 0) {
            if (mstate->rep > config->burnin_iterations) {
                if (config->nobayesCpi && config->mixture) {
                    for(j=0; j<gdata->num_loci; j++) {
                        jj = gdata->current_distribution[j]; // 1-based stored? yes
                        if (jj > 0) mstore->sum_distribution_counts[j][jj-1] += 1.0;
                        
                        jj = gdata->current_category[j]; // 0-based index of category
                        mstore->sum_category_counts[j][jj] += 1.0;
                    }
                } else if (config->nobayesCpi && !config->mixture) {
                    for(j=0; j<gdata->num_loci; j++) {
                        for(jj=0; jj<config->num_categories; jj++) {
                             if (gdata->distribution_per_category[j][jj] > 0) {
                                 int idx = gdata->distribution_per_category[j][jj] - 1;
                                 mstore->sum_distribution_counts[j][idx] += 1.0;
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
    for(int k=0; k<gdata->num_loci; k++) {
        mstore->sum_snp_effects[k] /= div;
        mstore->varustore[k] /= div;
        mstore->varistore[k] /= div;
        for(int n=0; n<config->num_distributions; n++) mstore->sum_distribution_counts[k][n] /= div;
        for(int n=0; n<config->num_categories; n++) mstore->sum_category_counts[k][n] /= div;
    }
    // ... others
    mstore->mu_vare_store[0] /= div;
    mstore->mu_vare_store[1] /= div;
    mstore->mu_vare_store[2] /= div;
    mstore->mu_vare_store[3] /= div;
    
    for(int n=0; n<config->num_distributions; n++) {
        for(int c=0; c<config->num_categories; c++) {
             mstore->sum_mixture_proportions[n][c] /= div;
             mstore->sum_snps_per_distribution[n][c] /= div;
             mstore->sum_variance_per_distribution[n][c] /= div;
        }
    }
    
    if (config->nobayesCpi && !config->mixture) {
         for(int k=0; k<gdata->num_loci; k++) {
             if (gdata->annotations_per_locus[k] > 1) {
                 for(int n=0; n<config->num_distributions; n++) mstore->sum_distribution_counts[k][n] /= (double)gdata->annotations_per_locus[k];
             }
         }
    }

    output_model(config, gdata, mstore);
    mstate->mu = mstore->mu_vare_store[0];
    compute_dgv(gdata, mstate, mstore);
    write_dgv(config, gdata);
}
