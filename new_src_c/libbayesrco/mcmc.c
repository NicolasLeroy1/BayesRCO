#include "mcmc.h"
#include "mcmc_utils.h"
#include "mcmc_mixture.h"
#include "mcmc_additive.h"
#include "mcmc_bayesCpi.h"
#include "utils.h"
#include "mcmc_core_io.h"

// -------------------------------------------------------------------------
// Main Driver
// -------------------------------------------------------------------------

int run_mcmc(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore, prng_state *rs) {
    int i, j, jj;
    int nloci = gdata->num_loci;
    int ncat = config->num_categories;
    int ndist = config->num_distributions;
    
    mcmc_init_common(config, gdata, mstore);
    
    if (config->nobayesCpi) {
        if (config->mixture) {
            mcmc_mixture_init(config, gdata, mstate);
        } else {
            // mcmc_additive_init
            for(int k=0; k<nloci; k++) {
                 for(int cat=0; cat<ncat; cat++) gdata->effects_per_category[IDX2(k, cat, ncat)] = 0.0;
            }
            mcmc_bayesCpi_init(config, gdata, mstate, mstore); 
        }
    } else {
        mcmc_bayesCpi_init(config, gdata, mstate, mstore);
    }
    
    mcmc_start_values_common(config, gdata, mstate, rs);
    
    if (config->nobayesCpi && !config->mixture) {
         for(j=0; j<ncat; j++) gdata->permannot[j] = j; 
    }
    
    for (i = 1; i <= config->num_iterations; i++) {
        mstate->rep = i;
        mcmc_iteration_pre_common(config, mstate, rs);
        
        if (config->permute) {
            permutate(rs, nloci, gdata->permvec);
        }
        
        if (config->nobayesCpi && !config->mixture) {
            permutate(rs, ncat, gdata->permannot);
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
        
        mcmc_update_hypers_common(ndist, config, gdata, mstate, rs);
        
        if (mstate->rep % config->thinning_interval == 0) {
            if (mstate->rep > config->burnin_iterations) {
                if (config->nobayesCpi && config->mixture) {
                    for(j=0; j<nloci; j++) {
                        jj = gdata->current_distribution[j]; // 1-based stored? yes
                        if (jj > 0) mstore->sum_distribution_counts[IDX2(j, jj-1, ndist)] += 1.0;
                        
                        jj = gdata->current_category[j]; // 0-based index of category
                        mstore->sum_category_counts[IDX2(j, jj, ncat)] += 1.0;
                    }
                } else if (config->nobayesCpi && !config->mixture) {
                    for(j=0; j<nloci; j++) {
                        for(jj=0; jj<ncat; jj++) {
                             if (gdata->distribution_per_category[IDX2(j, jj, ncat)] > 0) {
                                  int idx = gdata->distribution_per_category[IDX2(j, jj, ncat)] - 1;
                                  mstore->sum_distribution_counts[IDX2(j, idx, ndist)] += 1.0;
                             }
                        }
                    }
                }
                
                mcmc_save_samples_common(config, gdata, mstate, mstore);
            }
        }
    }
    
    // Post process
    mstate->counter = (mstate->counter > 0) ? mstate->counter : 1;
    double div = (double)mstate->counter;
    for(int k=0; k<nloci; k++) {
        mstore->sum_snp_effects[k] /= div;
        mstore->varustore[k] /= div;
        mstore->varistore[k] /= div;
        for(int n=0; n<ndist; n++) mstore->sum_distribution_counts[IDX2(k, n, ndist)] /= div;
        for(int n=0; n<ncat; n++) mstore->sum_category_counts[IDX2(k, n, ncat)] /= div;
    }
    
    mstore->mu_vare_store[0] /= div;
    mstore->mu_vare_store[1] /= div;
    mstore->mu_vare_store[2] /= div;
    mstore->mu_vare_store[3] /= div;
    
    for(int n=0; n<ndist; n++) {
        for(int c=0; c<ncat; c++) {
             mstore->sum_mixture_proportions[IDX2(n, c, ncat)] /= div;
             mstore->sum_snps_per_distribution[IDX2(n, c, ncat)] /= div;
             mstore->sum_variance_per_distribution[IDX2(n, c, ncat)] /= div;
        }
    }
    
    if (config->nobayesCpi && !config->mixture) {
         for(int k=0; k<nloci; k++) {
             if (gdata->annotations_per_locus[k] > 1) {
                 for(int n=0; n<ndist; n++) mstore->sum_distribution_counts[IDX2(k, n, ndist)] /= (double)gdata->annotations_per_locus[k];
             }
         }
    }

    output_model(config, gdata, mstore);
    mstate->mu = mstore->mu_vare_store[0];
    compute_dgv(gdata, mstate, mstore);
    write_dgv(config, gdata);

    return SUCCESS;
}
