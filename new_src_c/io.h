#ifndef IO_H
#define IO_H

#include "bayesRCO.h"

void get_size(ModelConfig *config, GenomicData *gdata);
void load_phenos_plink(ModelConfig *config, GenomicData *gdata);
void load_snp_binary(ModelConfig *config, GenomicData *gdata);
void init_random_seed_custom(ModelConfig *config, prng_state *rs); // Renamed to avoid confusion
void allocate_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore);
void load_param(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore, MCMCState *mstate);
void load_categories(ModelConfig *config, GenomicData *gdata);
void write_dgv(ModelConfig *config, GenomicData *gdata);
void output_model(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore);
void output_beta(ModelConfig *config, MCMCState *mstate, GenomicData *gdata);

#endif // IO_H
