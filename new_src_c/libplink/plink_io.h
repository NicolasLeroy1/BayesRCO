#ifndef PLINK_IO_H
#define PLINK_IO_H

#include "../libbayesrco/bayesRCO.h"
#include <stdlib.h>
#include <stdio.h>

void get_size(ModelConfig *config, GenomicData *gdata);
void load_phenos_plink(ModelConfig *config, GenomicData *gdata);
void load_snp_binary(ModelConfig *config, GenomicData *gdata);

#endif // PLINK_IO_H
