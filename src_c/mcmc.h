#ifndef MCMC_H
#define MCMC_H

#include "bayes_structs.h"

void init_model(const BayesConfig *config, const BayesData *data, BayesModel *model);
void run_mcmc(const BayesConfig *config, const BayesData *data, BayesModel *model);
void free_bayes_model(BayesModel *model);

#endif // MCMC_H
