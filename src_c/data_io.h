#ifndef DATA_IO_H
#define DATA_IO_H

#include "bayes_structs.h"

int get_data_size(const char *phenfil, const char *bimfil, int *nind, int *nloci);
int load_phenotypes(const char *phenfil, int nind, int trait_pos, double *why, int *trains, int *nt);
int load_genotypes(const char *genfil, int nind, int nloci, const int *trains, int nt, double *X);
int load_categories(const char *catRC, int nloci, int ncat, int *C);
void compute_frequencies_and_center(int nt, int nloci, double *X, double *freq);
void free_bayes_data(BayesData *data);

#endif // DATA_IO_H
