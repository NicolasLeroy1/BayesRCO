#ifndef BAYES_STRUCTS_H
#define BAYES_STRUCTS_H

#include <stdbool.h>

/**
 * Configuration parameters for the BayesRCO run.
 */
typedef struct {
    char inprefix[256];
    char outprefix[256];
    char catRC[256];
    
    int numit;
    int burnin;
    int thin;
    int seed;
    int trait_pos;
    int ndist;
    int ncat;
    
    double vara;
    double vare;
    double dfvara;
    double dfvare;
    
    double *gpin; // [ndist]
    double *delta; // [ndist] Dirichlet prior
} BayesConfig;

/**
 * Loaded data (Genotypes, Phenotypes, Annotations).
 */
typedef struct {
    int nind;
    int nloci;
    int nt; // training size (nind where trains[i] == 0)
    
    double *why;        // Phenotypes [nind]
    int *trains;        // Training/Test indicator [nind] (0=train, 1=test)
    double *X;          // Genotype matrix [nloci * nt] (column-major-like flat array: snp-major)
                        // X[snp_idx * nt + ind_idx]
    
    int *C;             // Categories [nloci * ncat] (X[snp_idx * ncat + cat_idx])
    double *freq;       // Allele frequencies [nloci]
} BayesData;

/**
 * MCMC state and model parameters.
 */
typedef struct {
    double mu;
    double vara;
    double vare;
    double vary;
    
    double *g;          // SNP effects [nloci]
    double *p;          // Mix probabilities [ndist * ncat] (p[dist_idx * ncat + cat_idx])
    double *gp;         // [ndist]
    double *yadj;       // Residuals [nt]
    double *xpx;        // SNP diagonals [nloci]
    
    int *snptracker;    // [nloci * ncat]
    int *mc;            // [ndist * ncat]
    double *gh;         // [ndist * ncat]
} BayesModel;

#endif // BAYES_STRUCTS_H
