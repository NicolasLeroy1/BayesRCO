#ifndef BAYESRCO_H
#define BAYESRCO_H

#include <stdio.h>
#include <stdbool.h>
#include "rng.h"

// Constants from mod_defs.f90
#define PI 3.141592653589793238462
#define LOG_UPPER_LIMIT 700.0
#define MISSING_VALUE -9999.0
#define GENOTYPE_MISSING_THRESHOLD 3.0

// Data structures from mod_data.f90

typedef struct {
    int numit, burnin, thin, ndist, seed1, trait_pos, ncat, msize, mrep;
    int indistflag, burn, annotflag;
    int unit_log, unit_hyp, unit_loc, unit_cat, unit_beta; // File descriptors/IDs usually, but here we might use FILE* separately or keep IDs if needed for logic
    bool mcmc, snpout, permute, cat, beta, mixture, nobayesCpi, VCE;
    double dfvara, dfvare, vara_ap, vare_ap;
    char genfil[200], phenfil[200], bimfil[200], inprefix[200], outprefix[200];
    char logfil[200], freqfil[200], mbvfil[200], hypfil[200], locfil[200];
    char modfil[200], paramfil[200], betafil[200], catfil[200], catRC[200];
    
    // FILE pointers for better C handling
    FILE *fp_log, *fp_hyp, *fp_loc, *fp_cat, *fp_beta;
} ModelConfig;

typedef struct {
    int nloci, nind, nt, ntrain, ntest, nref;
    double *why;        // allocatable
    double *pred;       // allocatable
    double *freqstore;  // allocatable
    double *xpx;        // allocatable
    double *includedloci; // allocatable
    double *X;          // allocatable rank 2 (flattened 1D array)
                        // actually mod_data defines X as real(dp), allocatable(:,:). 
    
    int **C;            // allocatable rank 2. 'integer, dimension(:,:), allocatable :: C'
    
    int *nannot;        // allocatable
    int *a;             // allocatable
    int *vsnptrack;     // allocatable
    int *trains;        // allocatable
    int *permvec;       // allocatable
    int *permannot;     // allocatable
    int *atemp;         // allocatable
    
    int **snptracker;   // allocatable rank 2
    
    double **gannot;    // allocatable rank 2
    
    double msep, bhat, ahat, corr;
} GenomicData;

typedef struct {
    double mu, vara, vare, scale, yhat, vary, nnind;
    int included, counter, rep;
    
    double *gp;         // allocatable
    double *gpin;       // allocatable
    double *delta;      // allocatable
    double *g;          // allocatable
    double *yadj;       // allocatable
    
    double *dirx;       // allocatable
    double *pia;        // allocatable
    double *ytemp;      // allocatable
    double *dira;       // allocatable
    double *s;          // allocatable
    double *stemp;      // allocatable
    double *sstemp;     // allocatable
    double *ss;         // allocatable
    
    double **varindist; // allocatable rank 2
    double **p;         // allocatable rank 2
    double **log_p;     // allocatable rank 2
    
    int **snpindist;    // allocatable rank 2
    
    double *log_gp;     // allocatable rank 2 in Fortran? No, rank 1 in mod_data: 'real(dp), dimension(:), allocatable :: log_gp'
    double *vare_gp;    // allocatable rank 1
    
    double zz, rhs, lhs, v1, gk, zz_vare;
    double logdetV, uhat, total_ssq, detV, maxs, maxtemp;
    double xhat, sk, skk, r, ssculm, clike;
    
    double *z;          // pointer
} MCMCState;

typedef struct {
    double *gstore;        // allocatable
    double *mu_vare_store; // allocatable
    double *varustore;     // allocatable
    double *varistore;     // allocatable
    
    double **snpstore;     // allocatable rank 2
    double **varstore;     // allocatable rank 2
    double **indiststore;  // allocatable rank 2
    double **pstore;       // allocatable rank 2
    double **annotstore;   // allocatable rank 2
} MCMCStorage;

// Global RNG state (mimicking implicit state or passed state)
// mod_random uses random_number which is intrinsic. We will use a passed state or global. 
// Given the C structure, let's put RNG state in ModelConfig or a global if we want strictly to match 'mod_random' being a module.
// But passed state is cleaner. Let's add it to MCMCState or ModelConfig? 
// Actually, `call init_random_seed(config)` suggests config might hold seed, but mod_random functions don't take state.
// To fully emulate Fortran's `random_number` which is global, we might need a global instance in `rng.c`.
// However, the user provided `libgfortran_rng.c` has a struct `prng_state`.
// Let's create a global instance in `rng.c` and expose functions that use it, OR pass it around.
// Passing it around is better for C, but might require changing function signatures compared to Fortran.
// Fortran: `rand_uniform(a,b)`. C: `rand_uniform(&rng, a, b)`. This is fine.

#endif // BAYESRCO_H
