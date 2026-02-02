#include "utils.h"
#include <stdlib.h>
#include <math.h>

// Helper for dot product
double dot_product_row(double *X_row, double *vec, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += X_row[i] * vec[i];
    }
    return sum;
}

// From mod_stats.f90

void permutate(prng_state *rs, int n, int *p) {
    int i, j, k, ipj, itemp, m;
    double u[100];
    
    for (i = 0; i < n; i++) {
        p[i] = i; // 0-based permutation vs 1-based in Fortran? 
                  // Fortran: p(i) = i (1..n)
                  // Let's use 0-based p[i] = i (0..n-1) if usage is for array indices. 
                  // If usage is logical ID, we might need check. 
                  // Usage in mod_mcmc: usually for shuffling order of loci.
                  // Assume 0-based for C.
    }
    
    // generate up to 100 u(0,1) numbers at a time.
    for (i = 0; i < n; i += 100) {
        m = (n - i < 100) ? (n - i) : 100;
        for (int x = 0; x < 100; x++) u[x] = rand_uniform(rs, 0.0, 1.0); 
// Assuming rand_uniform takes state
        
        for (j = 0; j < m; j++) {
            ipj = i + j; // 0-based index
            // k = int(u(j) * (n - ipj + 1)) + ipj (Fortran) -> u is 0..1
            // Fortran: u(j) * (n - indices remaining) + current_index
            // C: k = (int)(u[j] * (n - ipj)) + ipj; 
            // Wait, range is [ipj, n-1]. Size is n - ipj.
            k = (int)(u[j] * (n - ipj)) + ipj;
            if (k >= n) k = n - 1; 

            itemp = p[ipj];
            p[ipj] = p[k];
            p[k] = itemp;
        }
    }
}

void compute_dgv(GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    int i, tr = 0;
    
    // Reset pred
    for(i=0; i<gdata->nind; i++) gdata->pred[i] = MISSING_VALUE;

    for (i = 0; i < gdata->nind; i++) {
        if (gdata->trains[i] == 0) { // Using 0 for training set as per load_phenos
             // Wait, compute_dgv in Fortran: if (gdata%trains(i) == 0) then ...
             // BUT mod_stats says:
             // 38:             if (gdata%trains(i) == 0) then
             // 39:                 tr = tr + 1
             // 40:                 gdata%pred(i) = mstate%mu + dot_product(gdata%X(tr, 1:gdata%nloci), mstore%gstore(1:gdata%nloci))
             
             // This implies X has rows only for training individuals?
             // allocate_data: allocate(gdata%X(gdata%nt, gdata%nloci)
             // Yes, X only stores training genotypes.
             // So tr index matches X rows.
             
             // In C, 0-based tr.
             
             // Note: dot_product of row tr of X with gstore.
             // X is (nt x nloci).
             // Assume row-major X: X[tr * nloci + ...]
             
             double dp = dot_product_row(&gdata->X[tr * gdata->nloci], mstore->gstore, gdata->nloci);
             gdata->pred[i] = mstate->mu + dp;
             
             tr++;
        }
    }
}

void compute_residuals(GenomicData *gdata, MCMCState *mstate) {
    int i, tr = 0;
    
    for (i = 0; i < gdata->nind; i++) {
        if (gdata->trains[i] == 0) {
            double dp = dot_product_row(&gdata->X[tr * gdata->nloci], mstate->g, gdata->nloci);
            mstate->yadj[tr] = gdata->why[i] - dp - mstate->mu;
            tr++;
        }
    }
}

// From mod_standardize.f90

void xcenter(ModelConfig *config, GenomicData *gdata) {
    int j, cnt;
    double q, qtest;
    FILE *fp;
    
    // For temp storage of column data if needed, or iterate
    // X is row-major (nt x nloci) -> accessing column j is strided.
    // X[i * nloci + j]
    
    if (config->mcmc) {
        for (j = 0; j < gdata->nloci; j++) {
            // First pass: Calculate q
            double sum = 0.0;
            cnt = 0;
            for (int tr = 0; tr < gdata->nt; tr++) {
                double val = gdata->X[tr * gdata->nloci + j];
                if (val < GENOTYPE_MISSING_THRESHOLD) {
                    sum += val;
                    cnt++;
                }
            }
            
            if (cnt > 0)
                q = sum / (2.0 * cnt);
            else 
                q = 0.0; // Should not handle?

            if (q == 1.0 || q == 0.0) {
                 for (int tr = 0; tr < gdata->nt; tr++) {
                     gdata->X[tr * gdata->nloci + j] = 0.0;
                 }
            } else {
                 double denom = sqrt(2.0 * q * (1.0 - q));
                 for (int tr = 0; tr < gdata->nt; tr++) {
                     double val = gdata->X[tr * gdata->nloci + j];
                     if (val > 2.0) val = 2.0 * q; // Handle missing replacement
                     gdata->X[tr * gdata->nloci + j] = (val - 2.0 * q) / denom;
                 }
            }
            gdata->freqstore[j] = q;
        }
        
        fp = fopen(config->freqfil, "w");
        if (fp) {
            for (j = 0; j < gdata->nloci; j++) {
                fprintf(fp, "%10.6f\n", gdata->freqstore[j]);
            }
            fclose(fp);
        }
        
    } else {
        fp = fopen(config->freqfil, "r");
        if (fp) {
            for (j = 0; j < gdata->nloci; j++) {
                // read E15.7
                if (fscanf(fp, "%lf", &gdata->freqstore[j]) != 1) break;
            }
            fclose(fp);
        }
        
        for (j = 0; j < gdata->nloci; j++) {
            q = gdata->freqstore[j];
            if (q == 1.0 || q == 0.0) {
                 for (int tr = 0; tr < gdata->nt; tr++) {
                     gdata->X[tr * gdata->nloci + j] = 0.0;
                 }
            } else {
                 // Recalculate qtest for missing replacement logic?
                 double sum = 0.0;
                 cnt = 0;
                 for (int tr = 0; tr < gdata->nt; tr++) {
                     double val = gdata->X[tr * gdata->nloci + j];
                     if (val < GENOTYPE_MISSING_THRESHOLD) {
                         sum += val;
                         cnt++;
                     }
                 }
                 if (cnt > 0) qtest = sum / (2.0 * cnt); else qtest = 0.0;
                 
                 double denom = sqrt(2.0 * q * (1.0 - q));
                 for (int tr = 0; tr < gdata->nt; tr++) {
                     double val = gdata->X[tr * gdata->nloci + j];
                     if (val > 2.0) val = 2.0 * qtest; 
                     gdata->X[tr * gdata->nloci + j] = (val - 2.0 * q) / denom;
                 }
            }
        }
    }
}
