#include "io.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

void get_size(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    char buffer[1024];
    
    gdata->nind = 0;
    fp = fopen(config->phenfil, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->nind++;
        }
        fclose(fp);
    } else {
        printf("Error opening phen file: %s\n", config->phenfil);
        exit(1);
    }

    gdata->nloci = 0;
    fp = fopen(config->bimfil, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->nloci++;
        }
        fclose(fp);
    } else {
        printf("Error opening bim file: %s\n", config->bimfil);
        exit(1);
    }
}

void load_phenos_plink(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    char line[4096];
    
    gdata->trains = (int*)calloc(gdata->nind, sizeof(int));
    gdata->why = (double*)calloc(gdata->nind, sizeof(double));
    
    fp = fopen(config->phenfil, "r");
    if (!fp) exit(1);
    
    for (int i = 0; i < gdata->nind; i++) {
        if (!fgets(line, sizeof(line), fp)) break;
        
        // Parse line to find config->trait_pos (1-based index) column
        // Standard PLINK pheno file: Family ID, Individual ID, Phenotype...
        // Fortran code skips to column trait_pos + 5 ??
        // Fortran: 
        // pos2 = index(str(pos1:), " ") ... 
        // if (n == config%trait_pos + 5) exit
        // Wait, trait_pos usually is 1 for the first phenotype.
        // PLINK FAM structure: FID IID PID MID SEX PHENOTYPE
        // So default trait_pos=1 means the 6th column (PHENOTYPE).
        // If config%trait_pos is 1, Fortran loop looks for n == 6. Correct.
        
        char *ptr = line;
        int col = 0;
        char *token;
        char *saveptr;
        double val = MISSING_VALUE;
        int found = 0;
        
        // Manual tokenization to handle multiple spaces
        char *p = line;
        while (*p && isspace(*p)) p++; // skip leading
        
        int current_col = 0;
        int target_col = 5 + (config->trait_pos - 1); // 0-based index of target column (6th col is index 5)
        
        // Naive parsing assuming space delimited
        // We'll walk the string
        while (*p) {
            if (current_col == target_col) {
                // parse value
                if (strncmp(p, "NA", 2) == 0) {
                    val = MISSING_VALUE;
                } else {
                    char temp[100];
                    int k=0;
                    while (*p && !isspace(*p) && k<99) temp[k++] = *p++;
                    temp[k] = '\0';
                    if (strcmp(temp, "NA") == 0) {
                         val = MISSING_VALUE;
                    } else {
                         val = atof(temp);
                    }
                }
                found = 1;
                break;
            }
            
            // Skip word
            while (*p && !isspace(*p)) p++;
            // Skip space
            while (*p && isspace(*p)) p++;
            current_col++;
        }
        
        if (found && val != MISSING_VALUE) { 
             gdata->trains[i] = 0;
             gdata->why[i] = val;
        } else {
             gdata->trains[i] = 1;
             gdata->why[i] = MISSING_VALUE;
        }
    }
    fclose(fp);
}

void load_snp_binary(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    unsigned char b1, b2, b3;
    double igen[4] = {0.0, 3.0, 1.0, 2.0}; // HomRef, Missing, Het, HomAlt (PLINK binary format)
    
    // allocate X
    // X is (nt x nloci) double
    // Fortran: gdata%X(gdata%nt, gdata%nloci)
    gdata->X = (double*)calloc(gdata->nt * gdata->nloci, sizeof(double));
    
    fp = fopen(config->genfil, "rb");
    if (!fp) exit(1);
    
    // Header
    fread(&b1, 1, 1, fp);
    fread(&b2, 1, 1, fp);
    fread(&b3, 1, 1, fp);
    // Check magic? (0x6c 0x1b 0x01)
    
    // Loop over loci (columns in PLINK binary usually? No, PLINK is SNP-major usually)
    // Fortran: loop j=1..nloci (outer), i=1..nind (inner)
    // PLINK bed (snp-major):
    // For each SNP:
    //   ceil(N/4) bytes
    
    int nbytes = (gdata->nind + 3) / 4;
    unsigned char *buffer = (unsigned char*)malloc(nbytes);
    
    for (int j = 0; j < gdata->nloci; j++) {
        if (fread(buffer, 1, nbytes, fp) != nbytes) break;
        
        int tr = 0; // index in X rows
        for (int i = 0; i < gdata->nind; i++) {
            int byte_idx = i / 4;
            int bit_idx = (i % 4) * 2;
            unsigned char b = buffer[byte_idx];
            unsigned char code = (b >> bit_idx) & 3;
            double val = igen[code];
            
            if (gdata->trains[i] == 0) {
                // X[tr, j]
                gdata->X[tr * gdata->nloci + j] = val;
                tr++;
            }
        }
    }
    
    free(buffer);
    fclose(fp);
}

void init_random_seed_custom(ModelConfig *config, prng_state *rs) {
    uint64_t seed[4];
    if (config->seed1 != 0) {
        for(int i=0; i<4; i++) seed[i] = abs(config->seed1) + i;
    } else {
        time_t t = time(NULL);
        for(int i=0; i<4; i++) seed[i] = t + 37 * i;
    }
    manual_seed(rs, seed);
}

void allocate_data(ModelConfig *config, GenomicData *gdata, MCMCState *mstate, MCMCStorage *mstore) {
    // Determine nt, ntrain etc.
    if (!config->mcmc) {
        for(int i=0; i<gdata->nind; i++) {
            if (gdata->trains[i] == 0) gdata->trains[i] = 3;
            else if (gdata->trains[i] == 1) gdata->trains[i] = 0;
        }
        for(int i=0; i<gdata->nind; i++) {
            if (gdata->trains[i] == 3) gdata->trains[i] = 1;
        }
    }
    
    gdata->nt = 0;
    for(int i=0; i<gdata->nind; i++) if(gdata->trains[i]==0) gdata->nt++;
    
    // Allocations
    gdata->pred = (double*)calloc(gdata->nind, sizeof(double));
    
    mstate->gpin = (double*)calloc(config->ndist, sizeof(double));
    mstate->gp = (double*)calloc(config->ndist, sizeof(double));
    mstate->p = (double**)malloc(config->ndist * sizeof(double*));
    mstate->log_p = (double**)malloc(config->ndist * sizeof(double*));
    mstate->snpindist = (int**)malloc(config->ndist * sizeof(int*));
    mstate->varindist = (double**)malloc(config->ndist * sizeof(double*));
    
    for(int i=0; i<config->ndist; i++) {
        mstate->p[i] = (double*)calloc(config->ncat, sizeof(double));
        mstate->log_p[i] = (double*)calloc(config->ncat, sizeof(double));
        mstate->snpindist[i] = (int*)calloc(config->ncat, sizeof(int));
        mstate->varindist[i] = (double*)calloc(config->ncat, sizeof(double));
    }

    gdata->permannot = (int*)calloc(config->ncat, sizeof(int));    
    mstate->delta = (double*)calloc(config->ndist, sizeof(double));
    mstate->dirx = (double*)calloc(config->ndist, sizeof(double));
    mstate->g = (double*)calloc(gdata->nloci, sizeof(double));
    mstate->yadj = (double*)calloc(gdata->nt, sizeof(double));
    mstate->z = (double*)calloc(gdata->nt, sizeof(double));
    
    mstate->s = (double*)calloc(config->ndist, sizeof(double));
    mstate->stemp = (double*)calloc(config->ndist, sizeof(double));
    mstate->sstemp = (double*)calloc(config->ncat, sizeof(double));
    
    gdata->xpx = (double*)calloc(gdata->nloci, sizeof(double));
    
    mstore->gstore = (double*)calloc(gdata->nloci, sizeof(double));
    
    // Rank 2 stores
    mstore->snpstore = (double**)malloc(config->ndist * sizeof(double*));
    mstore->varstore = (double**)malloc(config->ndist * sizeof(double*));
    mstore->pstore = (double**)malloc(config->ndist * sizeof(double*));
    for(int i=0; i<config->ndist; i++) {
        mstore->snpstore[i] = (double*)calloc(config->ncat, sizeof(double));
        mstore->varstore[i] = (double*)calloc(config->ncat, sizeof(double));
        mstore->pstore[i] = (double*)calloc(config->ncat, sizeof(double));
    }
    
    // mstore%indiststore(gdata%nloci, config%ndist) - swapped indices?
    // Fortran: (nloci, ndist). C: [nloci][ndist]
    mstore->indiststore = (double**)malloc(gdata->nloci * sizeof(double*));
    for(int i=0; i<gdata->nloci; i++) mstore->indiststore[i] = (double*)calloc(config->ndist, sizeof(double));
    
    mstore->mu_vare_store = (double*)calloc(4, sizeof(double));
    gdata->freqstore = (double*)calloc(gdata->nloci, sizeof(double));
    gdata->permvec = (int*)calloc(gdata->nloci, sizeof(int));
    
    gdata->snptracker = (int**)malloc(gdata->nloci * sizeof(int*));
    for(int i=0; i<gdata->nloci; i++) gdata->snptracker[i] = (int*)calloc(config->ncat, sizeof(int));
    
    mstate->log_gp = (double*)calloc(config->ndist, sizeof(double));
    mstate->vare_gp = (double*)calloc(config->ndist, sizeof(double));
    
    mstore->varustore = (double*)calloc(gdata->nloci, sizeof(double));
    mstore->varistore = (double*)calloc(gdata->nloci, sizeof(double));
    
    gdata->vsnptrack = (int*)calloc(gdata->nloci, sizeof(int));
    gdata->nannot = (int*)calloc(gdata->nloci, sizeof(int)); // nannot is dim(nloci)
    gdata->a = (int*)calloc(gdata->nloci, sizeof(int));
    
    mstate->ss = (double*)calloc(config->ncat, sizeof(double));
    gdata->includedloci = (double*)calloc(gdata->nloci, sizeof(double));
    mstate->pia = (double*)calloc(config->ncat, sizeof(double));
    mstate->ytemp = (double*)calloc(gdata->nloci, sizeof(double));
    mstate->dira = (double*)calloc(config->ncat, sizeof(double));
    gdata->atemp = (int*)calloc(config->ncat, sizeof(int));
    
    mstore->annotstore = (double**)malloc(gdata->nloci * sizeof(double*));
    gdata->gannot = (double**)malloc(gdata->nloci * sizeof(double*));
    for(int i=0; i<gdata->nloci; i++) {
        mstore->annotstore[i] = (double*)calloc(config->ncat, sizeof(double));
        gdata->gannot[i] = (double*)calloc(config->ncat, sizeof(double));
    }
}

void load_param(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore, MCMCState *mstate) {
    FILE *fp;
    char buffer[1024];
    int nc = config->ndist + 1;
    double *gtemp = (double*)calloc(nc, sizeof(double));
    
    fp = fopen(config->paramfil, "r");
    if(!fp) return;
    fgets(buffer, sizeof(buffer), fp); // Header
    
    for(int i=0; i<gdata->nloci; i++) {
        for(int k=0; k<nc; k++) {
            fscanf(fp, "%lf", &gtemp[k]);
        }
        mstore->gstore[i] = gtemp[nc-1];
    }
    fclose(fp);
    free(gtemp);
    
    fp = fopen(config->modfil, "r");
    if(fp) {
        char dum[100];
        fscanf(fp, "%s %lf", dum, &mstate->mu);
        fclose(fp);
    }
}

void load_categories(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    fp = fopen(config->catRC, "r");
    if (!fp) {
        printf("Cannot open cat file %s\n", config->catRC);
        exit(1);
    }
    
    gdata->C = (int**)malloc(gdata->nloci * sizeof(int*));
    for(int i=0; i<gdata->nloci; i++) {
        gdata->C[i] = (int*)calloc(config->ncat + 1, sizeof(int));
        for(int j=1; j<=config->ncat; j++) { // Using 1-based index in struct for C? 
            // Fortran: C(i, j), j=1..ncat.
            // C: C[i][j-1]? 
            // C[i] is size ncat+1. Let's use 0-based for internal C array but loop matches file.
            // "read(unit_cat, *) (gdata%C(i, j), j = 1, config%ncat)"
            // Fortran reads columns.
            int val;
            fscanf(fp, "%d", &val);
            gdata->C[i][j-1] = val; // Store at 0..ncat-1
        }
    }
    fclose(fp);
}

void write_dgv(ModelConfig *config, GenomicData *gdata) {
    FILE *fp = fopen(config->mbvfil, "w");
    if(!fp) return;
    
    for(int i=0; i<gdata->nind; i++) {
        if(gdata->trains[i] == 0) {
            fprintf(fp, "%15.7E\n", gdata->pred[i]);
        } else {
            fprintf(fp, "NA\n");
        }
    }
    fclose(fp);
}

void output_model(ModelConfig *config, GenomicData *gdata, MCMCStorage *mstore) {
    FILE *fp;
    
    fp = fopen(config->paramfil, "w");
    if (!fp) exit(1);
    
    // Header
    fprintf(fp, "  ");
    for(int i=1; i<=config->ndist; i++) fprintf(fp, " PIP%-4d", i);
    fprintf(fp, "  beta");
    for(int i=1; i<=config->ncat; i++) fprintf(fp, " PAIP%-4d", i);
    fprintf(fp, " Vbeta   Vi\n");
    
    for (int i=0; i<gdata->nloci; i++) {
        for(int k=0; k<config->ndist; k++) fprintf(fp, "%15.7E ", mstore->indiststore[i][k]);
        fprintf(fp, "%15.7E ", mstore->gstore[i]);
        for(int k=0; k<config->ncat; k++) fprintf(fp, "%15.7E ", mstore->annotstore[i][k]);
        fprintf(fp, "%15.7E %15.7E\n", mstore->varustore[i], mstore->varistore[i]);
    }
    fclose(fp);
    
    fp = fopen(config->modfil, "w");
    if(!fp) exit(1);
    
    fprintf(fp, "Mean      %15.7E\n", mstore->mu_vare_store[0]);
    fprintf(fp, "Nsnp      %15.7E\n", mstore->mu_vare_store[1]);
    fprintf(fp, "Va        %15.7E\n", mstore->mu_vare_store[2]);
    fprintf(fp, "Ve        %15.7E\n", mstore->mu_vare_store[3]);
    
    // Loops matching Fortran
    for(int j=1; j<=config->ncat; j++) {
        for(int i=1; i<=config->ndist; i++) {
             fprintf(fp, "Nk%d_%d      %15.7E\n", i, j, mstore->snpstore[i-1][j-1]);
        }
    }
    for(int j=1; j<=config->ncat; j++) {
        for(int i=1; i<=config->ndist; i++) {
             fprintf(fp, "Pk%d_%d      %15.7E\n", i, j, mstore->pstore[i-1][j-1]);
        }
    }
    for(int j=1; j<=config->ncat; j++) {
        for(int i=1; i<=config->ndist; i++) {
             fprintf(fp, "Vk%d_%d      %15.7E\n", i, j, mstore->varstore[i-1][j-1]);
        }
    }
    fclose(fp);
}

void output_beta(ModelConfig *config, MCMCState *mstate, GenomicData *gdata) {
    if (!config->fp_beta) return;
    for(int i=0; i<gdata->nloci; i++) {
        fprintf(config->fp_beta, " %15.6E", mstate->g[i]*mstate->g[i]);
    }
    fprintf(config->fp_beta, "\n");
}
