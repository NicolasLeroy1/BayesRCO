/**
 * @file plink_io.c
 * @brief PLINK format file I/O functions.
 * 
 * Provides reading of PLINK binary format files (.bed, .bim, .fam).
 */

#include "bayesrco_io.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* =========================================================================
 * PLINK Size Functions
 * ========================================================================= */

/**
 * Get dimensions from PLINK files (count lines).
 */
int plink_get_size(const char *prefix, int *num_ind, int *num_snps) {
    FILE *fp;
    char buffer[1024];
    char path[PATH_MAX_LENGTH];
    
    /* Count individuals from .fam file */
    *num_ind = 0;
    snprintf(path, sizeof(path), "%s.fam", prefix);
    fp = fopen(path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            (*num_ind)++;
        }
        fclose(fp);
    } else {
        fprintf(stderr, "Error opening fam file: %s\n", path);
        return ERR_FILE_IO;
    }

    /* Count SNPs from .bim file */
    *num_snps = 0;
    snprintf(path, sizeof(path), "%s.bim", prefix);
    fp = fopen(path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            (*num_snps)++;
        }
        fclose(fp);
    } else {
        fprintf(stderr, "Error opening bim file: %s\n", path);
        return ERR_FILE_IO;
    }
    return SUCCESS;
}

/**
 * Load phenotypes from PLINK .fam file.
 */
int plink_load_phenotypes(const char *prefix, double *phenotypes, int num_ind) {
    FILE *fp;
    char line[4096];
    char path[PATH_MAX_LENGTH];
    
    snprintf(path, sizeof(path), "%s.fam", prefix);
    fp = fopen(path, "r");
    if (!fp) return ERR_FILE_IO;
    
    for (int i = 0; i < num_ind; i++) {
        if (!fgets(line, sizeof(line), fp)) break;
        
        /* Parse 6th column (phenotype) from .fam file */
        char *p = line;
        int col = 0;
        while (*p && col < 5) {
            while (*p && !isspace(*p)) p++;
            while (*p && isspace(*p)) p++;
            col++;
        }
        
        if (*p) {
            if (strncmp(p, "NA", 2) == 0 || strncmp(p, "-9", 2) == 0) {
                phenotypes[i] = MISSING_VALUE;
            } else {
                phenotypes[i] = atof(p);
            }
        } else {
            phenotypes[i] = MISSING_VALUE;
        }
    }
    fclose(fp);
    return SUCCESS;
}

/**
 * Load genotypes from PLINK .bed file.
 */
int plink_load_genotypes(const char *prefix, double *genotypes, int num_ind, int num_snps) {
    FILE *fp;
    unsigned char b1, b2, b3;
    double igen[4] = {0.0, 3.0, 1.0, 2.0}; /* HomRef, Missing, Het, HomAlt (PLINK binary format) */
    char path[PATH_MAX_LENGTH];
    
    snprintf(path, sizeof(path), "%s.bed", prefix);
    fp = fopen(path, "rb");
    if (!fp) return ERR_FILE_IO;
    
    /* Read magic bytes */
    if (fread(&b1, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    if (fread(&b2, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    if (fread(&b3, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    
    /* Verify magic number (0x6c, 0x1b for PLINK format) */
    if (b1 != 0x6c || b2 != 0x1b) {
        fclose(fp);
        return ERR_FILE_IO;
    }
    
    int nbytes = (num_ind + 3) / 4;
    unsigned char *buffer = (unsigned char*)malloc(nbytes);
    if (!buffer) { fclose(fp); return ERR_MEMORY; }
    
    for (int j = 0; j < num_snps; j++) {
        if (fread(buffer, 1, nbytes, fp) != (size_t)nbytes) break;
        
        for (int i = 0; i < num_ind; i++) {
            int byte_idx = i / 4;
            int bit_idx = (i % 4) * 2;
            unsigned char b = buffer[byte_idx];
            unsigned char code = (b >> bit_idx) & 3;
            genotypes[j * num_ind + i] = igen[code];
        }
    }
    
    free(buffer);
    fclose(fp);
    return SUCCESS;
}

/* =========================================================================
 * Legacy IOConfig-based Functions (for backwards compatibility)
 * ========================================================================= */

/**
 * Get size using IOConfig paths.
 */
int io_get_size(IOConfig *ioconfig, GenomicData *gdata) {
    FILE *fp;
    char buffer[1024];
    
    gdata->num_individuals = 0;
    fp = fopen(ioconfig->pheno_file_path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->num_individuals++;
        }
        fclose(fp);
    } else {
        fprintf(stderr, "Error opening phen file: %s\n", ioconfig->pheno_file_path);
        return ERR_FILE_IO;
    }

    /* Use .bim path derived from geno path by changing extension */
    char bim_path[PATH_MAX_LENGTH];
    strncpy(bim_path, ioconfig->geno_file_path, PATH_MAX_LENGTH-1);
    char *ext = strrchr(bim_path, '.');
    if (ext) strcpy(ext, ".bim");
    
    gdata->num_loci = 0;
    fp = fopen(bim_path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->num_loci++;
        }
        fclose(fp);
    } else {
        fprintf(stderr, "Error opening bim file: %s\n", bim_path);
        return ERR_FILE_IO;
    }
    return SUCCESS;
}

/**
 * Load phenotypes using IOConfig.
 */
int io_load_phenotypes(IOConfig *ioconfig, GenomicData *gdata) {
    FILE *fp;
    char line[4096];
    
    gdata->trains = (int*)calloc(gdata->num_individuals, sizeof(int));
    if (!gdata->trains) return ERR_MEMORY;
    gdata->phenotypes = (double*)calloc(gdata->num_individuals, sizeof(double));
    if (!gdata->phenotypes) return ERR_MEMORY;
    
    fp = fopen(ioconfig->pheno_file_path, "r");
    if (!fp) return ERR_FILE_IO;
    
    for (int i = 0; i < gdata->num_individuals; i++) {
        if (!fgets(line, sizeof(line), fp)) break;
        
        char *p = line;
        while (*p && isspace(*p)) p++;
        
        /* Use trait column from config (default 1 means 6th column in .fam) */
        int col = 0;
        int target_col = 5 + (ioconfig->trait_column_index - 1);
        double val = MISSING_VALUE;
        
        while (*p && col < target_col) {
            while (*p && !isspace(*p)) p++;
            while (*p && isspace(*p)) p++;
            col++;
        }
        
        if (*p) {
            if (strncmp(p, "NA", 2) == 0 || strncmp(p, "-9", 2) == 0) {
                val = MISSING_VALUE;
            } else {
                val = atof(p);
            }
        }
        
        if (val != MISSING_VALUE) { 
            gdata->trains[i] = 0;
            gdata->phenotypes[i] = val;
        } else {
            gdata->trains[i] = 1;
            gdata->phenotypes[i] = MISSING_VALUE;
        }
    }
    fclose(fp);
    return SUCCESS;
}

/**
 * Load genotypes using IOConfig.
 */
int io_load_genotypes(IOConfig *ioconfig, GenomicData *gdata) {
    FILE *fp;
    unsigned char b1, b2, b3;
    double igen[4] = {0.0, 3.0, 1.0, 2.0};
    
    gdata->genotypes = (double*)calloc(gdata->num_loci * gdata->num_phenotyped_individuals, sizeof(double));
    if (!gdata->genotypes) return ERR_MEMORY;
    
    fp = fopen(ioconfig->geno_file_path, "rb");
    if (!fp) return ERR_FILE_IO;
    
    /* Read magic bytes */
    if (fread(&b1, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    if (fread(&b2, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    if (fread(&b3, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    
    int nbytes = (gdata->num_individuals + 3) / 4;
    unsigned char *buffer = (unsigned char*)malloc(nbytes);
    if (!buffer) { fclose(fp); return ERR_MEMORY; }
    
    for (int j = 0; j < gdata->num_loci; j++) {
        if (fread(buffer, 1, nbytes, fp) != (size_t)nbytes) break;
        
        int tr = 0;
        for (int i = 0; i < gdata->num_individuals; i++) {
            int byte_idx = i / 4;
            int bit_idx = (i % 4) * 2;
            unsigned char b = buffer[byte_idx];
            unsigned char code = (b >> bit_idx) & 3;
            double val = igen[code];
            
            if (gdata->trains[i] == 0) {
                gdata->genotypes[j * gdata->num_phenotyped_individuals + tr] = val;
                tr++;
            }
        }
    }
    
    free(buffer);
    fclose(fp);
    return SUCCESS;
}
