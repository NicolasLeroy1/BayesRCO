#include "plink_io.h"
#include <string.h>
#include <ctype.h>

int get_size(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    char buffer[1024];
    
    gdata->num_individuals = 0;
    fp = fopen(config->phenotype_file_path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->num_individuals++;
        }
        fclose(fp);
    } else {
        printf("Error opening phen file: %s\n", config->phenotype_file_path);
        return ERR_FILE_IO;
    }

    gdata->num_loci = 0;
    fp = fopen(config->bim_file_path, "r");
    if (fp) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            gdata->num_loci++;
        }
        fclose(fp);
    } else {
        printf("Error opening bim file: %s\n", config->bim_file_path);
        return ERR_FILE_IO;
    }
    return SUCCESS;
}

int load_phenos_plink(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    char line[4096];
    
    gdata->trains = (int*)calloc(gdata->num_individuals, sizeof(int));
    if (!gdata->trains) return ERR_MEMORY;
    gdata->phenotypes = (double*)calloc(gdata->num_individuals, sizeof(double));
    if (!gdata->phenotypes) return ERR_MEMORY;
    
    fp = fopen(config->phenotype_file_path, "r");
    if (!fp) return ERR_FILE_IO;
    
    for (int i = 0; i < gdata->num_individuals; i++) {
        if (!fgets(line, sizeof(line), fp)) break;
        
        char *p = line;
        while (*p && isspace(*p)) p++; // skip leading
        
        int current_col = 0;
        int target_col = 5 + (config->trait_column_index - 1); // 0-based index of target column
        double val = MISSING_VALUE;
        int found = 0;
        
        while (*p) {
            if (current_col == target_col) {
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
            while (*p && !isspace(*p)) p++;
            while (*p && isspace(*p)) p++;
            current_col++;
        }
        
        if (found && val != MISSING_VALUE) { 
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

int load_snp_binary(ModelConfig *config, GenomicData *gdata) {
    FILE *fp;
    unsigned char b1, b2, b3;
    double igen[4] = {0.0, 3.0, 1.0, 2.0}; // HomRef, Missing, Het, HomAlt (PLINK binary format)
    
    gdata->genotypes = (double*)calloc(gdata->num_loci * gdata->num_phenotyped_individuals, sizeof(double));
    if (!gdata->genotypes) return ERR_MEMORY;
    
    fp = fopen(config->genotype_file_path, "rb");
    if (!fp) return ERR_FILE_IO;
    
    if (fread(&b1, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    if (fread(&b2, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    if (fread(&b3, 1, 1, fp) != 1) { fclose(fp); return ERR_FILE_IO; }
    
    int nbytes = (gdata->num_individuals + 3) / 4;
    unsigned char *buffer = (unsigned char*)malloc(nbytes);
    if (!buffer) { fclose(fp); return ERR_MEMORY; }
    
    for (int j = 0; j < gdata->num_loci; j++) {
        if (fread(buffer, 1, nbytes, fp) != nbytes) break;
        
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
