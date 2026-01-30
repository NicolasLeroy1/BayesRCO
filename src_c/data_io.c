#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_io.h"

int get_data_size(const char *phenfil, const char *bimfil, int *nind, int *nloci) {
    FILE *fp;
    char buffer[4096];

    *nind = 0;
    fp = fopen(phenfil, "r");
    if (!fp) {
        perror("Error opening phenotype file");
        return -1;
    }
    while (fgets(buffer, sizeof(buffer), fp)) {
        if (strlen(buffer) > 1) (*nind)++;
    }
    fclose(fp);

    *nloci = 0;
    fp = fopen(bimfil, "r");
    if (!fp) {
        perror("Error opening bim file");
        return -1;
    }
    while (fgets(buffer, sizeof(buffer), fp)) {
        if (strlen(buffer) > 1) (*nloci)++;
    }
    fclose(fp);

    return 0;
}

int load_phenotypes(const char *phenfil, int nind, int trait_pos, double *why, int *trains, int *nt) {
    FILE *fp = fopen(phenfil, "r");
    if (!fp) return -1;

    char line[8192];
    *nt = 0;
    int target_col = trait_pos + 5;

    for (int i = 0; i < nind; i++) {
        if (!fgets(line, sizeof(line), fp)) break;

        char *ptr = line;
        char *token;
        int col = 0;
        char *val_str = NULL;
        char *saveptr;

        while ((token = strtok_r(ptr, " \t\n\r", &saveptr))) {
            ptr = NULL; // for subsequent calls
            col++;
            if (col == target_col) {
                val_str = token;
                break;
            }
        }

        if (val_str) {
            if (strcmp(val_str, "NA") == 0) {
                trains[i] = 1;
                why[i] = -999.0;
            } else {
                why[i] = atof(val_str);
                trains[i] = 0;
                (*nt)++;
            }
        } else {
            trains[i] = 1;
            why[i] = -999.0;
        }
    }
    fclose(fp);
    return 0;
}

int load_genotypes(const char *genfil, int nind, int nloci, const int *trains, int nt, double *X) {
    FILE *fp = fopen(genfil, "rb");
    if (!fp) return -1;

    unsigned char magic[3];
    if (fread(magic, 1, 3, fp) != 3 || magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
        fprintf(stderr, "Error: Invalid PLINK bed file (must be snp-major)\n");
        fclose(fp);
        return -1;
    }

    double igen[4] = {0.0, 3.0, 1.0, 2.0}; // 3.0 for missing
    int nbytes = (nind + 3) / 4;
    unsigned char *buf = malloc(nbytes);

    for (int j = 0; j < nloci; j++) {
        if (fread(buf, 1, nbytes, fp) != nbytes) {
            fprintf(stderr, "Error reading genotypes at locus %d\n", j);
            free(buf);
            fclose(fp);
            return -1;
        }

        int tr = 0;
        for (int i = 0; i < nind; i++) {
            if (trains[i] == 0) {
                int byte_idx = i / 4;
                int bit_idx = (i % 4) * 2;
                int b = (buf[byte_idx] >> bit_idx) & 3;
                X[j * nt + tr] = igen[b];
                tr++;
            }
        }
    }

    free(buf);
    fclose(fp);
    return 0;
}

int load_categories(const char *catRC, int nloci, int ncat, int *C) {
    if (!catRC || strlen(catRC) == 0) {
        memset(C, 0, nloci * ncat * sizeof(int));
        return 0;
    }

    FILE *fp = fopen(catRC, "r");
    if (!fp) return -1;

    for (int i = 0; i < nloci; i++) {
        for (int j = 0; j < ncat; j++) {
            if (fscanf(fp, "%d", &C[i * ncat + j]) != 1) {
                C[i * ncat + j] = 0;
            }
        }
    }

    fclose(fp);
    return 0;
}

void compute_frequencies_and_center(int nt, int nloci, double *X, double *freq) {
    for (int j = 0; j < nloci; j++) {
        double sum = 0.0;
        int nomiss = 0;
        double *snp_ptr = &X[j * nt];

        for (int i = 0; i < nt; i++) {
            if (snp_ptr[i] < 3.0) {
                sum += snp_ptr[i];
                nomiss++;
            }
        }

        double q = (nomiss > 0) ? (sum / (2.0 * nomiss)) : 0.5;
        freq[j] = q;

        if (q <= 0.0 || q >= 1.0) {
            for (int i = 0; i < nt; i++) snp_ptr[i] = 0.0;
        } else {
            double mean = 2.0 * q;
            double sd = sqrt(2.0 * q * (1.0 - q));
            for (int i = 0; i < nt; i++) {
                if (snp_ptr[i] > 2.0) snp_ptr[i] = mean;
                snp_ptr[i] = (snp_ptr[i] - mean) / sd;
            }
        }
    }
}

void free_bayes_data(BayesData *data) {
    if (!data) return;
    free(data->why);
    free(data->trains);
    free(data->X);
    free(data->C);
    free(data->freq);
}
