#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mcmc.h"
#include "stats.h"

static void shuffle(int *array, int n) {
    if (n > 1) {
        for (int i = n - 1; i > 0; i--) {
            int j = rand() % (i + 1);
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

void init_model(const BayesConfig *config, const BayesData *data, BayesModel *model) {
    int nloci = data->nloci;
    int nt = data->nt;
    int ndist = config->ndist;
    int ncat = config->ncat;

    model->mu = 1.0;
    model->vara = config->vara;
    model->vare = config->vare;

    // Estimate vary (variance of why in training)
    double sum_y = 0.0;
    double ssq_y = 0.0;
    for (int i = 0; i < data->nind; i++) {
        if (data->trains[i] == 0) {
            sum_y += data->why[i];
            ssq_y += data->why[i] * data->why[i];
        }
    }
    double mean_y = sum_y / nt;
    model->vary = (ssq_y - sum_y * mean_y) / (nt - 1.0);

    if (config->dfvara < -2.0) {
        model->vara = config->vara * model->vary;
    }

    model->g = calloc(nloci, sizeof(double));
    model->p = calloc(ndist * ncat, sizeof(double));
    model->gp = calloc(ndist, sizeof(double));
    model->yadj = calloc(nt, sizeof(double));
    model->xpx = calloc(nloci, sizeof(double));
    model->snptracker = calloc(nloci * ncat, sizeof(int));
    model->mc = calloc(ndist * ncat, sizeof(int));
    model->gh = calloc(ndist * ncat, sizeof(double));

    // Initialize gp
    for (int i = 0; i < ndist; i++) {
        model->gp[i] = config->gpin[i] * model->vara;
    }

    // Initialize p (probabilities)
    for (int j = 0; j < ncat; j++) {
        model->p[0 * ncat + j] = 0.5;
        double sum_p = 0.0;
        for (int i = 1; i < ndist; i++) {
            model->p[i * ncat + j] = 1.0 / (config->gpin[i] > 0 ? config->gpin[i] : 1.0);
            sum_p += model->p[i * ncat + j];
        }
        for (int i = 1; i < ndist; i++) {
            model->p[i * ncat + j] = 0.5 * model->p[i * ncat + j] / sum_p;
        }
    }

    // Initialize g (SNP effects)
    double g_init = sqrt(model->vara / (0.5 * nloci));
    for (int i = 0; i < nloci; i++) model->g[i] = g_init;

    // Initial residuals: yadj = why - Xg - mu
    for (int i = 0; i < nt; i++) {
        model->yadj[i] = data->why[i] - model->mu;
    }

    // Precompute xpx
    for (int j = 0; j < nloci; j++) {
        model->xpx[j] = dot_product(&data->X[j * nt], &data->X[j * nt], nt);
    }
}

static void update_vara(const BayesConfig *config, BayesModel *model) {
    int total_mc = 0;
    double total_gh = 0.0;
    for (int i = 1; i < config->ndist; i++) {
        for (int j = 0; j < config->ncat; j++) {
            total_mc += model->mc[i * config->ncat + j];
            total_gh += model->gh[i * config->ncat + j];
        }
    }

    double shape = 0.5 * (total_mc + config->dfvara);
    double scale = 0.5 * (config->vara * (config->dfvara > 0 ? config->dfvara : 0.0) + total_gh);

    if (config->dfvara > 0.0) {
        model->vara = rand_inverse_gamma(shape, scale);
    } else {
        if (total_mc > 0) model->vara = rand_inverse_gamma(0.5 * total_mc, 0.5 * total_gh);
    }

    for (int i = 0; i < config->ndist; i++) model->gp[i] = config->gpin[i] * model->vara;
}

static void update_pi(const BayesConfig *config, BayesModel *model) {
    for (int j = 0; j < config->ncat; j++) {
        double *alpha = malloc(config->ndist * sizeof(double));
        double *res = malloc(config->ndist * sizeof(double));
        for (int i = 0; i < config->ndist; i++) {
            alpha[i] = config->delta[i] + model->mc[i * config->ncat + j];
        }
        rdirichlet(config->ndist, alpha, res);
        for (int i = 0; i < config->ndist; i++) {
            model->p[i * config->ncat + j] = res[i];
        }
        free(alpha);
        free(res);
    }
}

void run_mcmc(const BayesConfig *config, const BayesData *data, BayesModel *model) {
    int nloci = data->nloci;
    int nt = data->nt;
    int ndist = config->ndist;
    int ncat = config->ncat;

    int *permvec = malloc(nloci * sizeof(int));
    for (int i = 0; i < nloci; i++) permvec[i] = i;

    double *s = malloc(ndist * sizeof(double));
    double *stemp = malloc(ndist * sizeof(double));
    double *vare_gp = malloc(ndist * sizeof(double));

    char hyp_file[512];
    sprintf(hyp_file, "%s.hyp", config->outprefix);
    FILE *fp_hyp = fopen(hyp_file, "w");
    if (fp_hyp) fprintf(fp_hyp, "Iteration mu vara vare\n");

    for (int rep = 1; rep <= config->numit; rep++) {
        shuffle(permvec, nloci);
        int included = 0;
        memset(model->mc, 0, ndist * ncat * sizeof(int));
        memset(model->gh, 0, ndist * ncat * sizeof(double));

        // 1. Update vare
        double ssq_res = dot_product(model->yadj, model->yadj, nt);
        model->vare = ssq_res / rand_chi_square((double)nt + 3.0);

        // 2. Update mu
        for (int i = 0; i < nt; i++) model->yadj[i] += model->mu;
        double sum_yadj = 0.0;
        for (int i = 0; i < nt; i++) sum_yadj += model->yadj[i];
        
        model->mu = rand_normal(sum_yadj / nt, sqrt(model->vare / nt));
        for (int i = 0; i < nt; i++) model->yadj[i] -= model->mu;

        // Update vare_gp ratios
        for (int i = 1; i < ndist; i++) {
            vare_gp[i] = model->vare / (model->gp[i] > 0 ? model->gp[i] : 1e-10);
        }

        // 3. Update SNP effects
        for (int ii = 0; ii < nloci; ii++) {
            int k = permvec[ii];
            double *z = &data->X[k * nt];
            double zz = model->xpx[k];
            double gk = model->g[k];

            // Compute rhs without modifying yadj first
            double rhs = dot_product(model->yadj, z, nt) + zz * gk;

            // Identify active categories for this SNP
            int active_cats[100]; 
            int n_active = 0;
            for(int j=0; j<ncat; j++) {
                if(data->C[k*ncat + j]) {
                    if (n_active < 100) active_cats[n_active++] = j;
                }
            }
            if (n_active == 0) { // Fallback if no category assigned
                active_cats[0] = 0; n_active = 1;
            }

            // Calculate joint probabilities P(A=a, G=d)
            double log_probs[400];
            int map_cat[400];
            int map_dist[400];
            int n_opts = 0;

            for(int i_ac=0; i_ac < n_active; i_ac++) {
                int cat_idx = active_cats[i_ac];
                
                // d=0 (Null)
                log_probs[n_opts] = log(model->p[0 * ncat + cat_idx]);
                map_cat[n_opts] = cat_idx;
                map_dist[n_opts] = 0;
                n_opts++;

                // d > 0
                for (int i = 1; i < ndist; i++) {
                    double logdetV = log(model->gp[i] * zz / model->vare + 1.0);
                    double uhat = rhs / (zz + vare_gp[i]);
                    double log_lik = -0.5 * (logdetV - (rhs * uhat / model->vare));
                    
                    log_probs[n_opts] = log_lik + log(model->p[i * ncat + cat_idx]);
                    map_cat[n_opts] = cat_idx;
                    map_dist[n_opts] = i;
                    n_opts++;
                }
            }

            // Log-sum-exp
            double max_log = log_probs[0];
            for(int i=1; i<n_opts; i++) if(log_probs[i] > max_log) max_log = log_probs[i];
            
            double sum_prob = 0.0;
            for(int i=0; i<n_opts; i++) {
                log_probs[i] = exp(log_probs[i] - max_log);
                sum_prob += log_probs[i];
            }

            // Sample
            double r = rand_uniform(0.0, 1.0) * sum_prob;
            double cum = 0.0;
            int selected_idx = n_opts - 1; // Default to last if precision issues
            for(int i=0; i<n_opts; i++) {
                cum += log_probs[i];
                if(r < cum) { selected_idx = i; break; }
            }

            int choice = map_dist[selected_idx];
            int chosen_cat = map_cat[selected_idx];

            double gk_new = 0.0;
            if (choice == 0) {
                gk_new = 0.0;
            } else {
                double v1 = zz + vare_gp[choice];
                gk_new = rand_normal(rhs / v1, sqrt(model->vare / v1));
                included++;
                model->gh[choice * ncat + chosen_cat] += gk_new * gk_new / config->gpin[choice];
            }

            // Update yadj only if gk changed
            double diff = gk - gk_new;
            if (diff != 0.0) {
                 for (int i = 0; i < nt; i++) model->yadj[i] += z[i] * diff;
            }

            model->g[k] = gk_new;
            // Record the choice for the chosen active category
            model->snptracker[k * ncat + chosen_cat] = choice;
            model->mc[choice * ncat + chosen_cat]++;
        }

        // 4. Update hyperparameters
        update_vara(config, model);
        update_pi(config, model);

        if (fp_hyp) fprintf(fp_hyp, "%d %f %f %f\n", rep, model->mu, model->vara, model->vare);

        if (rep % 10 == 0 || rep == config->numit) {
            printf("Iteration %d: included %d, vare %f\n", rep, included, model->vare);
        }
    }

    if (fp_hyp) fclose(fp_hyp);
    free(permvec);
    free(s);
    free(stemp);
    free(vare_gp);
}

void free_bayes_model(BayesModel *model) {
    if (!model) return;
    free(model->g);
    free(model->p);
    free(model->gp);
    free(model->yadj);
    free(model->xpx);
    free(model->snptracker);
    free(model->mc);
    free(model->gh);
}
