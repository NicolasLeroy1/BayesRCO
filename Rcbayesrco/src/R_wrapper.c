#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "../../new_src_c/libbayesrco_stats/bayesrco_stats.h"
#include <stdint.h>

/* Helper to convert R numeric matrix to C MCMCParams and run model */
static SEXP run_model_generic(int model_type,
                            SEXP genotypes, SEXP phenotypes, SEXP categories,
                            SEXP num_iterations, SEXP burnin_iterations, SEXP thinning_interval,
                            SEXP num_distributions, SEXP variance_scaling_factors, SEXP dirichlet_priors,
                            SEXP vara_ap, SEXP vare_ap, SEXP dfvara, SEXP dfvare, SEXP seed) {
    
    int nloci = (int)ncols(genotypes);
    int nt = (int)nrows(genotypes);
    int ncat = (int)ncols(categories);
    int nind = nt; /* For now, nind = nt in R wrapper unless we add prediction */

    MCMCParams params;
    params.num_iterations = asInteger(num_iterations);
    params.burnin_iterations = asInteger(burnin_iterations);
    params.thinning_interval = asInteger(thinning_interval);
    params.num_distributions = asInteger(num_distributions);
    params.num_categories = ncat;
    params.permute = true;
    params.VCE = true;
    params.dfvara = asReal(dfvara);
    params.dfvare = asReal(dfvare);
    params.vara_ap = asReal(vara_ap);
    params.vare_ap = asReal(vare_ap);
    params.variance_genetic = params.vara_ap;
    params.variance_residual = params.vare_ap;
    
    params.variance_scaling_factors = REAL(variance_scaling_factors);
    params.dirichlet_priors = REAL(dirichlet_priors);

    /* trains array: all 0 (training) for now */
    int *c_trains = (int*)calloc(nt, sizeof(int));
    if (!c_trains) error("Failed to allocate trains array");

    MCMCResults results;
    if (allocate_results(&results, nloci, params.num_distributions, ncat, nt) != SUCCESS) {
        free(c_trains);
        error("Failed to allocate MCMC results");
    }

    int ret;
    uint64_t c_seed = (uint64_t)asReal(seed);

    if (model_type == 0) {
        ret = run_bayesrco_mixture(&params, REAL(genotypes), REAL(phenotypes), INTEGER(categories), 
                                  nloci, nt, ncat, nt, c_trains, c_seed, &results);
    } else if (model_type == 1) {
        ret = run_bayesrco_additive(&params, REAL(genotypes), REAL(phenotypes), INTEGER(categories), 
                                   nloci, nt, ncat, nt, c_trains, c_seed, &results);
    } else {
        ret = run_bayesrco_bayesCpi(&params, REAL(genotypes), REAL(phenotypes), INTEGER(categories), 
                                   nloci, nt, ncat, nt, c_trains, c_seed, &results);
    }

    free(c_trains);

    if (ret != SUCCESS) {
        free_results(&results);
        error("MCMC execution failed with code %d", ret);
    }

    /* Create R list to return results */
    SEXP res_list = PROTECT(allocVector(VECSXP, 5));
    SEXP names = PROTECT(allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, mkChar("snp_effects"));
    SET_STRING_ELT(names, 1, mkChar("variance_genetic"));
    SET_STRING_ELT(names, 2, mkChar("variance_residual"));
    SET_STRING_ELT(names, 3, mkChar("pi"));
    SET_STRING_ELT(names, 4, mkChar("predicted_values"));
    setAttrib(res_list, R_NamesSymbol, names);

    /* SNP effects */
    SEXP r_snp_effects = allocVector(REALSXP, nloci);
    for (int i = 0; i < nloci; i++) REAL(r_snp_effects)[i] = results.posterior_means[i];
    SET_VECTOR_ELT(res_list, 0, r_snp_effects);

    /* Hyperparameters */
    SET_VECTOR_ELT(res_list, 1, ScalarReal(results.vara));
    SET_VECTOR_ELT(res_list, 2, ScalarReal(results.vare));

    /* Pi (mixture proportions) */
    SEXP r_pi = allocMatrix(REALSXP, params.num_distributions, ncat);
    for (int j = 0; j < ncat; j++) {
        for (int i = 0; i < params.num_distributions; i++) {
            /* 
               Results are stored as ndist * ncat flattened in ROW-major order internally? 
               Wait, bayesrco_types.h IDX2(i, j, nj) is (i*nj + j).
               In state and storage:
               p: size ndist * ncat. indexed with IDX2(i, j, ncat) -> i*ncat + j.
               So it is effectively a row-major ndist x ncat matrix.
               R matrices are column-major. To get an ndist x ncat matrix in R,
               we fill the R REAL(r_pi) pointer.
            */
            REAL(r_pi)[j * params.num_distributions + i] = results.distribution_probs[i * ncat + j];
        }
    }
    SET_VECTOR_ELT(res_list, 3, r_pi);
    
    /* Predicted values */
    SEXP r_pred = allocVector(REALSXP, nt);
    for (int i = 0; i < nt; i++) REAL(r_pred)[i] = results.predicted_values[i];
    SET_VECTOR_ELT(res_list, 4, r_pred);

    free_results(&results);
    UNPROTECT(2);
    return res_list;
}

SEXP C_run_bayesrco_mixture(SEXP genotypes, SEXP phenotypes, SEXP categories,
                          SEXP num_iterations, SEXP burnin_iterations, SEXP thinning_interval,
                          SEXP num_distributions, SEXP variance_scaling_factors, SEXP dirichlet_priors,
                          SEXP vara_ap, SEXP vare_ap, SEXP dfvara, SEXP dfvare, SEXP seed) {
    return run_model_generic(0, genotypes, phenotypes, categories, num_iterations, burnin_iterations, thinning_interval,
                           num_distributions, variance_scaling_factors, dirichlet_priors, vara_ap, vare_ap, dfvara, dfvare, seed);
}

SEXP C_run_bayesrco_additive(SEXP genotypes, SEXP phenotypes, SEXP categories,
                           SEXP num_iterations, SEXP burnin_iterations, SEXP thinning_interval,
                           SEXP num_distributions, SEXP variance_scaling_factors, SEXP dirichlet_priors,
                           SEXP vara_ap, SEXP vare_ap, SEXP dfvara, SEXP dfvare, SEXP seed) {
    return run_model_generic(1, genotypes, phenotypes, categories, num_iterations, burnin_iterations, thinning_interval,
                           num_distributions, variance_scaling_factors, dirichlet_priors, vara_ap, vare_ap, dfvara, dfvare, seed);
}

SEXP C_run_bayesrco_bayesCpi(SEXP genotypes, SEXP phenotypes, SEXP categories,
                           SEXP num_iterations, SEXP burnin_iterations, SEXP thinning_interval,
                           SEXP num_distributions, SEXP variance_scaling_factors, SEXP dirichlet_priors,
                           SEXP vara_ap, SEXP vare_ap, SEXP dfvara, SEXP dfvare, SEXP seed) {
    return run_model_generic(2, genotypes, phenotypes, categories, num_iterations, burnin_iterations, thinning_interval,
                           num_distributions, variance_scaling_factors, dirichlet_priors, vara_ap, vare_ap, dfvara, dfvare, seed);
}

static const R_CallMethodDef CallEntries[] = {
    {"C_run_bayesrco_mixture", (DL_FUNC) &C_run_bayesrco_mixture, 14},
    {"C_run_bayesrco_additive", (DL_FUNC) &C_run_bayesrco_additive, 14},
    {"C_run_bayesrco_bayesCpi", (DL_FUNC) &C_run_bayesrco_bayesCpi, 14},
    {NULL, NULL, 0}
};

void R_init_Rcbayesrco(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
