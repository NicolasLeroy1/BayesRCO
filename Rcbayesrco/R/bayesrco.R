#' Run BayesRCO Mixture Model
#' @export
run_bayesrco_mixture <- function(genotypes, phenotypes, categories,
                                 num_iterations = 50000,
                                 burnin_iterations = 20000,
                                 thinning_interval = 10,
                                 num_distributions = 4,
                                 variance_scaling_factors = c(0.0, 0.0001, 0.001, 0.01),
                                 dirichlet_priors = rep(1.0, num_distributions),
                                 vara_ap = 0.01,
                                 vare_ap = 0.01,
                                 dfvara = -2.0,
                                 dfvare = -2.0,
                                 seed = 0) {
  
  # Validation
  if (!is.matrix(genotypes)) stop("genotypes must be a matrix")
  if (length(phenotypes) != nrow(genotypes)) stop("phenotypes length must match genotypes rows")
  if (!is.matrix(categories)) stop("categories must be a matrix")
  
  # Call C wrapper
  .Call("C_run_bayesrco_mixture",
        genotypes, phenotypes, categories,
        as.integer(num_iterations),
        as.integer(burnin_iterations),
        as.integer(thinning_interval),
        as.integer(num_distributions),
        as.double(variance_scaling_factors),
        as.double(dirichlet_priors),
        as.double(vara_ap),
        as.double(vare_ap),
        as.double(dfvara),
        as.double(dfvare),
        as.double(seed),
        PACKAGE = "Rcbayesrco")
}

#' Run BayesRCO Additive Model
#' @export
run_bayesrco_additive <- function(genotypes, phenotypes, categories,
                                  num_iterations = 50000,
                                  burnin_iterations = 20000,
                                  thinning_interval = 10,
                                  num_distributions = 4,
                                  variance_scaling_factors = c(0.0, 0.0001, 0.001, 0.01),
                                  dirichlet_priors = rep(1.0, num_distributions),
                                  vara_ap = 0.01,
                                  vare_ap = 0.01,
                                  dfvara = -2.0,
                                  dfvare = -2.0,
                                  seed = 0) {
  
  .Call("C_run_bayesrco_additive",
        genotypes, phenotypes, categories,
        as.integer(num_iterations),
        as.integer(burnin_iterations),
        as.integer(thinning_interval),
        as.integer(num_distributions),
        as.double(variance_scaling_factors),
        as.double(dirichlet_priors),
        as.double(vara_ap),
        as.double(vare_ap),
        as.double(dfvara),
        as.double(dfvare),
        as.double(seed),
        PACKAGE = "Rcbayesrco")
}

#' Run BayesCpi Model
#' @export
run_bayesrco_bayesCpi <- function(genotypes, phenotypes, categories,
                                  num_iterations = 50000,
                                  burnin_iterations = 20000,
                                  thinning_interval = 10,
                                  num_distributions = 4,
                                  variance_scaling_factors = c(0.0, 0.0001, 0.001, 0.01),
                                  dirichlet_priors = rep(1.0, num_distributions),
                                  vara_ap = 0.01,
                                  vare_ap = 0.01,
                                  dfvara = -2.0,
                                  dfvare = -2.0,
                                  seed = 0) {
  
  .Call("C_run_bayesrco_bayesCpi",
        genotypes, phenotypes, categories,
        as.integer(num_iterations),
        as.integer(burnin_iterations),
        as.integer(thinning_interval),
        as.integer(num_distributions),
        as.double(variance_scaling_factors),
        as.double(dirichlet_priors),
        as.double(vara_ap),
        as.double(vare_ap),
        as.double(dfvara),
        as.double(dfvare),
        as.double(seed),
        PACKAGE = "Rcbayesrco")
}
