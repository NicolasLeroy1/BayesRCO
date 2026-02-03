library(Rcbayesrco)

# Create dummy data
nind <- 10
nloci <- 5
ncat <- 2

set.seed(123)
genotypes <- matrix(rnorm(nind * nloci), nrow = nind, ncol = nloci)
phenotypes <- rnorm(nind)
categories <- matrix(sample(0:1, nloci * ncat, replace = TRUE), nrow = nloci, ncol = ncat)

# Run mixture model
print("Testing mixture model...")
results <- run_bayesrco_mixture(
  genotypes = genotypes,
  phenotypes = phenotypes,
  categories = categories,
  num_iterations = 100,
  burnin_iterations = 20,
  thinning_interval = 1,
  num_distributions = 4,
  seed = 42
)

print("Mixture Model Results:")
print(names(results))
print(paste("Vara:", results$variance_genetic))
print(paste("Vare:", results$variance_residual))
print("Pi matrix dimension:")
print(dim(results$pi))

# Run additive model
print("Testing additive model...")
results_add <- run_bayesrco_additive(
  genotypes = genotypes,
  phenotypes = phenotypes,
  categories = categories,
  num_iterations = 100,
  burnin_iterations = 20,
  thinning_interval = 1,
  num_distributions = 4,
  seed = 42
)
print("Additive Model Success!")

# Run BayesCpi
print("Testing BayesCpi...")
results_cpi <- run_bayesrco_bayesCpi(
  genotypes = genotypes,
  phenotypes = phenotypes,
  categories = categories,
  num_iterations = 100,
  burnin_iterations = 20,
  thinning_interval = 1,
  num_distributions = 4,
  seed = 42
)
print("BayesCpi Success!")
