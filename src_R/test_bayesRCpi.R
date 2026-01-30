source("bayesRCpi.R")

# Generate Synthetic Data
set.seed(42)
n <- 200
p <- 100
X <- matrix(rnorm(n * p), n, p)
C <- sample(1:2, p, replace = TRUE) # 2 categories

# True effects
g_true <- rep(0, p)
g_true[sample(1:p, 10)] <- rnorm(10, 0, 0.5)
y <- X %*% g_true + rnorm(n, 0, 0.2)

# Config
config <- list(
  num_iterations = 100,
  num_distributions = 4,
  mixture_variance_scales = c(0, 0.0001, 0.001, 0.01),
  initial_genetic_variance = 0.01,
  initial_residual_variance = 0.01,
  genetic_variance_df = -2.0,
  dirichlet_prior_counts = rep(1.0, 4)
)

# Run Kernel
cat("Starting BayesRCpi Kernel Verification (Readable Version)...\n")
results <- run_bayesRCpi_mcmc(y, X, C, config)

cat("\nVerification Complete.\n")
cat("Final residual variance sample:", tail(results$residual_variance_chain, 1), "\n")
cat("Final genetic variance sample:", tail(results$genetic_variance_chain, 1), "\n")
cat("Number of non-zero SNPs:", sum(results$final_snp_effects != 0), "out of", p, "\n")
