#!/usr/bin/env Rscript
source("src_R_mimic/bayesR_mimic.R")

# Help/Arg parser
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

# --- Parse Arguments ---
in_prefix  <- get_arg("-bfile", "input")
out_prefix <- get_arg("-out", "output")
num_it     <- as.integer(get_arg("-numit", 100))
burn_in    <- as.integer(get_arg("-burnin", 0))
seed       <- as.integer(get_arg("-seed", 42))
cat_file   <- get_arg("-catfile", "")

set.seed(seed)

# --- Load Data ---
cat("Loading PLINK data using genio...\n")
if (!requireNamespace("genio", quietly = TRUE)) {
  stop("Package 'genio' is required. Run: install.packages('genio')")
}

data <- genio::read_plink(in_prefix)
genotypes  <- data$X # SNPs are rows, inds are cols

# genio stores phenotypes in data$fam$pheno
phenotypes <- data$fam$pheno

# IMPORTANT: C and Fortran do not handle PLINK's -9 missing code and treat it as a valid trait.
# To match their frequencies, we must also include all individuals.
# genio might have converted -9 to NA, so we put it back if needed for consistency.
phenotypes[is.na(phenotypes)] <- -9 

# Target all individuals to match C/Fortran nt=nind behavior
train_mask <- rep(TRUE, length(phenotypes))

y_train <- phenotypes[train_mask]
X_train <- genotypes[, train_mask]
nloci   <- nrow(genotypes)
nind    <- sum(train_mask)

cat(sprintf("Loaded %d individuals (%d training) and %d SNPs.\n", length(phenotypes), nind, nloci))

# --- Preprocessing (Center & Scale) ---
# Project logic: q = sum(X)/(2*n_no_miss)
# frequency = q
# mean = 2*q
# sd = sqrt(2*q*(1-q))
# scaled_X = (X - mean) / sd

cat("Preprocessing and centering...\n")
frequencies <- numeric(nloci)
X_processed <- matrix(0, nrow = length(y_train), ncol = nloci)

for (j in 1:nloci) {
  # Flip from genio's default (A1 dosage) to project default (A2 dosage)
  snp_vals <- X_train[j, ] 
  
  no_miss_mask <- !is.na(snp_vals)
  
  if (sum(no_miss_mask) > 0) {
    q <- sum(snp_vals[no_miss_mask]) / (2 * sum(no_miss_mask))
  } else {
    q <- 0.5
  }
  frequencies[j] <- q
  
  if (q <= 0 || q >= 1) {
    X_processed[, j] <- 0
  } else {
    mean_val <- 2 * q
    sd_val   <- sqrt(2 * q * (1 - q))
    # Replace NAs with mean
    snp_vals[is.na(snp_vals)] <- mean_val
    X_processed[, j] <- (snp_vals - mean_val) / sd_val
  }
}

# Save frequencies to match verification
freq_file <- paste0(out_prefix, ".frq")
write.table(data.frame(sprintf("%10.6f", frequencies)), 
            file = freq_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("Frequency file saved to:", freq_file, "\n")

# --- Load Categories ---
# Default: List of integer 1s
snp_categories <- replicate(nloci, as.integer(c(1)), simplify = FALSE)
if (cat_file != "" && file.exists(cat_file)) {
    cat_data <- as.matrix(read.table(cat_file, header = FALSE))
    if (nrow(cat_data) == nloci) {
        # Match C logic: find all columns with a 1 (1-based index)
        # Store as a list of vectors for BayesRCpi
        snp_categories <- apply(cat_data, 1, function(row) {
            ones <- which(row == 1)
            if (length(ones) > 0) return(as.integer(ones)) else return(as.integer(c(1)))
        })
        # If result is a matrix (because all rows same length), force to list
        if (is.matrix(snp_categories)) {
          snp_categories <- split(snp_categories, col(snp_categories))
        } else if (!is.list(snp_categories)) {
           # If vector (all single annot)
           snp_categories <- as.list(snp_categories)
        }
        cat("Loaded categories from matrix file (List format).\n")
    } else {
        warning("Category file row count does not match SNP count. Using default (1).")
        snp_categories <- replicate(nloci, as.integer(c(1)), simplify = FALSE)
    }
} else {
    # Default list
    snp_categories <- replicate(nloci, as.integer(c(1)), simplify = FALSE)
}

# --- Run Mimic Kernel ---
config <- list(
  num_iterations = num_it,
  num_distributions = 4,
  mixture_variance_scales = c(0, 0.0001, 0.001, 0.01),
  initial_genetic_variance = 0.01,
  initial_residual_variance = 0.01,
  genetic_variance_df = -2.0,
  dirichlet_prior_counts = rep(1.0, 4)
)

results <- run_bayesR_mimic_mcmc(y_train, X_processed, snp_categories, config)

# Save hyperparameter chains to match verify.sh expected format
hyp_file <- paste0(out_prefix, ".hyp")
hyp_data <- data.frame(
  Iteration = 1:num_it,
  mu = results$intercept_chain,
  vara = results$genetic_variance_chain,
  vare = results$residual_variance_chain
)
write.table(hyp_data, file = hyp_file, quote = FALSE, row.names = FALSE, sep = " ")
cat("Hyperparameter chains saved to:", hyp_file, "\n")
cat("R Mimic version run complete.\n")
