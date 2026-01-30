#' BayesR Mimic Kernel (Mimics Old Fortran)
#'
#' This file implements a Gibbs sampler designed to replicate the specific
#' statistical behavior of the 'Old Fortran' implementation of BayesRCO.
#'
#' KEY DIFFERENCES from standard BayesRC:
#' 1. Genetic Variance (Va) update uses SIMPLE Sum of Squares (sum(g^2))
#'    instead of Weighted Sum of Squares (sum(g^2 / gamma)).
#' 2. Residual Variance (Ve) update uses df = N - 2 instead of N + 3.
#'

# --- Statistical Helper Functions ---

sample_inverse_gamma <- function(n, shape, scale) {
  1.0 / rgamma(n, shape, rate = scale)
}

sample_dirichlet <- function(alpha_counts) {
  gamma_samples <- rgamma(length(alpha_counts), alpha_counts, 1.0)
  return(gamma_samples / sum(gamma_samples))
}

# --- Parameter Update Functions (MIMIC MODE) ---

#' Update Residual Variance (Mimic)
#' Law: sigma_e^2 | ... ~ Scaled-Inv-Chi2( n_ind + df_fixed, RSS )
#' Fortran dfvare defaults to -2.
update_residual_variance_mimic <- function(residuals, num_individuals) {
  sum_squared_residuals <- sum(residuals^2)
  # Mimic Fortran: df = n_ind + dfvare (where dfvare = -2)
  df_eff <- num_individuals - 2.0
  return(sum_squared_residuals / rchisq(1, df = df_eff))
}

#' Update Intercept (Standard)
update_intercept <- function(residuals, intercept, residual_variance, num_individuals) {
  y_adj <- residuals + intercept
  intercept_mean <- mean(y_adj)
  intercept_sd   <- sqrt(residual_variance / num_individuals)
  new_intercept  <- rnorm(1, mean = intercept_mean, sd = intercept_sd)
  return(list(value = new_intercept, residuals = y_adj - new_intercept))
}

#' Update SNP Effects (Standard BayesRCpi logic)
#' Note: We keep the "Overlapping" sampling logic from R because it's superior/cleaner,
#' provided that with the Mimic Variance updates, the results align.
update_snp_effects <- function(residuals, snp_effects, genotype_matrix, snp_category_list, 
                               snp_diagonal_xpx, mixture_proportions, mixture_variances, 
                               residual_variance, noise_to_signal_ratios, config) {
  
  num_snps        <- ncol(genotype_matrix)
  num_dists       <- length(mixture_variances)
  num_categories  <- ncol(mixture_proportions)
  
  snp_counts_per_dist_cat <- matrix(0, nrow = num_dists, ncol = num_categories)
  # For Mimic: We track SIMPLE Sum of Squares (g^2) 
  sum_sq_simple_effects  <- matrix(0, nrow = num_dists, ncol = num_categories)
  
  shuffled_snps <- sample(1:num_snps)
  
  for (k in shuffled_snps) {
    diagonal_zz <- snp_diagonal_xpx[k]
    current_gk  <- snp_effects[k]
    possible_annotations <- snp_category_list[[k]]
    snp_column  <- genotype_matrix[, k]
    
    if (current_gk != 0.0) {
      residuals <- residuals + snp_column * current_gk
    }
    
    # Optimized dot product using BLAS
    rhs_dot_product <- drop(crossprod(snp_column, residuals))
    
    # Vectorized calculation for dists 2:num_dists
    # log_likelihoods[1] = 0 (Null model)
    # Avoid loop by calculating for all k > 1 at once
    
    d_indices <- 2:num_dists
    
    # Calculate conditional parameters for all active distributions
    denom <- diagonal_zz + noise_to_signal_ratios[d_indices]
    cond_means_active <- rhs_dot_product / denom
    cond_vars_active  <- residual_variance / denom
    log_det_active    <- log(mixture_variances[d_indices] * diagonal_zz / residual_variance + 1.0)
    
    # L(d) calculation
    log_lik_active <- -0.5 * (log_det_active - (rhs_dot_product * cond_means_active / residual_variance))
    
    # Combine with Null model
    log_likelihoods <- c(0, log_lik_active)
    cond_means      <- c(0, cond_means_active)
    cond_vars       <- c(0, cond_vars_active)
    
    # Flatten probability calculation for (d, a) pairs
    # log_probs = L(d) + log(pi_{d,a})
    # We can compute this using outer product-like logic or simpler vectorized addition
    
    # Pre-allocate if num_opts is large, but usually small (dists * n_annot)
    # Using matrix formulation for selecting
    
    probs_matrix <- matrix(nrow = num_dists, ncol = length(possible_annotations))
    
    for (i in seq_along(possible_annotations)) {
      a <- possible_annotations[i]
      probs_matrix[, i] <- log_likelihoods + log(mixture_proportions[, a])
    }
    
    # Flatten for sampling
    log_probs_flat <- as.vector(probs_matrix)
    
    max_log_prob <- max(log_probs_flat)
    stable_probs <- exp(log_probs_flat - max_log_prob)
    sampled_flat <- sample.int(length(log_probs_flat), 1, prob = stable_probs)
    
    # Decode index
    # Matrix was (num_dists x num_annots), traversed column-major? 
    # as.vector flattens COLUMNS first (d varies fast, a slow)
    # So index = (a_idx - 1)*num_dists + d_idx
    
    selected_dist_idx <- (sampled_flat - 1) %% num_dists + 1
    selected_annot_idx <- (sampled_flat - 1) %/% num_dists + 1
    selected_annot <- possible_annotations[selected_annot_idx]
    
    if (selected_dist_idx == 1) {
      new_gk <- 0.0
    } else {
      new_gk <- rnorm(1, mean = cond_means[selected_dist_idx], sd = sqrt(cond_vars[selected_dist_idx]))
      residuals <- residuals - snp_column * new_gk
      
      # MIMIC: Accumulate SIMPLE g^2, NOT weighted g^2/gamma
      sum_sq_simple_effects[selected_dist_idx, selected_annot] <- sum_sq_simple_effects[selected_dist_idx, selected_annot] + new_gk^2
    }
    
    snp_effects[k] <- new_gk
    snp_counts_per_dist_cat[selected_dist_idx, selected_annot] <- snp_counts_per_dist_cat[selected_dist_idx, selected_annot] + 1
  }
  
  return(list(
    snp_effects = snp_effects,
    residuals = residuals,
    counts = snp_counts_per_dist_cat,
    simple_ss = sum_sq_simple_effects
  ))
}

#' Update Genetic Variance (Mimic)
#' Law: sigma_a^2 | ... ~ Scaled-Inv-Chi2
#' Fortran logic:
#' scale = ( included * sum(g**2 / 1.0??) ) / (included + df)
#' Wait, Fortran code:
#' scale=(dble(included)*sum(g**2) + vara_ap*dfvara)/(dfvara+dble(included))
#' vara=rand_scaled_inverse_chi_square(dble(included)+dfvara,scale)
#' 
#' This implies it treats ALL effects as if they come from the SAME distribution with variance 'vara',
#' effectively ignoring the gamma scaling in the posterior update of vara itself.
update_genetic_variance_mimic <- function(snp_counts_per_dist_cat, sum_sq_simple_effects, config) {
  num_dists <- nrow(snp_counts_per_dist_cat)
  
  total_snps_included <- sum(snp_counts_per_dist_cat[2:num_dists, ]) # Sum of non-null counts
  total_simple_ss     <- sum(sum_sq_simple_effects[2:num_dists, ])   # Sum of g^2
  
  df_prior <- -2.0
  prior_scale <- 0.0 # because df=-2 implies uninformative
  
  # Effective DF and Scale
  # If df=-2, let's assume it behaves like N-2.
  df_post <- total_snps_included + df_prior
  
  numerator <- total_snps_included * total_simple_ss # The 'bug' replication

  
  if (df_post > 0) {
    return(numerator / rchisq(1, df_post))
  } else {
    return(config$initial_genetic_variance) # Fallback if no SNPs included
  }
}

update_mixture_proportions <- function(snp_counts_per_dist_cat, dirichlet_prior_counts) {
  num_categories <- ncol(snp_counts_per_dist_cat)
  num_dists      <- nrow(snp_counts_per_dist_cat)
  new_pi <- matrix(0, nrow = num_dists, ncol = num_categories)
  for (cat_idx in 1:num_categories) {
    posterior_alpha <- dirichlet_prior_counts + snp_counts_per_dist_cat[, cat_idx]
    new_pi[, cat_idx] <- sample_dirichlet(posterior_alpha)
  }
  return(new_pi)
}

# --- Main Driver Function ---

run_bayesR_mimic_mcmc <- function(phenotypes, genotype_matrix, snp_category_list, config) {
  num_individuals <- length(phenotypes)
  num_snps        <- ncol(genotype_matrix)
  num_dists       <- config$num_distributions
  
  all_cats <- unlist(snp_category_list)
  num_categories <- if (length(all_cats) == 0) 1 else max(all_cats)
  
  intercept        <- mean(phenotypes)
  snp_effects      <- rep(0.0, num_snps)
  genetic_variance <- config$initial_genetic_variance
  residual_variance <- config$initial_residual_variance
  # Initialize mixture proportions (pi) to match Old Fortran
  # Fortran: p(1)=0.5, others inversely proportional to gpin
  mixture_proportions <- matrix(0, nrow = num_dists, ncol = num_categories)
  # Assuming Row 1 is Null (scale 0)
  mixture_proportions[1, ] <- 0.5
  
  if (num_dists > 1) {
    # Scales for 2:num_dists
    scales <- config$mixture_variance_scales[2:num_dists]
    # Check for zeros to avoid Inf
    scales[scales == 0] <- 1e-9 
    
    inv_scales <- 1.0 / scales
    # Normalize to sum to 0.5
    norm_factors <- 0.5 * inv_scales / sum(inv_scales)
    
    for (i in 1:(num_dists - 1)) {
      mixture_proportions[i + 1, ] <- norm_factors[i]
    }
  }
  
  # Initialize SNP effects (g) to match Old Fortran
  # g = sqrt(vara / (0.5 * nloci))
  initial_g_sd <- sqrt(genetic_variance / (0.5 * num_snps))
  snp_effects  <- rep(initial_g_sd, num_snps)
  
  # Calculate initial residuals consistent with initialized effects
  # residuals = y - mu - Xg
  residuals <- phenotypes - intercept - drop(genotype_matrix %*% snp_effects)
  
  snp_diagonal_xpx <- colSums(genotype_matrix^2)
  
  chain_intercept        <- numeric(config$num_iterations)
  chain_genetic_variance <- numeric(config$num_iterations)
  chain_residual_variance <- numeric(config$num_iterations)
  
  for (iteration in 1:config$num_iterations) {
    intercept_update <- update_intercept(residuals, intercept, residual_variance, num_individuals)
    intercept <- intercept_update$value
    residuals <- intercept_update$residuals
    
    mixture_variances      <- config$mixture_variance_scales * genetic_variance
    noise_to_signal_ratios <- residual_variance / mixture_variances
    
    snp_update <- update_snp_effects(
      residuals, snp_effects, genotype_matrix, snp_category_list,
      snp_diagonal_xpx, mixture_proportions, mixture_variances,
      residual_variance, noise_to_signal_ratios, config
    )
    
    snp_effects <- snp_update$snp_effects
    residuals   <- snp_update$residuals
    
    genetic_variance <- update_genetic_variance_mimic(snp_update$counts, snp_update$simple_ss, config)
    mixture_proportions <- update_mixture_proportions(snp_update$counts, config$dirichlet_prior_counts)
    
    # Move residual variance update to END of loop to match Old Fortran behavior
    # This leads to a "hot start" (low shrinkage) in the first iteration if init var is low.
    residual_variance <- update_residual_variance_mimic(residuals, num_individuals)
    
    chain_intercept[iteration]        <- intercept
    chain_genetic_variance[iteration] <- genetic_variance
    chain_residual_variance[iteration] <- residual_variance
    
    if (iteration %% 10 == 0) {
      cat(sprintf("Mimic Iter %d: Va=%.6f, Ve=%.6f, Poly=%d\n", 
                  iteration, genetic_variance, residual_variance, sum(snp_effects != 0.0)))
    }
  }
  
  return(list(
    intercept_chain = chain_intercept,
    genetic_variance_chain = chain_genetic_variance,
    residual_variance_chain = chain_residual_variance,
    final_snp_effects = snp_effects,
    final_mixture_proportions = mixture_proportions
  ))
}
