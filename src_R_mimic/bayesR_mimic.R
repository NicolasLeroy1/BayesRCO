#' BayesR Mimic Kernel (Mimics Old Fortran)
#'
#' This file implements a Gibbs sampler designed to replicate the specific
#' statistical behavior of the 'Old Fortran' implementation of BayesRCO.

# --- Statistical Helper Functions ---

# rand_scaled_inverse_chi_square(df, scale) in Fortran:
# res = (df * scale) / rchisq(df)
# Note: For Va, Fortran uses scale = (included * sum(g^2) + vara_ap*dfvara)/(dfvara+included)
# So vara = (included * sum(g^2) + ...) / rchisq(included + dfvara)

sample_scaled_inv_chi_sq <- function(df, scale) {
  if (df <= 0) return(0)
  return((df * scale) / rchisq(1, df))
}

sample_dirichlet <- function(alpha_counts) {
  gamma_samples <- rgamma(length(alpha_counts), alpha_counts, 1.0)
  if (sum(gamma_samples) == 0) return(rep(1/length(alpha_counts), length(alpha_counts)))
  return(gamma_samples / sum(gamma_samples))
}

# --- Parameter Update Functions (MIMIC MODE) ---

#' Update Residual Variance (Mimic)
#' Fortran VCE=true logic:
#' vare=(dot_product(yadj,yadj)+vare_ap*dfvare)/ (nnind+dfvare)
#' vare=rand_scaled_inverse_chi_square(nnind+dfvare,vare)
update_residual_variance_mimic <- function(residuals, num_individuals, df_prior, vare_ap) {
  df_post <- num_individuals + df_prior
  ss_post <- sum(residuals^2) + vare_ap * df_prior
  scale_post <- ss_post / df_post
  return(sample_scaled_inv_chi_sq(df_post, scale_post))
}

#' Update Intercept (Standard)
#' Fortran: mu=rand_normal(sum(yadj)/nnind, dsqrt(vare/nnind))
update_intercept <- function(residuals, intercept, residual_variance, num_individuals) {
  y_adj <- residuals + intercept
  intercept_mean <- sum(y_adj) / num_individuals
  intercept_sd   <- sqrt(residual_variance / num_individuals)
  new_intercept  <- rnorm(1, mean = intercept_mean, sd = intercept_sd)
  return(list(value = new_intercept, residuals = y_adj - new_intercept))
}

#' Update SNP Effects (Standard BayesRCpi logic)
update_snp_effects <- function(residuals, snp_effects, genotype_matrix, snp_category_list, 
                               snp_diagonal_xpx, mixture_proportions, mixture_variances, 
                               residual_variance, noise_to_signal_ratios, config) {
  
  num_snps        <- ncol(genotype_matrix)
  num_dists       <- length(mixture_variances)
  num_categories  <- ncol(mixture_proportions)
  
  snp_counts_per_dist_cat <- matrix(0, nrow = num_dists, ncol = num_categories)
  sum_sq_simple_effects  <- matrix(0, nrow = num_dists, ncol = num_categories)
  
  # Fortran uses a custom permutation logic if permute=T
  shuffled_snps <- sample(1:num_snps)
  
  for (k in shuffled_snps) {
    diagonal_zz <- snp_diagonal_xpx[k]
    current_gk  <- snp_effects[k]
    # For mimic, we assume only 1 annotation per SNP as in simple BayesR
    # but we handle the list format for compatibility
    possible_annotations <- snp_category_list[[k]]
    snp_col <- genotype_matrix[, k]
    
    if (current_gk != 0.0) {
      residuals <- residuals + snp_col * current_gk
    }
    
    rhs <- sum(snp_col * residuals) # dot_product
    
    d_indices <- 2:num_dists
    
    # Calculate conditional parameters
    denom <- diagonal_zz + noise_to_signal_ratios[d_indices]
    cond_means_active <- rhs / denom
    cond_vars_active  <- residual_variance / denom
    
    # Fortran logdetV = dlog(gp*zz_vare + 1.0)
    # R: log_det = log(mixture_variances * diagonal_zz / residual_variance + 1.0)
    log_det_active <- log(mixture_variances[d_indices] * diagonal_zz / residual_variance + 1.0)
    
    # Fortran: s(kk) = -0.5*(logdetV - (rhs*uhat/vare)) + log_p
    log_lik_active <- -0.5 * (log_det_active - (rhs * cond_means_active / residual_variance))
    
    log_likelihoods <- c(0, log_lik_active)
    
    # FLAT SAMPLING (d, a)
    # We only take the first annotation for mimic if multiple provided, 
    # but let's stick to the logic for robustness.
    probs_matrix <- matrix(nrow = num_dists, ncol = length(possible_annotations))
    for (i in seq_along(possible_annotations)) {
      a <- possible_annotations[i]
      probs_matrix[, i] <- log_likelihoods + log(mixture_proportions[, a])
    }
    
    log_probs_flat <- as.vector(probs_matrix)
    max_log_prob <- max(log_probs_flat)
    stable_probs <- exp(log_probs_flat - max_log_prob)
    sampled_flat <- sample.int(length(log_probs_flat), 1, prob = stable_probs)
    
    selected_dist_idx <- (sampled_flat - 1) %% num_dists + 1
    selected_annot_idx <- (sampled_flat - 1) %/% num_dists + 1
    selected_annot <- possible_annotations[selected_annot_idx]
    
    if (selected_dist_idx == 1) {
      new_gk <- 0.0
    } else {
      v1 <- diagonal_zz + residual_variance / mixture_variances[selected_dist_idx]
      new_gk <- rnorm(1, mean = rhs / v1, sd = sqrt(residual_variance / v1))
      residuals <- residuals - snp_col * new_gk
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
#' Fortran: scale=(dble(included)*sum(g**2) + vara_ap*dfvara)/(dfvara+dble(included))
#' vara=rand_scaled_inverse_chi_square(dble(included)+dfvara,scale)
update_genetic_variance_mimic <- function(snp_counts_per_dist_cat, sum_sq_simple_effects, df_prior, vara_ap) {
  num_dists <- nrow(snp_counts_per_dist_cat)
  included <- sum(snp_counts_per_dist_cat[2:num_dists, ])
  total_g2 <- sum(sum_sq_simple_effects[2:num_dists, ])
  
  df_post <- included + df_prior
  # scale = (included * total_g2 + vara_ap * df_prior) / (df_prior + included)
  # vara = (df_post * scale) / rchisq(df_post) = (included * total_g2 + ...) / rchisq(df_post)
  
  ss_post <- included * total_g2 + vara_ap * df_prior
  if (df_post > 0) {
    return(ss_post / rchisq(1, df_post))
  } else {
    return(0.01) # Fallback
  }
}

update_mixture_proportions <- function(snp_counts_per_dist_cat, dirichlet_prior_counts) {
  num_categories <- ncol(snp_counts_per_dist_cat)
  num_dists      <- nrow(snp_counts_per_dist_cat)
  new_pi <- matrix(0, nrow = num_dists, ncol = num_categories)
  # Fortran: dirx=dble(snpindist(:,j))+delta; p(:,j)=rdirichlet(ndist,dirx)
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
  
  # --- Step 0: Initialization (Align with Fortran) ---
  # Fortran: mu=1.0, yadj=0.0
  # Note: yadj in Fortran is residuals, but computed as y - Xg - mu
  mu <- 1.0
  genetic_variance <- config$initial_genetic_variance
  residual_variance <- config$initial_residual_variance
  
  # Fortran: g=dsqrt(vara/(0.5*dble(nloci)))
  initial_g_val <- sqrt(genetic_variance / (0.5 * num_snps))
  snp_effects <- rep(initial_g_val, num_snps)
  
  # Fortran: p(1)=0.5, p(2:ndist)=1.0/gpin(2:ndist) ... normalized to 0.5
  mixture_proportions <- matrix(0, nrow = num_dists, ncol = num_categories)
  mixture_proportions[1, ] <- 0.5
  if (num_dists > 1) {
    inv_scales <- 1.0 / config$mixture_variance_scales[2:num_dists]
    inv_scales[is.infinite(inv_scales)] <- 0
    norm_factors <- 0.5 * inv_scales / sum(inv_scales)
    for (i in 2:num_dists) {
      mixture_proportions[i, ] <- norm_factors[i-1]
    }
  }
  
  # Calculate initial residuals: y - mu - Xg
  residuals <- phenotypes - mu - drop(genotype_matrix %*% snp_effects)
  
  snp_diagonal_xpx <- colSums(genotype_matrix^2)
  
  chain_intercept        <- numeric(config$num_iterations)
  chain_genetic_variance <- numeric(config$num_iterations)
  chain_residual_variance <- numeric(config$num_iterations)
  
  vara_ap <- if (config$genetic_variance_df == -2) 0 else genetic_variance
  vare_ap <- if (config$residual_variance_df == -2) 0 else residual_variance
  
  # --- MCMC Loop ---
  for (iteration in 1:config$num_iterations) {
    
    # 1. Sample Mu (Matches Fortran)
    mu_update <- update_intercept(residuals, mu, residual_variance, num_individuals)
    mu <- mu_update$value
    residuals <- mu_update$residuals
    
    # 2. Update SNP Effects
    mixture_variances      <- config$mixture_variance_scales * genetic_variance
    noise_to_signal_ratios <- residual_variance / mixture_variances
    
    snp_update <- update_snp_effects(
      residuals, snp_effects, genotype_matrix, snp_category_list,
      snp_diagonal_xpx, mixture_proportions, mixture_variances,
      residual_variance, noise_to_signal_ratios, config
    )
    
    snp_effects <- snp_update$snp_effects
    residuals   <- snp_update$residuals
    
    # 3. Update Genetic Variance (VCE=true)
    genetic_variance <- update_genetic_variance_mimic(snp_update$counts, snp_update$simple_ss, 
                                                       config$genetic_variance_df, vara_ap)
    
    # 4. Update Residual Variance (VCE=true)
    residual_variance <- update_residual_variance_mimic(residuals, num_individuals, 
                                                        config$residual_variance_df, vare_ap)
    
    # 5. Update Mixture Proportions
    mixture_proportions <- update_mixture_proportions(snp_update$counts, config$dirichlet_prior_counts)
    
    # Store Samples
    chain_intercept[iteration]        <- mu
    chain_genetic_variance[iteration] <- genetic_variance
    chain_residual_variance[iteration] <- residual_variance
    
    if (iteration %% 50 == 0) {
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
