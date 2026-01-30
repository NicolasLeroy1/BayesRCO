#' Pure Statistical Kernel for BayesRCpi (Modular & Documented)
#' 
#' --- Total Prior Law ---
#' The model is: y = mu + Xg + e, with e ~ N(0, sigma_e^2 * I)
#' 
#' 1. Intercept: mu ~ Flat
#' 2. SNP Effects: g_k | d_k = d, sigma_a^2 ~ Normal(0, sigma_a^2 * scale_d)
#' 3. Distribution Assignment: d_k | category = c ~ Categorical(pi_c)
#' 4. Mixture Proportions: pi_c ~ Dirichlet(delta)
#' 5. Genetic Variance: sigma_a^2 ~ Inverse-Gamma(alpha_a, beta_a)
#' 6. Residual Variance: sigma_e^2 ~ Scaled-Inverse-Chi2(df_e, S_e)
#' -------------------------
#' 
#' This file implements a Gibbs sampler for Bayesian Multiple Regression
#' with Annotation-Specific Mixture Proportions (RCO).

# --- Statistical Helper Functions ---

# Sample from Inverse Gamma distribution
sample_inverse_gamma <- function(n, shape, scale) {
  1.0 / rgamma(n, shape, rate = scale)
}

# Sample from Dirichlet distribution
sample_dirichlet <- function(alpha_counts) {
  gamma_samples <- rgamma(length(alpha_counts), alpha_counts, 1.0)
  return(gamma_samples / sum(gamma_samples))
}

# --- Parameter Update Functions ---

#' Update Residual Variance (sigma_e^2)
#' @description 
#' Law: sigma_e^2 | y, mu, g ~ Scaled Inverse Chi-Square(n_ind + 3, RSS)
#' where RSS = sum( (y - mu - Xg)^2 )
update_residual_variance <- function(residuals, num_individuals) {
  sum_squared_residuals <- sum(residuals^2)
  return(sum_squared_residuals / rchisq(1, df = num_individuals + 3))
}

#' Update Intercept (mu)
#' @description 
#' Law: mu | y, g, sigma_e^2 ~ Normal( mean(y - Xg), sigma_e^2 / n_ind )
update_intercept <- function(residuals, intercept, residual_variance, num_individuals) {
  y_adj <- residuals + intercept
  
  intercept_mean <- mean(y_adj)
  intercept_sd   <- sqrt(residual_variance / num_individuals)
  new_intercept  <- rnorm(1, mean = intercept_mean, sd = intercept_sd)
  
  return(list(
    value = new_intercept,
    residuals = y_adj - new_intercept
  ))
}

#' Update SNP Effects (g_k) and Distribution Assignments (d_k)
#' @description 
#' Handles Overlapping Annotations (BayesRCpi):
#' 1. For SNP k, retrieve set of possible annotations S_k.
#' 2. Calculate Joint Probability P(A_k=a, C_k=c | Data) proportional to:
#'    Likelihood(Data | C_k=c) * pi_{a,c}
#' 3. Sample pair (a_new, c_new).
#' 
#' @param snp_category_list List where the k-th element is an integer vector of annotation indices for SNP k.
update_snp_effects <- function(residuals, snp_effects, genotype_matrix, snp_category_list, 
                               snp_diagonal_xpx, mixture_proportions, mixture_variances, 
                               residual_variance, noise_to_signal_ratios, config) {
  
  num_snps        <- ncol(genotype_matrix)
  num_individuals <- nrow(genotype_matrix)
  num_dists       <- length(mixture_variances)
  num_categories  <- ncol(mixture_proportions)
  
  # Tracks counts for pi update: rows=dist, cols=category
  snp_counts_per_dist_cat <- matrix(0, nrow = num_dists, ncol = num_categories)
  # Tracks weighted SS for var_g update (aggregated by component, but we track by cat too for consistency if needed later)
  # Actually update_genetic_variance sums over categories to get total per component.
  sum_sq_weighted_effects  <- matrix(0, nrow = num_dists, ncol = num_categories)
  
  shuffled_snps <- sample(1:num_snps)
  
  for (k in shuffled_snps) {
    diagonal_zz <- snp_diagonal_xpx[k]
    current_gk  <- snp_effects[k]
    
    # List of possible annotations for this SNP
    possible_annotations <- snp_category_list[[k]]
    
    snp_column  <- genotype_matrix[, k]
    
    # Remove current effect from residuals
    if (current_gk != 0.0) {
      residuals <- residuals + snp_column * current_gk
    }
    
    rhs_dot_product <- sum(residuals * snp_column)
    
    # --- Calculate Likelihoods for each Mixture Component (d) ---
    # L(d) depends only on d, not on annotation a
    log_likelihoods <- numeric(num_dists)
    
    # d=1 (Null)
    log_likelihoods[1] <- 0 # Ref level, actually usually we handle log prob differences.
    # But wait, the full form is:
    # L(d) \propto 1/sqrt(V_y) * exp( ... )
    # For d=1 (var=0), V_y = sigma_e2. LogDet = log(1) = 0. Quad = 0.
    
    cond_means <- numeric(num_dists)
    cond_vars  <- numeric(num_dists)
    
    for (d in 2:num_dists) {
      log_det_variance <- log(mixture_variances[d] * diagonal_zz / residual_variance + 1.0)
      conditioned_mean <- rhs_dot_product / (diagonal_zz + noise_to_signal_ratios[d])
      
      log_likelihoods[d] <- -0.5 * (log_det_variance - (rhs_dot_product * conditioned_mean / residual_variance))
      
      cond_means[d] <- conditioned_mean
      cond_vars[d]  <- residual_variance / (diagonal_zz + noise_to_signal_ratios[d])
    }
    
    # --- Calculate Joint Posteriors for (Annotation a, Component d) ---
    # P(a, d) \propto L(d) * pi_{a, d}
    # We flatten this into a single vector of length |Annotations| * |Components|
    # But we can just sample directly.
    
    num_opts <- length(possible_annotations) * num_dists
    log_probs_flat <- numeric(num_opts)
    
    # Map back structure
    # flat_index -> (annot_idx_in_list, dist_idx)
    
    flat_idx <- 1
    for (i in 1:length(possible_annotations)) {
      a <- possible_annotations[i]
      for (d in 1:num_dists) {
        log_probs_flat[flat_idx] <- log_likelihoods[d] + log(mixture_proportions[d, a])
        flat_idx <- flat_idx + 1
      }
    }
    
    # --- Sample (a, d) ---
    max_log_prob <- max(log_probs_flat)
    stable_probs <- exp(log_probs_flat - max_log_prob)
    sum_probs    <- sum(stable_probs)
    
    sampled_flat <- sample(1:num_opts, 1, prob = stable_probs / sum_probs)
    
    # Decode sampled index
    # (sampled_flat - 1) // num_dists gives annot index (0-based)
    # (sampled_flat - 1) %% num_dists gives dist index (0-based)
    
    sampled_annot_vec_idx <- (sampled_flat - 1) %/% num_dists + 1
    selected_annot <- possible_annotations[sampled_annot_vec_idx]
    
    selected_dist_idx0 <- (sampled_flat - 1) %% num_dists
    selected_dist <- selected_dist_idx0 + 1
    
    # --- Update Effect Size ---
    if (selected_dist == 1) {
      new_gk <- 0.0
    } else {
      new_gk <- rnorm(1, mean = cond_means[selected_dist], sd = sqrt(cond_vars[selected_dist]))
      residuals <- residuals - snp_column * new_gk
      
      sum_sq_weighted_effects[selected_dist, selected_annot] <- sum_sq_weighted_effects[selected_dist, selected_annot] + 
        (new_gk^2 / config$mixture_variance_scales[selected_dist])
    }
    
    snp_effects[k] <- new_gk
    snp_counts_per_dist_cat[selected_dist, selected_annot] <- snp_counts_per_dist_cat[selected_dist, selected_annot] + 1
  }
  
  return(list(
    snp_effects = snp_effects,
    residuals = residuals,
    counts = snp_counts_per_dist_cat,
    weighted_ss = sum_sq_weighted_effects
  ))
}

#' Update Genetic Variance (sigma_a^2)
#' @description 
#' Law: sigma_a^2 | g, d ~ Inverse-Gamma( shape_post, scale_post )
update_genetic_variance <- function(snp_counts_per_dist_cat, sum_sq_weighted_effects, config) {
  num_dists <- nrow(snp_counts_per_dist_cat)
  
  total_snps_included <- sum(snp_counts_per_dist_cat[2:num_dists, ])
  total_weighted_ss   <- sum(sum_sq_weighted_effects[2:num_dists, ])
  
  if (config$genetic_variance_df > 0.0) {
    ig_shape <- 0.5 * (total_snps_included + config$genetic_variance_df)
    ig_scale <- 0.5 * (config$initial_genetic_variance * config$genetic_variance_df + total_weighted_ss)
    return(sample_inverse_gamma(1, ig_shape, ig_scale))
  } else {
    ig_shape <- 0.5 * total_snps_included
    ig_scale <- 0.5 * total_weighted_ss
    if (ig_shape > 0.0) {
      return(sample_inverse_gamma(1, ig_shape, ig_scale))
    } else {
      return(config$initial_genetic_variance)
    }
  }
}

#' Update Mixture Proportions (pi_c)
#' @description 
#' Law: pi_c | d ~ Dirichlet( delta + counts_c )
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

run_bayesRCpi_mcmc <- function(phenotypes, genotype_matrix, snp_category_list, config) {
  num_individuals <- length(phenotypes)
  num_snps        <- ncol(genotype_matrix)
  num_dists       <- config$num_distributions
  
  # Find max category index to size matrices
  # Flatten list
  all_cats <- unlist(snp_category_list)
  if (length(all_cats) == 0) {
    num_categories <- 1 
  } else {
    num_categories <- max(all_cats)
  }
  
  intercept        <- mean(phenotypes)
  snp_effects      <- rep(0.0, num_snps)
  genetic_variance <- config$initial_genetic_variance
  residual_variance <- config$initial_residual_variance
  mixture_proportions <- matrix(1/num_dists, nrow = num_dists, ncol = num_categories)
  
  residuals <- phenotypes - intercept
  snp_diagonal_xpx <- colSums(genotype_matrix^2)
  
  chain_intercept        <- numeric(config$num_iterations)
  chain_genetic_variance <- numeric(config$num_iterations)
  chain_residual_variance <- numeric(config$num_iterations)
  
  for (iteration in 1:config$num_iterations) {
    residual_variance <- update_residual_variance(residuals, num_individuals)
    
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
    
    genetic_variance <- update_genetic_variance(snp_update$counts, snp_update$weighted_ss, config)
    mixture_proportions <- update_mixture_proportions(snp_update$counts, config$dirichlet_prior_counts)
    
    chain_intercept[iteration]        <- intercept
    chain_genetic_variance[iteration] <- genetic_variance
    chain_residual_variance[iteration] <- residual_variance
    
    if (iteration %% 10 == 0) {
      cat(sprintf("Iteration %d: Genetic Var=%.6f, Residual Var=%.6f, Polygenic SNPs=%d\n", 
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
