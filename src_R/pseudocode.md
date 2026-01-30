# Pseudocode: BayesRCpi (R Implementation)

This version provides a clean, modular R implementation of the BayesRCpi algorithm.

## Algorithm Overview

```
1. Initialize:
   - intercept = mean(phenotypes)
   - genetic_variance = initial_val
   - mixture_proportions = 1/ndist (per category)
   - residuals = phenotypes - intercept
   - snp_diagonal_xpx = colSums(genotype_matrix^2)

2. Main MCMC Loop (iterations):
   a. Update Residual Variance (Ve):
      - SSR = sum(residuals^2)
      - Ve = SSR / rchisq(n_ind + 3)
   
   b. Update Intercept (mu):
      - y_adj = residuals + intercept
      - mu = rnorm(1, mean(y_adj), sqrt(Ve/n_ind))
      - residuals = y_adj - mu
   
   c. Update SNP Effects (g) - Random Order:
      - For each SNP k:
         i.   Remove current effect: residuals += Xk * g[k]
         ii.  Get possible categories S_k for SNP k
         iii. Vectorized Likelihood Calculation:
              - Compute log-likelihoods L(d) for all d in 1:ndist
         iv.  Joint Sampling:
              - probs_matrix[d, a] = L(d) + log(pi[d, a])
              - Sample (selected_dist, selected_annot)
         v.   If selected_dist > 1:
              - Sample new_g from Normal posterior
              - residuals -= Xk * new_g
         vi.  Else:
              - new_g = 0

   d. Update Genetic Variance (Va):
      - Sample from Inverse-Gamma using sum(g^2 / scale_d)
   
   e. Update Mixture Proportions (pi):
      - For each category c:
         - Sample from Dirichlet(delta + count[c])
```

## Key Features
- **Modular Design**: Uses separate functions for each parameter update.
- **R-Idiomatic**: Leverages `crossprod` and vectorized operations where possible.
- **Joint Sampling**: Correctly handles multiple/overlapping annotations via joint distribution/annotation sampling.
