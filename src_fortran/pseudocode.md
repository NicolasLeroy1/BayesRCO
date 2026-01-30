# Pseudocode: BayesRCpi (Fortran Implementation)

This version implements the **BayesRCpi** algorithm with joint sampling of annotations and mixture distributions.

## Algorithm Overview

```
1. Initialize:
   - mu = 1.0
   - vara = initial_vara (scaled by vary if dfvara < -2)
   - vare = initial_vare
   - p (mixture proportions) = 0.5 for Null, 0.5 split among others (per category)
   - g (SNP effects) = sqrt(vara / (0.5 * nloci))

2. Main MCMC Loop (numit iterations):
   a. Update Residual Variance (vare):
      - Sample from Inverse-Gamma using SSR (Sum of Squared Residuals)
   
   b. Update Intercept (mu):
      - Sample from Normal using Mean(residuals + mu)
   
   c. Update SNP Effects (g) - Random order:
      - For each SNP k:
         i.   Remove current effect: yadj = yadj + Xk * gk
         ii.  Identify active categories (C[k, i] == 1)
         iii. For each (category c, distribution d) pair:
              - Compute log-likelihood L(d)
              - Compute log-probability: P(c, d) = L(d) + log(pi[d, c])
         iv.  Sample (chosen_cat, chosen_dist) via log-sum-exp & categorical sampling
         v.   If chosen_dist is NOT Null:
              - Sample new gk from Normal posterior
              - Update residuals: yadj = yadj - Xk * gk
         vi.  Else:
              - gk = 0

   d. Update Genetic Variance (vara):
      - Sample from Inverse-Gamma using sum(gk^2 / scale_k)
   
   e. Update Mixture Proportions (pi):
      - For each category c:
         - Sample new pi vector from Dirichlet(prior + counts[c])
```

## Key Features
- **Joint Sampling**: Samples the annotation and distribution class in a single step.
- **Inverse-Gamma Priors**: Uses Inverse-Gamma for both variance components.
- **RCO Support**: Mixture proportions $\pi$ are specific to each annotation category.
