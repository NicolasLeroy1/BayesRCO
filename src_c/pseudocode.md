# Pseudocode: BayesRCpi (C Implementation)

This version is a high-performance C port of the BayesRCpi algorithm.

## Algorithm Overview

```
1. Initialize:
   - mu = 1.0, vara = config->vara, vare = config->vare
   - g = 1.0 (approximated initially)
   - Precompute XX (xpx) for each SNP
   - Residuals (yadj) = phenotypes - mu

2. Main MCMC Loop (config->numit iterations):
   a. Update Residual Variance (vare):
      - Sample from Chi-square: vare = SumSquaredResiduals / rchisq(nt + 3)
   
   b. Update Intercept (mu):
      - mu = Normal(Mean(yadj + mu), sqrt(vare/nt))
      - Adjust residuals: yadj -= mu
   
   c. Update SNP Effects (g) - Shuffled order:
      - For each SNP k:
         i.   Calculate Current RHS: rhs = dot_product(yadj, Xk) + xpx[k] * g[k]
         ii.  Identify Active Categories (C[k, j] == 1)
         iii. For each (annotation a, distribution d) pair:
              - Compute log-likelihood L(d) using precomputed ratios
              - Compute log-probability: P(a, d) = L(d) + log(pi[d, a])
         iv.  Sample (chosen_a, chosen_d) using categorical sampling
         v.   If chosen_d is Null (0):
              - new_gk = 0
         vi.  Else:
              - Sample new_gk from Normal(rhs / denom, sqrt(vare / denom))
         vii. Efficiency Step: Update residuals yadj ONLY if gk changed:
              - yadj += Xk * (old_gk - new_gk)

   d. Update Genetic Variance (vara):
      - Sample from Inverse-Gamma using accumulated sum_sq_weighted_effects
   
   e. Update Mixture Proportions (pi):
      - For each category j:
         - Sample from Dirichlet(delta + count[j])
```

## Key Features
- **Joint Sampling**: Samples annotation and distribution class jointly.
- **Efficient Residual Update**: Only modifies the residual vector when a SNP effect changes.
- **Precomputed ZZ**: Precomputes diagonal of X'X for speed.
