# Pseudocode: BayesR Mimic (R Implementation)

This version is specifically designed to replicate the statistical behavior of the **Old Fortran** BayesRCO implementation.

## Algorithm Overview

```
1. Initialize (Mimic Old Fortran):
   - mu = 1.0
   - genetic_variance (Va) = initial_val
   - mixture_proportions (pi): 0.5 for Null, 0.5 split among others (mimic logic)
   - initial g = sqrt(Va / (0.5 * nloci))
   - residuals = phenotypes - mu - (X %*% g)
   - vara_ap = Va, vare_ap = initial_Ve

2. Main MCMC Loop:
   a. Update Intercept (mu):
      - mu = Normal(Mean(residuals + mu), sqrt(Ve/N))
      - residuals adjusted
   
   b. Update SNP Effects (g) - Shuffled:
      - For each SNP k:
         i.   Remove current effect: residuals += Xk * gk
         ii.  Calculate current RHS = dot_product(Xk, residuals)
         iii. Calculate joint probabilities (for sampling annotation & distribution)
         iv.  Sample (selected_dist, selected_annot)
         v.   If selected_dist > 1:
              - Sample new_gk from Normal(rhs / v1, sqrt(Ve/v1))
              - residuals -= Xk * new_gk
         vi.  Else:
              - new_gk = 0

   c. Update Genetic Variance (Va) - Scaled Inverse Chi-square:
      - df_post = included + df_prior
      - ss_post = (included * sum(g^2)) + (vara_ap * df_prior)
      - Va = ss_post / rchisq(df_post)
   
   d. Update Residual Variance (Ve) - Scaled Inverse Chi-square:
      - df_post = N + df_prior
      - ss_post = sum(residuals^2) + (vare_ap * df_prior)
      - Ve = ss_post / rchisq(df_post)

   e. Update Mixture Proportions (pi):
      - Update per category using Dirichlet distribution.
```

## Key Features
- **Statistical Mimicry**: Uses Scaled Inverse Chi-square distributions for variance components instead of direct Inverse-Gamma, matching legacy Fortran behavior.
- **Old-Style Initialization**: Follows exactly how the old Fortran initialized $g$ and $\pi$.
- **Validation Tool**: Used to ensure that modern refactors haven't drifted from the original statistical model.
