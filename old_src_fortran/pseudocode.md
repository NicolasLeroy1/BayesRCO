# Pseudocode: Old BayesRCO (Legacy Fortran)

This is the original implementation, featuring several sampling "strategies" and using legacy Fortran styles (global variables, etc.).

## Mixture Strategy (Standard BayesRCO)

```
1. Initialize:
   - mu = 1.0, g = sqrt(Va / (0.5 * nloci))
   - Mixture proportions (pi): 0.5 for Null, 0.5 among others.
   - Initial residuals: yadj = why - Xg - mu

2. Main MCMC Loop:
   a. Update Residual Variance (Ve) - Scaled Inverse Chi-square:
      - Ve = (sum(yadj^2) + vare_ap * dfvare) / rchisq(N + dfvare)
   
   b. Update Intercept (mu):
      - mu = Normal(Mean(yadj + mu), sqrt(Ve/N))
      - yadj adjusted
   
   c. Update SNP Effects (g) - Two-Step Sampling:
      - For each SNP k:
         i.   Remove effect: yadj += Xk * gk
         ii.  Sample Annotation (Two-step logic):
              - If multiple annotations (nannot > 1):
                - Calculate probabilities for each annotation c.
                - Sample annotation 'a' for this iteration.
         iii. Sample Distribution Class (conditional on 'a'):
              - Compute log-likelihoods L(d).
              - Sample distribution class 'indistflag'.
         iv.  If indistflag > 1:
              - Sample new_gk from Normal posterior.
              - yadj -= Xk * new_gk
         v.   Else:
              - new_gk = 0
         vi.  g[k] = new_gk

   d. Update Genetic Variance (Va) - Scaled Inverse Chi-square:
      - Va = (included * sum(g^2) + vara_ap * dfvara) / rchisq(included + dfvara)

   e. Update Mixture Proportions (pi):
      - Update for each category using Dirichlet distribution.
```

## Other Strategies
- **Additive**: Iterates through each category for each SNP, samples distributions independently, and sums the effects.
- **BayesCpi**: Category-major loop where mixture proportions/distributions are updated within categories.

## Key Features
- **Two-Step Sampling**: For SNPs with multiple categories, it first samples which category to "assign" the SNP to for that iteration, then samples the effect size distribution. (Modern versions use joint sampling).
- **Scaled Inverse Chi-square**: Uses this distribution for variance parameter updates.
