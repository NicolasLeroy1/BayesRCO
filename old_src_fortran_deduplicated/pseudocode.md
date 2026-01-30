# Pseudocode: Deduplicated BayesRCO (Refactored Legacy Fortran)

This version is a refactored and modularized version of the **Old BayesRCO** code. It maintains exact algorithmic equivalence with the original legacy version while improving code structure.

## Core Logic (Run MCMC)

```
1. Initialize Common Parameters (mcmc_init_common)
2. Select Strategy:
   a. If Mixture: 
      - Initialize Mixture (mcmc_init_mixture)
      - Run MCMC Loop with Mixture Strategy
   b. Else If Additive:
      - Initialize Additive (mcmc_init_additive)
      - Run MCMC Loop with Additive Strategy
   c. Else If BayesCpi:
      - Initialize BayesCpi (mcmc_init_bayesCpi)
      - Run MCMC Loop with BayesCpi Strategy
3. Finalize and Save Results (mcmc_finalize)
```

## MCMC Loop (General Template)

```
For each iteration (1 to numit):
   1. Sample Residual Variance (Ve)
   2. Sample Intercept (mu)
   3. Update Variance Ratios (Ve/Vj)
   4. Permute SNP order
   5. Update SNP Effects (Strategy-specific):
      - Strategy 1 (Mixture): Sample annotation 'a' then sample distribution 'd'
      - Strategy 2 (Additive Loci-Major): Sample for each category per SNP
      - Strategy 3 (Additive Category-Major): Loop over categories then SNPs
   6. Compute Iteration Stats (Va, pi)
   7. Save Samples (if thin/burnin conditions met)
```

## Key Differences from Legacy Version
- **Modularity**: Logic is split into `mod_mcmc_core`, `mod_mcmc_helpers`, and `mod_mcmc_strategies`.
- **Maintainability**: Reduced use of global variables in favor of structured data (where possible in Fortran).
- **Functional Equivalence**: Uses the same Scaled Inverse Chi-square and two-step sampling logic as the `old_src_fortran` version.
