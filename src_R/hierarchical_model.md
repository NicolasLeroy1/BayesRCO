# Hierarchical Bayesian Model (R Implementation)

This document describes the statistical model implemented in the R code (`src_R/bayesRCpi.R`). The model differs slightly in parameterization from the Fortran version but represents the same underlying statistical framework.

## 1. The Linear Model

The phenotype vector $\mathbf{y}$ is modeled as:

$$ \mathbf{y} = \mu \mathbf{1} + \mathbf{X}\mathbf{g} + \mathbf{e} $$

## 2. Distributions and Priors

### 2.1 Residual Error
$$ \mathbf{e} \sim \mathcal{N}(0, \sigma^2_e \mathbf{I}) $$

**Prior**: Scaled Inverse Chi-Square (implemented via Chi-Square sampling).
$$ \sigma^2_e \sim \chi^{-2}(N + 3, RSS) $$

### 2.2 Global Mean
**Prior**: Flat (Improper Uniform).
**Update**: Normal distribution conditional on $\sigma^2_e$.

### 2.3 Genetic Effects (Mixture Model)
$$ g_j \sim \sum_{k=1}^{K} \pi_{k,c(j)} \mathcal{N}(0, \sigma^2_k) $$

*   Component 1: $\sigma^2_1 = 0$
*   Component $k > 1$: $\sigma^2_k = \gamma_k \sigma^2_a$

### 2.4 Genetic Variance ($\sigma^2_a$)
**Prior**: Inverse-Gamma $(\alpha_a, \beta_a)$.
*   Note: The Fortran code uses Scaled Inverse Chi-Square. These are related, but the R implementation explicitly uses the Gamma parameterization `rgamma` then inverts.

$$ \sigma^2_a \sim IG\left( \frac{N_{inc} + \nu_a}{2}, \frac{S_{prior} + \sum g^2_{weighted}}{2} \right) $$

where $\sum g^2_{weighted} = \sum_{k>1} \frac{\sum g_{j \in k}^2}{\gamma_k}$.

### 2.5 Mixture Proportions ($\pi$)
**Prior**: Dirichlet Distribution.
$$ \mathbf{\pi}_c \sim \text{Dirichlet}(\alpha + \mathbf{n}_c) $$

## 3. Full Conditional Posteriors (Gibbs Sampling)

### 3.1 Residual Variance ($\sigma^2_e$)
$$ \sigma^2_e | \dots \sim \frac{RSS}{\chi^2_{df}} $$
where $df = N + 3$.

### 3.2 Intercept ($\mu$)
$$ \mu | \dots \sim \mathcal{N}(\bar{y}_{adj}, \frac{\sigma^2_e}{N}) $$

### 3.3 Genetic Variance ($\sigma^2_a$)
Updated using the weighted sum of squares of effects.
$$ \sigma^2_a | \dots \sim IG \left( \frac{\nu_{prior} + N_{inc}}{2}, \frac{S_{prior}\nu_{prior} + \sum (\frac{g^2}{\gamma})}{2} \right) $$

### 3.4 SNP Effects and Allocation (BayesRCpi)
The R implementation supports **Overlapping Annotations**.
For a SNP $j$ with a set of possible annotations $S_j$:

1.  **Joint Posterior Probability**:
    For each pair $(a, d)$ where $a \in S_j$ (annotation) and $d \in \{1 \dots K\}$ (mixture component):
    $$ P(a, d | \text{data}) \propto \mathcal{L}(d) \times \pi_{d,a} $$
    
    where $\mathcal{L}(d) = \mathcal{N}(y^* | 0, \sigma^2_e + z'z\sigma^2_d)$ is the marginal likelihood of the residual given component $d$.

2.  **Sampling**:
    Sample the pair $(a, d)$ directly from the flattened probability vector.

3.  **Effect Update**:
    Sample $g_j$ based on the chosen component $d$.

### 3.5 Mixture Proportions ($\pi$)
Standard Dirichlet update based on counts of SNPs assigned to each component within each annotation category.
