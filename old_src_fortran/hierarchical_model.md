# Hierarchical Bayesian Model (BayesRCO)

This document describes the statistical model implemented in the deduplicated Fortran code. The model is a hierarchical Bayesian linear regression model (BayesR) extended with annotations (BayesRC), estimating SNP effects using a mixture of normal distributions.

## 1. The Linear Model

The phenotype vector $\mathbf{y}$ of size $N$ is modeled as:

$$ \mathbf{y} = \mu \mathbf{1} + \mathbf{X}\mathbf{g} + \mathbf{e} $$

Where:
*   $\mathbf{y}$: Vector of phenotypes (centered/adjusted).
*   $\mu$: Global intercept (mean).
*   $\mathbf{1}$: Vector of ones.
*   $\mathbf{X}$: $N \times M$ matrix of centered and scaled genotypes.
*   $\mathbf{g}$: Vector of size $M$ containing genetic effects for each SNP.
*   $\mathbf{e}$: Vector of residual errors.

---

## 2. Distributions and Priors

### 2.1 Residual Error
The residuals are assumed to differ normally with variance $\sigma^2_e$:
$$ \mathbf{e} \sim \mathcal{N}(0, \sigma^2_e \mathbf{I}) $$

**Prior for $\sigma^2_e$**: Scaled Inverse Chi-Square distribution.
$$ \sigma^2_e \sim \chi^{-2}(\nu_e, S_e) $$
Where $\nu_e$ are degrees of freedom and $S_e$ is the scale parameter.

### 2.2 Global Mean
The intercept $\mu$ has an implicit flat prior.

### 2.3 Genetic Effects (Mixture Model)
Each SNP effect $g_j$ is drawn from a mixture of $K$ (typically 4) Normal distributions:

$$ g_j \sim \sum_{k=1}^{K} \pi_{k,c(j)} \mathcal{N}(0, \sigma^2_k) $$

Where:
*   $c(j)$: Category annotation for SNP $j$.
*   $\pi_{k,c}$: Probability (proportion) that a SNP in category $c$ belongs to mixture component $k$.
*   $\sigma^2_k$: Variance of component $k$, defined relative to the total genetic variance $\sigma^2_a$:
    *   Component 1 (Null): $\sigma^2_1 = 0$
    *   Component $k > 1$: $\sigma^2_k = \gamma_k \sigma^2_a$ (where $\gamma_k$ are fixed constants, e.g., $0.0001, 0.001, 0.01$).

### 2.4 Genetic Variance ($\sigma^2_a$)
**Prior**: Scaled Inverse Chi-Squared.
$$ \sigma^2_a \sim \chi^{-2}(\nu_a, S_a) $$

### 2.5 Mixture Proportions ($\pi$)
The vector of proportions $\mathbf{\pi}_c = (\pi_{1,c}, \dots, \pi_{K,c})$ for category $c$ follows a Dirichlet distribution:
$$ \mathbf{\pi}_c \sim \text{Dirichlet}(\alpha + \mathbf{n}_c) $$
Where $\mathbf{n}_c$ is the vector of counts of SNPs in category $c$ assigned to each component.

---

## 3. Full Conditional Posteriors (Gibbs Sampling)

The parameters are updated iteratively using a Gibbs sampler.

### 3.1 Sampling $\sigma^2_e$ (Residual Variance)
$$ \sigma^2_e | \mathbf{y}, \mathbf{g}, \dots \sim \chi^{-2}(N + \nu_e, \mathbf{e}'\mathbf{e} + S_e) $$
In the code:
```fortran
vare = dot_product(yadj,yadj) / rand_chi_square(nnind + dfvare)
```

### 3.2 Sampling $\mu$ (Mean)
$$ \mu | \mathbf{y}, \sigma^2_e \sim \mathcal{N}\left( \frac{\sum (y_i - \mathbf{x}_i \mathbf{g})}{N}, \frac{\sigma^2_e}{N} \right) $$

### 3.3 Sampling $\sigma^2_a$ (Genetic Variance)
$$ \sigma^2_a | \mathbf{g}, \text{allocation} \sim \chi^{-2}(N_{inc} + \nu_a, \sum g^2 + S_a) $$
Where $N_{inc}$ is the number of SNPs included in the model (non-zero effect).

### 3.4 Sampling SNP Effects $g_j$ and Component Allocation $k_j$

For each SNP $j$:

1.  **Calculate Marginal Likelihood** for each component $k$:
    Let $\mathbf{y}^* = \mathbf{y} - \mathbf{X}_{-j}\mathbf{g}_{-j}$ (residuals without SNP $j$).
    
    $$ L_k \propto \frac{1}{\sqrt{\mathbf{z}_j'\mathbf{z}_j \sigma^2_k + \sigma^2_e}} \exp \left( \frac{(\mathbf{z}_j'\mathbf{y}^*)^2}{2} \cdot \frac{\sigma^2_k}{\sigma^2_e (\mathbf{z}_j'\mathbf{z}_j \sigma^2_k + \sigma^2_e)} \right) $$

2.  **Calculate Posterior Probability**:
    $$ P(k_j=k | \text{data}) = \frac{\pi_k L_k}{\sum_{l} \pi_l L_l} $$

3.  **Sample Allocation**:
    Draw component $k$ from the categorical distribution defined by $P(k | \dots)$.

4.  **Sample Effect**:
    *   If $k=1$: $g_j = 0$
    *   If $k>1$: Draw $g_j$ from conditional normal:
        $$ g_j \sim \mathcal{N}\left( \frac{\mathbf{z}_j'\mathbf{y}^*}{\mathbf{z}_j'\mathbf{z}_j + \frac{\sigma^2_e}{\sigma^2_k}}, \frac{\sigma^2_e}{\mathbf{z}_j'\mathbf{z}_j + \frac{\sigma^2_e}{\sigma^2_k}} \right) $$

### 3.5 Sampling Mixture Proportions $\pi$
For each category $c$:
$$ \mathbf{\pi}_c \sim \text{Dirichlet}(\alpha_1 + n_{1,c}, \dots, \alpha_K + n_{K,c}) $$
Where $n_{k,c}$ is the number of SNPs in category $c$ assigned to component $k$ in the current iteration.
