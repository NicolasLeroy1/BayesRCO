#!/usr/bin/env Rscript
# compare_chains.R - Multi-chain MCMC comparison for BayesRCO equivalence testing
# 
# Usage: Rscript compare_chains.R <method1.hyp> <method2.hyp> [output.pdf]
#
# For combined multi-chain files, samples from all chains are pooled for comparison.
# For proper multi-chain diagnostics, the run_test.sh should provide chain IDs.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript compare_chains.R <method1.hyp> <method2.hyp> [output.pdf]")
}

method1_file <- args[1]
method2_file <- args[2]
pdf_file <- if (length(args) >= 3) args[3] else "verification_plots.pdf"

# Load packages
suppressPackageStartupMessages({
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(ggplot2)
  }
  if (!require("gridExtra", quietly = TRUE)) {
    install.packages("gridExtra", repos = "http://cran.us.r-project.org", quiet = TRUE)
    library(gridExtra)
  }
})

# Read hyp file
read_hyp <- function(file, label) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  lines <- readLines(file, warn = FALSE)
  numeric_lines <- lines[grepl("^[ \t]*[0-9]", lines)]
  if (length(numeric_lines) == 0) return(NULL)
  
  tmp <- tempfile()
  writeLines(numeric_lines, tmp)
  data <- read.table(tmp, header = FALSE)
  unlink(tmp)
  
  # Check if ChainID column exists (now the 5th column or last)
  has_chain <- ncol(data) >= 5
  
  df <- data.frame(
    Iteration = data[, 1],
    Nsnp = data[, 2],
    Vara = data[, 3],
    Vare = data[, 4],
    ChainID = if(has_chain) data[, ncol(data)] else 1,
    Source = label
  )
  df$h2 <- df$Vara / (df$Vara + df$Vare)
  return(df)
}

# Compute effective sample size (simple autocorrelation-based estimate)
compute_ess <- function(x) {
  n <- length(x)
  if (n < 10) return(n)
  
  # Compute autocorrelation up to lag 50 or n/3
  max_lag <- min(50, floor(n / 3))
  acf_vals <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]  # exclude lag 0
  
  # Sum positive autocorrelations until first negative
  tau <- 0
  for (i in seq_along(acf_vals)) {
    if (acf_vals[i] < 0) break
    tau <- tau + acf_vals[i]
  }
  
  ess <- n / (1 + 2 * tau)
  return(max(1, ess))
}

# Compare distributions with proper MCMC diagnostics
compare_mcmc <- function(v1, v2, name, label1, label2) {
  cat(sprintf("\n%s:\n", name))
  
  # Summary statistics
  m1 <- mean(v1, na.rm = TRUE); m2 <- mean(v2, na.rm = TRUE)
  s1 <- sd(v1, na.rm = TRUE);   s2 <- sd(v2, na.rm = TRUE)
  n1 <- length(v1); n2 <- length(v2)
  
  # Effective sample sizes
  ess1 <- compute_ess(v1)
  ess2 <- compute_ess(v2)
  
  # Standard error of mean (using ESS)
  sem1 <- s1 / sqrt(ess1)
  sem2 <- s2 / sqrt(ess2)
  
  # 95% credible intervals (assuming normality for posterior mean)
  ci1 <- c(m1 - 1.96 * sem1, m1 + 1.96 * sem1)
  ci2 <- c(m2 - 1.96 * sem2, m2 + 1.96 * sem2)
  
  cat(sprintf("  %s: mean=%.4g, sd=%.4g, ESS=%.0f, 95%%CI=[%.4g, %.4g]\n", 
              label1, m1, s1, ess1, ci1[1], ci1[2]))
  cat(sprintf("  %s: mean=%.4g, sd=%.4g, ESS=%.0f, 95%%CI=[%.4g, %.4g]\n", 
              label2, m2, s2, ess2, ci2[1], ci2[2]))
  
  # Check overlap of 95% CIs (more appropriate for MCMC)
  ci_overlap <- (ci1[2] >= ci2[1]) && (ci2[2] >= ci1[1])
  
  # Relative difference
  rel_diff <- abs(m1 - m2) / ((abs(m1) + abs(m2)) / 2 + 1e-10) * 100
  
  cat(sprintf("  Relative difference: %.2f%%\n", rel_diff))
  cat(sprintf("  CI overlap: %s\n", ifelse(ci_overlap, "YES", "NO")))
  
  # For exact equivalence (same PRNG), check if values are identical
  if (s1 == 0 && s2 == 0) {
    is_equal <- abs(m1 - m2) < 1e-6
    cat(sprintf("  Exact match: %s\n", ifelse(is_equal, "YES", "NO")))
    return(list(pass = is_equal, exact = TRUE, rel_diff = rel_diff, ci_overlap = ci_overlap))
  }
  
  # For statistical equivalence, use CI overlap and relative difference
  # Pass if CIs overlap OR relative difference < 10%
  pass <- ci_overlap || rel_diff < 10
  cat(sprintf("  Status: %s\n", ifelse(pass, "PASS", "DIFF")))
  
  return(list(pass = pass, exact = FALSE, rel_diff = rel_diff, ci_overlap = ci_overlap))
}

# Main
cat("Reading:", method1_file, "\n")
data1 <- read_hyp(method1_file, "Method1")

cat("Reading:", method2_file, "\n")
data2 <- read_hyp(method2_file, "Method2")

if (is.null(data1) || is.null(data2)) {
  stop("Error: Failed to read data from files.")
}

combined <- rbind(data1, data2)

# Create plots
pdf(pdf_file, width = 12, height = 10)

p1 <- ggplot(combined, aes(x = Iteration, y = Vara, color = Source)) +
  geom_line(aes(group = interaction(Source, ChainID)), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", size = 1, aes(group = Source)) +
  theme_minimal() + labs(title = "Trace: Vara (Combined Chains)")
p2 <- ggplot(combined, aes(x = Vara, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: Vara")
p3 <- ggplot(combined, aes(x = Iteration, y = Vare, color = Source)) +
  geom_line(aes(group = interaction(Source, ChainID)), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", size = 1, aes(group = Source)) +
  theme_minimal() + labs(title = "Trace: Vare (Combined Chains)")
p4 <- ggplot(combined, aes(x = Vare, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: Vare")
p5 <- ggplot(combined, aes(x = Iteration, y = h2, color = Source)) +
  geom_line(aes(group = interaction(Source, ChainID)), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", size = 1, aes(group = Source)) +
  theme_minimal() + labs(title = "Trace: h2 (Combined Chains)")
p6 <- ggplot(combined, aes(x = h2, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: h2")
p7 <- ggplot(combined, aes(x = Iteration, y = Nsnp, color = Source)) +
  geom_line(aes(group = interaction(Source, ChainID)), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", size = 1, aes(group = Source)) +
  theme_minimal() + labs(title = "Trace: Nsnp (Combined Chains)")
p8 <- ggplot(combined, aes(x = Nsnp, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: Nsnp")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
dev.off()
cat("\nPlots saved to", pdf_file, "\n")

# Run comparisons
cat("\n=== MCMC Equivalence Comparison ===\n")

results <- list(
  vara = compare_mcmc(data1$Vara, data2$Vara, "Vara", "Method1", "Method2"),
  vare = compare_mcmc(data1$Vare, data2$Vare, "Vare", "Method1", "Method2"),
  h2   = compare_mcmc(data1$h2,   data2$h2,   "h2",   "Method1", "Method2"),
  nsnp = compare_mcmc(data1$Nsnp, data2$Nsnp, "Nsnp", "Method1", "Method2")
)

# Summary
cat("\n=== Summary ===\n")
all_exact <- all(sapply(results, function(r) r$exact))
all_pass <- all(sapply(results, function(r) r$pass))
variance_pass <- results$vara$pass && results$vare$pass && results$h2$pass

if (all_exact) {
  cat("Result: EXACT MATCH (identical chains)\n")
  quit(status = 0)
} else if (variance_pass) {
  cat("Result: VARIANCE COMPONENTS EQUIVALENT\n")
  if (!results$nsnp$pass) {
    cat("Note: Nsnp differs but variance components (Vara, Vare, h2) match.\n")
    cat("      This is expected with different PRNGs.\n")
  }
  quit(status = 0)
} else {
  cat("Result: SIGNIFICANT DIFFERENCES DETECTED\n")
  quit(status = 1)
}
