#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript compare_models.R fortran.hyp c.hyp r.hyp [burnin]\n")
  q(status = 1)
}

fortran_file <- args[1]
c_file       <- args[2]
r_file       <- args[3]
burnin       <- if (length(args) >= 4) as.integer(args[4]) else 10

read_hyp <- function(file) {
  if (!file.exists(file)) return(NULL)
  data <- read.table(file, header = TRUE)
  if (nrow(data) <= burnin) return(NULL)
  return(data[(burnin + 1):nrow(data), ])
}

f_data <- read_hyp(fortran_file)
c_data <- read_hyp(c_file)
r_data <- read_hyp(r_file)

summarize <- function(name, data) {
  if (is.null(data)) {
    cat(sprintf("%-10s: FILE NOT FOUND or EMPTY\n", name))
    return()
  }
  
  params <- c("mu", "vara", "vare")
  cat(sprintf("\n--- %s summaries ---\n", name))
  cat(sprintf("%-6s %10s %10s %10s %10s\n", "Param", "Mean", "SD", "2.5%", "97.5%"))
  
  for (p in params) {
    if (p %in% names(data)) {
      val <- data[[p]]
      cat(sprintf("%-6s %10.4f %10.4f %10.4f %10.4f\n", 
                  p, mean(val), sd(val), quantile(val, 0.025), quantile(val, 0.975)))
    }
  }
}

summarize("Fortran", f_data)
summarize("C", c_data)
summarize("R", r_data)

# basic check
if (!is.null(f_data) && !is.null(c_data) && !is.null(r_data)) {
    cat("\nComparison Check:\n")
    # check vare means
    v_f <- mean(f_data$vare)
    v_c <- mean(c_data$vare)
    v_r <- mean(r_data$vare)
    
    diff_cf <- abs(v_c - v_f) / v_f
    diff_rf <- abs(v_r - v_f) / v_f
    
    cat(sprintf("Vare Rel Diff (C vs F): %.2f%%\n", diff_cf * 100))
    cat(sprintf("Vare Rel Diff (R vs F): %.2f%%\n", diff_rf * 100))
    
    if (diff_cf < 0.20 && diff_rf < 0.20) {
        cat("SUCCESS: All versions are statistically consistent.\n")
    } else {
        cat("WARNING: Significant discrepancy detected in posterior means!\n")
    }
}
