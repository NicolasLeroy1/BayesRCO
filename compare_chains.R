args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript compare_chains.R <ref.hyp> <new_fortran.hyp> [c.hyp] [output.pdf]")
}

ref_file <- args[1]
new_file <- args[2]
c_file   <- if (length(args) >= 3 && !grepl(".pdf$", args[3])) args[3] else NULL
pdf_file <- if (length(args) >= 4) args[4] else if (length(args) == 3 && grepl(".pdf$", args[3])) args[3] else "verification_plots.pdf"

if (!require("ggplot2", quietly = TRUE)) {
  cat("Installing ggplot2...\n")
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  library(ggplot2)
}
if (!require("gridExtra", quietly = TRUE)) {
  cat("Installing gridExtra...\n")
  install.packages("gridExtra", repos = "http://cran.us.r-project.org")
  library(gridExtra)
}

read_hyp <- function(file, label) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  lines <- readLines(file)
  # Filter for numeric lines (start with digit or space+digit)
  numeric_lines <- lines[grepl("^[ \t]*[0-9]", lines)]
  if (length(numeric_lines) == 0) return(NULL)
  
  tmp <- tempfile()
  writeLines(numeric_lines, tmp)
  data <- read.table(tmp, header = FALSE)
  unlink(tmp)
  
  # Columns: 1=Iter, 2=Included, 3=Vara, 4=Vare
  df <- data.frame(
    Iteration = data[,1],
    Nsnp = data[,2],
    Vara = data[,3],
    Vare = data[,4],
    Source = label
  )
  # Calculate Heritability
  df$h2 <- df$Vara / (df$Vara + df$Vare)
  
  return(df)
}

cat("Reading reference hyp file:", ref_file, "\n")
ref_data <- read_hyp(ref_file, "Ref (Old)")

cat("Reading new Fortran hyp file:", new_file, "\n")
new_data <- read_hyp(new_file, "Ref (New)")

c_data <- NULL
if (!is.null(c_file)) {
  cat("Reading C hyp file:", c_file, "\n")
  c_data <- read_hyp(c_file, "C Version")
}

if (is.null(ref_data) || is.null(new_data)) {
  stop("Error: Failed to read data from required files.")
}

combined <- rbind(ref_data, new_data)
if (!is.null(c_data)) {
    combined <- rbind(combined, c_data)
}

# Create plots
pdf(pdf_file, width = 12, height = 10)

# Plots
p1 <- ggplot(combined, aes(x = Iteration, y = Vara, color = Source)) +
  geom_line(alpha = 0.7) + theme_minimal() + labs(title = "Trace: Vara")
p2 <- ggplot(combined, aes(x = Vara, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: Vara")
p3 <- ggplot(combined, aes(x = Iteration, y = Vare, color = Source)) +
  geom_line(alpha = 0.7) + theme_minimal() + labs(title = "Trace: Vare")
p4 <- ggplot(combined, aes(x = Vare, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: Vare")
p5 <- ggplot(combined, aes(x = Iteration, y = h2, color = Source)) +
  geom_line(alpha = 0.7) + theme_minimal() + labs(title = "Trace: h2")
p6 <- ggplot(combined, aes(x = h2, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: h2")
p7 <- ggplot(combined, aes(x = Iteration, y = Nsnp, color = Source)) +
  geom_line(alpha = 0.7) + theme_minimal() + labs(title = "Trace: Nsnp")
p8 <- ggplot(combined, aes(x = Nsnp, fill = Source)) +
  geom_density(alpha = 0.4) + theme_minimal() + labs(title = "Density: Nsnp")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
dev.off()
cat("Plots saved to", pdf_file, "\n")

# Statistical Tests
compare_dist <- function(v1, v2, name, label) {
  cat(sprintf("\nComparing %s vs %s:\n", name, label))
  m1 <- mean(v1, na.rm=TRUE); m2 <- mean(v2, na.rm=TRUE)
  s1 <- sd(v1, na.rm=TRUE);   s2 <- sd(v2, na.rm=TRUE)
  cat(sprintf("  Ref: mean=%.6e, sd=%.6e\n", m1, s1))
  cat(sprintf("  %s: mean=%.6e, sd=%.6e\n", label, m2, s2))
  
  if (s1 == 0 && s2 == 0) {
      if (abs(m1 - m2) < 1e-6) return(TRUE)
      return(FALSE)
  }
  ks <- ks.test(v1, v2)
  cat(sprintf("  KS test p-value: %.4f\n", ks$p.value))
  return(ks$p.value > 0.001)
}

fail <- FALSE
cat("\n--- Checking New Fortran ---\n")
if (!compare_dist(ref_data$Vara, new_data$Vara, "Vara", "New Fortran")) fail <- TRUE
if (!compare_dist(ref_data$Vare, new_data$Vare, "Vare", "New Fortran")) fail <- TRUE
if (!compare_dist(ref_data$h2,   new_data$h2,   "h2",   "New Fortran")) fail <- TRUE
if (!compare_dist(ref_data$Nsnp, new_data$Nsnp, "Nsnp", "New Fortran")) fail <- TRUE

if (!is.null(c_data)) {
    cat("\n--- Checking C Version ---\n")
    if (!compare_dist(ref_data$Vara, c_data$Vara, "Vara", "C Version")) fail <- TRUE
    if (!compare_dist(ref_data$Vare, c_data$Vare, "Vare", "C Version")) fail <- TRUE
    if (!compare_dist(ref_data$h2,   c_data$h2,   "h2",   "C Version")) fail <- TRUE
    if (!compare_dist(ref_data$Nsnp, c_data$Nsnp, "Nsnp", "C Version")) fail <- TRUE
}

if (!fail) {
  cat("\nSUCCESS: Verification passed.\n")
  quit(status = 0)
} else {
  cat("\nFAILURE: Significant differences detected.\n")
  quit(status = 1)
}
