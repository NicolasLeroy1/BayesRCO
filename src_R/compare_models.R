#!/usr/bin/env Rscript

# Check for ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for visualization. Please install it.")
}
library(ggplot2)
library(grid)      # For layout
library(patchwork) # For combining plots
library(ggridges)  # For ridge plots

# --- Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript compare_models.R <file1.hyp> <file2.hyp> ... <burnin> [output.pdf]\n")
  q(status = 1)
}

# The last argument might be the PDF
last_arg <- args[length(args)]
if (grepl("\\.pdf$", last_arg, ignore.case = TRUE)) {
  pdf_path <- last_arg
  burnin_arg <- args[length(args) - 1]
  file_args <- args[1:(length(args) - 2)]
} else {
  pdf_path <- "similarity_visualization.pdf"
  burnin_arg <- last_arg
  file_args <- args[1:(length(args) - 1)]
}

burnin <- as.integer(burnin_arg)
if (is.na(burnin)) {
  stop("Burnin must be an integer.")
}

# Create named list of files
# Try to guess label from filename if possible, or just use indices
files <- list()
for (f in file_args) {
    label <- basename(f)
    label <- sub("\\.hyp$", "", label)
    # Nicer labels
    label <- sub("toy_", "", label)
    label <- sub("_complex", " (Complex)", label)
    label <- sub("_overl", " (Overlap)", label)
    label <- gsub("_", " ", label)
    # Capitalize
    label <- paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
    files[[label]] <- f
}

# --- Functions ---

read_hyp <- function(name, file) {
  if (!file.exists(file)) {
    cat(sprintf("WARNING: File not found for %s: %s\n", name, file))
    return(NULL)
  }
  
  # Strict reading: header=TRUE, no fill=TRUE
  data <- tryCatch(
    read.table(file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      cat(sprintf("ERROR: Failed to read %s (%s): %s\n", name, file, e$message))
      return(NULL)
    }
  )
  
  if (is.null(data)) return(NULL)
  
  # Normalize column names
  names(data) <- tolower(names(data))
  names(data)[names(data) == "va"] <- "vara"
  names(data)[names(data) == "ve"] <- "vare"
  
  needed <- c("mu", "vara", "vare")
  if (!all(needed %in% names(data))) {
    cat(sprintf("WARNING: %s missing required columns. Found: %s\n", name, paste(names(data), collapse=",")))
    return(NULL)
  }
  
  if (nrow(data) <= burnin) {
    cat(sprintf("WARNING: %s has insufficient rows (%d) for burnin (%d)\n", name, nrow(data), burnin))
    return(NULL)
  }
  
  final_data <- data[(burnin + 1):nrow(data), ]
  final_data$Model <- name
  return(final_data)
}

# --- Load Data ---
all_data_list <- lapply(names(files), function(n) read_hyp(n, files[[n]]))
names(all_data_list) <- names(files)
# Filter out NULLs
valid_data_list <- all_data_list[!sapply(all_data_list, is.null)]

if (length(valid_data_list) == 0) {
  stop("No valid data files loaded. Exiting.")
}

# Normalize columns for all valid data immediately
common_cols <- c("Model", "mu", "vara", "vare")
valid_data_list <- lapply(valid_data_list, function(d) {
  d_sub <- d[, common_cols]
  d_sub$Iteration_Rel <- 1:nrow(d_sub)
  d_sub
})

combined_df <- do.call(rbind, valid_data_list)
# Factor ordering for consistent plotting
combined_df$Model <- factor(combined_df$Model, levels = names(files))

params_to_check <- c("vare", "vara", "mu")
fail_verification <- FALSE

# --- Visualization ---
cat(sprintf("\nGenerating visualization to %s...\n", pdf_path))
pdf(pdf_path, width = 14, height = 12)

# Color Palette - Ensure keys match labels generated from files (line 48)
my_colors <- c(
  "Fortran"      = "#1f77b4", # New Fortran
  "Old fortran"  = "#d62728", 
  "Dedup"        = "#2ca02c", # Deduplicated
  "C"            = "#9467bd", 
  "R"            = "#8c564b",
  "R mimic"      = "#ff7f0e",
  # Keep legacy names just in case
  "New Fortran"  = "#1f77b4",
  "Old Fortran"  = "#d62728",
  "Deduplicated" = "#2ca02c",
  "R Mimic"      = "#ff7f0e"
)

# Reference model for automated pass/fail (optional, but keep for compatibility)
ref_name <- "Old Fortran"
if (!ref_name %in% names(valid_data_list)) ref_name <- names(valid_data_list)[1]

# Function to generate KS Heatmap for a variable
get_ks_heatmap <- function(var_name, data_list) {
  model_names <- names(data_list)
  n <- length(model_names)
  ks_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(model_names, model_names))
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        ks_matrix[i, j] <- 1.0
      } else {
        # Suppress warnings for ties
        ks_res <- suppressWarnings(ks.test(data_list[[model_names[i]]][[var_name]], 
                                           data_list[[model_names[j]]][[var_name]]))
        ks_matrix[i, j] <- ks_res$p.value
        
        # Simple heuristic for failure if comparing to reference
        if (model_names[j] == ref_name && ks_res$p.value < 0.0001) {
            # Check relative mean diff
            m1 <- mean(data_list[[model_names[i]]][[var_name]])
            m2 <- mean(data_list[[model_names[j]]][[var_name]])
            if (abs(m1-m2)/abs(m2) > 0.05) fail_verification <<- TRUE
        }
      }
    }
  }
  
  # Convert matrix to long format for ggplot
  ks_df <- as.data.frame(as.table(ks_matrix))
  names(ks_df) <- c("Model1", "Model2", "PValue")
  
  # Formatting p-values for display
  ks_df$Label <- ifelse(ks_df$PValue < 0.0001, sprintf("%.1e", ks_df$PValue), sprintf("%.4f", ks_df$PValue))
  
  ggplot(ks_df, aes(x = Model1, y = Model2, fill = PValue)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Label), size = 3, color = ifelse(ks_df$PValue < 0.01, "white", "black")) +
    scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850", 
                         midpoint = 0.05, limits = c(0, 1), 
                         name = "P-Value", 
                         trans = "sqrt") + 
    labs(title = paste("KS Test P-Values:", var_name),
         subtitle = "Red indicates significant distributional difference",
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# 1. Main Plots per Variable
for (p in params_to_check) {
  # Trace Plot
  trace_p <- ggplot(combined_df, aes(x = Iteration_Rel, y = .data[[p]], color = Model)) +
    geom_line(alpha = 0.8) +
    scale_color_manual(values = my_colors) +
    labs(title = paste("Trace Plot:", p), x = "Iteration (Post-Burnin)", y = p) +
    theme_minimal() +
    theme(legend.position = "none") # Hide legend here, show in density or top
  
  if (p == "vara") {
    trace_p <- trace_p + scale_y_log10() + labs(subtitle = "Y-axis: Log10 Scale")
  }

  # Density Ridge Plot
  ridge_p <- ggplot(combined_df, aes(x = .data[[p]], y = Model, fill = Model)) +
    geom_density_ridges(alpha = 0.7, scale = 1.5) +
    scale_fill_manual(values = my_colors) +
    labs(title = paste("Posterior Density:", p), x = p, y = "Model") +
    theme_minimal() +
    theme(legend.position = "none") # Redundant with y-axis
  
  if (p == "vara") {
    ridge_p <- ridge_p + scale_x_log10() + labs(subtitle = "X-axis: Log10 Scale")
  }

  # KS Heatmap
  ks_h <- get_ks_heatmap(p, valid_data_list)
  
  # Layout:
  # Top row: Trace | Ridge Plot
  # Bottom row: KS Heatmap
  
  combined_layout <- (trace_p | ridge_p) / ks_h + 
                     plot_annotation(title = paste("Comparison Summary for Variable:", p),
                                     theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))
  
  print(combined_layout)
}

# 2. Summary Table
# Create a text summary
summary_stats <- do.call(rbind, lapply(names(valid_data_list), function(m) {
  d <- valid_data_list[[m]]
  data.frame(
    Model = m,
    Vare_Mean = mean(d$vare),
    Vare_SD = sd(d$vare),
    Vara_Mean = mean(d$vara),
    Vara_SD = sd(d$vara),
    Mu_Mean = mean(d$mu)
  )
}))

# Normalize format
summary_stats$Vare_Mean <- sprintf("%.4f", summary_stats$Vare_Mean)
summary_stats$Vare_SD   <- sprintf("%.4f", summary_stats$Vare_SD)
summary_stats$Vara_Mean <- sprintf("%.4f", summary_stats$Vara_Mean)
summary_stats$Vara_SD   <- sprintf("%.4f", summary_stats$Vara_SD)
summary_stats$Mu_Mean   <- sprintf("%.4f", summary_stats$Mu_Mean)

# Plot table using grid
grid.newpage()
title_text <- textGrob("Model Comparison Summary", gp = gpar(fontsize = 18, fontface = "bold"), y = 0.95)
grid.draw(title_text)

# Convert DF to Grob table might be complex without extra pkgs, so simplistic text rendering
y_pos <- 0.85
row_height <- 0.05
cols <- names(summary_stats)
x_pos <- seq(0.1, 0.9, length.out = length(cols))

# Headers
for (i in seq_along(cols)) {
  grid.text(cols[i], x = x_pos[i], y = y_pos, gp = gpar(fontface = "bold"))
}

# Rows
for (r in 1:nrow(summary_stats)) {
  y_pos <- y_pos - row_height
  for (i in seq_along(cols)) {
    grid.text(summary_stats[r, i], x = x_pos[i], y = y_pos)
  }
}

# --- Console Output ---
cat("\n--- Comparison Summary Table ---\n")
print(summary_stats)
cat("--------------------------------\n")

dev.off()
cat("Visualization complete.\n")

if (fail_verification) {
  cat("\n!!! VERIFICATION FAILED: Significant differences detected.\n")
  q(status = 1)
} else {
  cat("\nVERIFICATION SUCCESSFUL.\n")
  q(status = 0)
}
