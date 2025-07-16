#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr) # For read_delim

# --- Configuration ---
# Path to the FST output file from Pixy
# This file is assumed to contain FST values for 'wild' vs 'commercial' populations
FST_FILE <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/results/pixy_output_commercial_strict_by_supplier/pixy_fst.txt"

# Output directory for results (tables and plots)
OUTPUT_DIR <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/jobscripts/pixy_output_origin_org/fst_zscore_outlier_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Population names for the specific comparison you're interested in
# These should match the 'pop1' and 'pop2' names in your FST file
POP1_NAME <- "wild"
POP2_NAME <- "commercial"

# Significance threshold for FDR corrected p-value
FDR_THRESHOLD <- 0.05

# --- Analysis ---

cat(paste0("Starting FST Z-score and p-value outlier analysis for ", POP1_NAME, " vs ", POP2_NAME, "...\n"))
cat(paste0("Reading FST data from: ", FST_FILE, "\n"))

# 1. Load the FST data
# Using read_delim for potentially faster reading of large files
fst_data <- read_delim(FST_FILE, delim = "\t", show_col_types = FALSE)

# 2. Filter for the specific population comparison
cat(paste0("Filtering data for FST between '", POP1_NAME, "' and '", POP2_NAME, "'...\n"))
fst_target_comparison <- fst_data %>%
  filter((pop1 == POP1_NAME & pop2 == POP2_NAME) | (pop1 == POP2_NAME & pop2 == POP1_NAME)) %>%
  distinct(chromosome, window_pos_1, window_pos_2, .keep_all = TRUE) # Remove potential duplicates if present

# Check if data exists for the target comparison
if (nrow(fst_target_comparison) == 0) {
  stop(paste0("No FST data found for comparison '", POP1_NAME, "' vs '", POP2_NAME, "'. Please check your FST_FILE and POP_NAMEs."))
}

# 3. Handle NA values in avg_wc_fst (important for Z-score calculation)
fst_target_comparison <- fst_target_comparison %>%
  filter(!is.na(avg_wc_fst))

if (nrow(fst_target_comparison) == 0) {
  stop(paste0("No valid FST values remaining after NA filtering for '", POP1_NAME, "' vs '", POP2_NAME, "'. Exiting."))
}

# 4. Calculate Z-scores
cat("Calculating Z-scores...\n")
mean_fst <- mean(fst_target_comparison$avg_wc_fst)
sd_fst <- sd(fst_target_comparison$avg_wc_fst)

# Handle cases where standard deviation might be zero (i.e., all FST values are identical)
if (sd_fst == 0) {
  warning("Standard deviation of FST values is zero. All Z-scores will be 0, and p-values will be 1.")
  fst_target_comparison$z_score <- 0
  fst_target_comparison$p_value <- 1 # All p-values become 1 if no variance
  fst_target_comparison$p_value_fdr <- 1
  fst_target_comparison$neglog10_p_value_fdr <- 0
} else {
  fst_target_comparison$z_score <- (fst_target_comparison$avg_wc_fst - mean_fst) / sd_fst

  # 5. Calculate p-values (one-tailed: interested in higher FST values)
  cat("Calculating p-values and applying FDR correction...\n")
  fst_target_comparison$p_value <- pnorm(fst_target_comparison$z_score, lower.tail = FALSE)

  # 6. Apply Benjamini-Hochberg (FDR) multiple testing correction
  fst_target_comparison$p_value_fdr <- p.adjust(fst_target_comparison$p_value, method = "fdr")

  # Create -log10(FDR p-value) for Manhattan plot
  fst_target_comparison$neglog10_p_value_fdr <- -log10(fst_target_comparison$p_value_fdr)
}

# 7. Identify outlier windows based on the FDR threshold
outlier_windows <- fst_target_comparison %>%
  filter(p_value_fdr < FDR_THRESHOLD)

cat(paste0("Found ", nrow(outlier_windows), " outlier windows at FDR < ", FDR_THRESHOLD, ".\n"))

# 8. Save results
output_results_path <- file.path(OUTPUT_DIR, paste0("fst_zscore_results_", POP1_NAME, "_vs_", POP2_NAME, ".csv"))
output_outliers_path <- file.path(OUTPUT_DIR, paste0("fst_zscore_outlier_windows_", POP1_NAME, "_vs_", POP2_NAME, ".csv"))

write_csv(fst_target_comparison, output_results_path)
write_csv(outlier_windows, output_outliers_path)

cat(paste0("Detailed results saved to: ", output_results_path, "\n"))
cat(paste0("Outlier windows saved to: ", output_outliers_path, "\n"))

# 9. Generate Manhattan Plot
cat("Generating Manhattan plot...\n")

# Prepare data for plotting (calculate cumulative positions for X-axis)
plot_data <- fst_target_comparison %>%
  mutate(chromosome = as.factor(chromosome),
         window_mid_pos = (window_pos_1 + window_pos_2) / 2)

# Calculate cumulative positions for chromosomes to arrange them on X-axis
data_cum <- plot_data %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(window_mid_pos, na.rm = TRUE)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(chromosome, bp_add)

plot_data_manhattan <- plot_data %>%
  left_join(data_cum, by = "chromosome") %>%
  mutate(bp_cum = window_mid_pos + bp_add)

# Prepare axis labels for the chromosomes
axis_set <- plot_data_manhattan %>%
  group_by(chromosome) %>%
  summarize(center = mean(bp_cum))

manhattan_plot <- ggplot(plot_data_manhattan, aes(x = bp_cum, y = neglog10_p_value_fdr, color = as.factor(chromosome))) +
  geom_point(alpha = 0.8, size = 1.3) +
  # Use a repeating color palette for chromosomes
  scale_color_manual(values = rep(c("grey", "skyblue"), length(unique(plot_data_manhattan$chromosome)))) +
  scale_x_continuous(label = axis_set$chromosome, breaks = axis_set$center) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", color = "red") + # Draw FDR threshold line
  labs(x = "Chromosome", y = expression(-log[10](FDR~P-value)), # LaTeX for P-value subscript
       title = paste0("Manhattan Plot of FST Outliers (", POP1_NAME, " vs ", POP2_NAME, ")")) +
  theme_minimal() +
  theme(
    legend.position = "none", # Hide legend for chromosome colors
    panel.grid.major.x = element_blank(), # Remove major vertical grid lines
    panel.grid.minor.x = element_blank(), # Remove minor vertical grid lines
    axis.text.x = element_text(angle = 45, hjust = 1) # Angle chromosome labels
  )

# Save the Manhattan plot
manhattan_plot_path <- file.path(OUTPUT_DIR, paste0("manhattan_plot_", POP1_NAME, "_vs_", POP2_NAME, ".png"))
ggsave(manhattan_plot_path, manhattan_plot, width = 12, height = 6, dpi = 300)
cat(paste0("Manhattan plot saved to: ", manhattan_plot_path, "\n"))

cat("\nFST Z-score and p-value outlier analysis complete.\n")
