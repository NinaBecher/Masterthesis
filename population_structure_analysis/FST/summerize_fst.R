#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Summarize FST values from pixy_fst.txt and save as CSV
# ------------------------------------------------------------------------------

# Load required packages
library(dplyr)
library(readr)

# Define input and output paths (absolute paths on Mogon)
input_file <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/results/pixy_output_commercial_strict_by_supplier/pixy_fst.txt"
output_file <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/results/fst_summary_output.csv"

# Read the FST data
fst <- read_delim(input_file, delim = "\t", show_col_types = FALSE)

# Remove rows with NA fst
fst_clean <- fst %>% filter(!is.na(avg_wc_fst))

# Create summary stats for each population pair
fst_summary <- fst_clean %>%
  mutate(pop_pair = paste(pop1, pop2, sep = "_")) %>%
  group_by(pop_pair) %>%
  summarise(
    mean_fst = mean(avg_wc_fst, na.rm = TRUE),
    median_fst = median(avg_wc_fst, na.rm = TRUE),
    sd_fst = sd(avg_wc_fst, na.rm = TRUE),
    n_windows = n(),
    .groups = "drop"
  )

# Write the summary to CSV
write.csv(fst_summary, output_file, row.names = FALSE)

# Confirmation message (shows in SLURM .out log)
cat("✅ Summary FST CSV written to:\n", output_file, "\n")

