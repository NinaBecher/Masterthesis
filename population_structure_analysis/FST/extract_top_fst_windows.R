#!/usr/bin/env Rscript

library(dplyr)
library(readr)

# Paths (adjust if needed)
fst_file <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/results/pixy_output_commercial_strict_by_supplier/pixy_fst.txt"
output_file <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/results/top1pct_fst_windows.csv"

# Read FST data
fst <- read_delim(fst_file, delim = "\t", show_col_types = FALSE)

# Clean and compute threshold
fst_clean <- fst %>%
  filter(!is.na(avg_wc_fst)) %>%
  mutate(pop_pair = paste(pop1, pop2, sep = "_"))

# Top 1% threshold per population pair
top_fst <- fst_clean %>%
  group_by(pop_pair) %>%
  mutate(threshold = quantile(avg_wc_fst, 0.99, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(avg_wc_fst >= threshold)

# Save top 1% FST windows
write.csv(top_fst, output_file, row.names = FALSE)

cat("✅ Top 1% FST windows written to:\n", output_file, "\n")

