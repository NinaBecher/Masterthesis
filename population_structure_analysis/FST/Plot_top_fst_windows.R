#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)

# Input
fst_file <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/results/pixy_output_commercial_strict_by_supplier/pixy_fst.txt"

# Read and prep
fst <- read_delim(fst_file, delim = "\t", show_col_types = FALSE) %>%
  filter(!is.na(avg_wc_fst)) %>%
  mutate(pop_pair = paste(pop1, pop2, sep = "_"))

# Calculate 99th percentile threshold
thresholds <- fst %>%
  group_by(pop_pair) %>%
  summarise(threshold = quantile(avg_wc_fst, 0.99, na.rm = TRUE))

# Join threshold and flag top 1%
fst <- fst %>%
  left_join(thresholds, by = "pop_pair") %>%
  mutate(top_1pct = avg_wc_fst >= threshold)

# Plot with outliers highlighted
plot <- ggplot(fst, aes(x = window_pos_1, y = avg_wc_fst, color = top_1pct)) +
  geom_point(alpha = 0.5, size = 0.8) +
  facet_wrap(~chromosome, scales = "free_x") +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "FST Genome Scan (Top 1% Highlighted)",
       x = "Genomic Window Start",
       y = "Average WC FST",
       color = "Top 1%")

# Save PDF
ggsave("/lustre/miifs01/project/m2_jgu-salmosex/nina/results/fst_outlier_plot.pdf",
       plot, width = 12, height = 8)

cat("✅ Plot saved to: fst_outlier_plot.pdf\n")

