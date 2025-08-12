############################################################
# SNP Density Plot
# Author: Nina
# Description: Calculates SNP density across the genome 
#              and visualizes it with a threshold line.
############################################################

# -------------------
# Load required libraries
# -------------------
library(ggplot2)
library(dplyr)

# -------------------
# Set working directory
# -------------------
setwd("C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/density_plot")

# -------------------
# Load SNP position data
# -------------------
# Expected format: two columns [CHROM, POS]
snp_data <- read.table(
  "snp_positions.txt",
  header = FALSE,
  col.names = c("CHROM", "POS")
)

# -------------------
# Parameters
# -------------------
window_size <- 10000  # Window size in base pairs (e.g., 10 kb)

# -------------------
# Calculate SNP density per window for each chromosome
# -------------------
snp_density <- snp_data %>%
  group_by(CHROM) %>%
  mutate(window = floor(POS / window_size)) %>%
  group_by(CHROM, window) %>%
  summarise(SNPs = n(), .groups = "drop") %>%
  mutate(
    start = window * window_size,
    end = start + window_size
  )

# -------------------
# Create cumulative positions for continuous x-axis
# -------------------
chr_lengths <- snp_density %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(end))

chr_offsets <- c(0, cumsum(head(chr_lengths$chr_len, -1)))
names(chr_offsets) <- chr_lengths$CHROM

snp_density <- snp_density %>%
  mutate(pos_cum = start + chr_offsets[CHROM])

# -------------------
# Plot SNP density
# -------------------
ggplot(snp_density, aes(x = pos_cum, y = SNPs)) +
  geom_point(size = 0.8) +
  theme_bw() +
  labs(
    title = "SNP Density Across the Genome",
    subtitle = "Reference Genome: Bter1.0",
    x = "Genomic position (bp)",
    y = paste("Number of SNPs (", window_size / 1000, " kb windows)", sep = "")
  ) +
  geom_hline(
    yintercept = mean(snp_density$SNPs) + 3 * sd(snp_density$SNPs),
    color = "red",
    linetype = "dashed"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),   # Bold & centered title
    plot.subtitle = element_text(face = "italic", hjust = 0.5) # Italic & centered subtitle
  )
