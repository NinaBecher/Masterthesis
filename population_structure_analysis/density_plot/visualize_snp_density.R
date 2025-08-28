############################################################
# SNP Density Plot for Wild Population (Bter1.0)
# Computes SNP density per 10 kb and plots rolling average
# Author: Nina Becher
#############################################################

# --- Load libraries ---
library(ggplot2)
library(dplyr)
library(data.table)
library(zoo)

# --- Set working directory ---
setwd("C:/Users/User/Documents/Masterthesis/Population_genomics")

# --- Read Wild SNP positions ---
wild <- fread("wild_snp_positions.txt", col.names = c("chrom", "pos"))

# --- Bin SNPs into 10 kb windows ---
window <- 10000
wild <- wild %>% mutate(bin = floor(pos / window) * window)

# --- Count SNPs per bin ---
density <- wild %>%
  group_by(chrom, bin) %>%
  summarise(count = n(), .groups = "drop")

# --- Compute cumulative positions for genome-wide plotting ---
chrom_lengths <- density %>%
  group_by(chrom) %>%
  summarise(max_pos = max(bin), .groups = "drop") %>%
  arrange(chrom) %>%
  mutate(cum_start = lag(cumsum(max_pos), default = 0))

density <- density %>%
  left_join(chrom_lengths %>% select(chrom, cum_start), by = "chrom") %>%
  mutate(cum_pos = bin + cum_start)

# --- Compute rolling average (500 bins ~5 Mb) ---
density <- density %>%
  arrange(cum_pos) %>%
  mutate(roll_avg = rollmean(count, k = 500, fill = NA, align = "center"))

# --- Compute chromosome midpoints for x-axis labels ---
chrom_midpoints <- density %>%
  group_by(chrom) %>%
  summarise(mid = (min(cum_pos) + max(cum_pos)) / 2, .groups = "drop")

# --- Plot SNP density ---
# --- Plot SNP density with rotated chromosome labels ---
ggplot(density, aes(x = cum_pos, y = count)) +
  geom_point(color = "#66a61e", alpha = 0.6, size = 1.5) +      
  geom_line(aes(y = roll_avg), color = "#d95f02", size = 1) +  
  geom_hline(yintercept = top1_threshold, color = "red", linetype = "dashed") + 
  scale_x_continuous(
    breaks = chrom_midpoints$mid,
    labels = chrom_midpoints$chrom
  ) +
  theme_bw() +
  labs(
    x = "Chromosome",
    y = "Number of SNPs per 10 kb",
    title = "Reference Genome Bter1.0"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    axis.text.y = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

