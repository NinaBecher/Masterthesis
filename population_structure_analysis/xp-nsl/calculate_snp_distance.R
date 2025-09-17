#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(ggpubr)

# Use a pipe to pre-process the file and select only the relevant columns.
# We skip the first column (the filename) and the second column (the dot).
# The input file is separated by multiple spaces, so we use `tr -s ' '` to
# squeeze multiple spaces into a single one before piping to `cut`.
file_stream <- pipe("cut -f2- /lustre/miifs01/project/m2_jgu-salmosex/nina/Population_genomics_results/xp-nsl/BomTerr1_2/results/combined_all_nsl_BomTerr1_2.txt")

# Read in data from the pipe with correct column names (8 columns from selscan)
nsl_data <- read.table(file = file_stream, header = FALSE,
                       col.names = c("chrom", "position", "freq1", "nsl1", "nsl0",
                                     "unstandard_nsl", "standard_nsl", "unusual_flag"))

# Top 1% cutoff (previously calculated)
top_cutoff <- 1.329084

# Identify top 1% SNPs using standardized nSL
top_one_percent <- subset(nsl_data, abs(standard_nsl) >= top_cutoff)

# Initialize vector to store distances
distance_1 <- vector()

# Loop over each top SNP
for (row in 1:nrow(top_one_percent)) {
  chrom_position <- top_one_percent[row, ]$chrom
  snp_position <- top_one_percent[row, ]$position
  snp_of_interest <- as.numeric(row.names(top_one_percent[row, ]))
  
  # Upstream and downstream window (±100 SNPs)
  downstream_snps <- max(snp_of_interest - 100, 1)
  upstream_snps <- min(snp_of_interest + 100, nrow(nsl_data))
  
  # Downstream SNPs
  downstream_snps_tmp <- nsl_data[downstream_snps:snp_of_interest, ]
  downstream_snps_tmp$diff <- snp_position - downstream_snps_tmp$position
  downstream_snps_min <- abs(min(subset(downstream_snps_tmp,
                                        chrom == chrom_position &
                                          diff != 0 &
                                          abs(standard_nsl) <= 1)$diff, na.rm = TRUE))
  
  # Upstream SNPs
  upstream_snps_tmp <- nsl_data[snp_of_interest:upstream_snps, ]
  upstream_snps_tmp$diff <- upstream_snps_tmp$position - snp_position
  upstream_snps_min <- abs(min(subset(upstream_snps_tmp,
                                      chrom == chrom_position &
                                        diff != 0 &
                                        abs(standard_nsl) <= 1)$diff, na.rm = TRUE))
  
  # Store minimum distance
  distance_1 <- c(distance_1, min(c(downstream_snps_min, upstream_snps_min), na.rm = TRUE))
}

# Remove any remaining infinite values
distance_1[is.infinite(distance_1)] <- NA
distance_1 <- na.omit(distance_1)

# Summary statistics of distances
summary_stats <- summary(distance_1)
print(summary_stats)

# Save distances to file
output_file <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/Population_genomics_results/xp-nsl/BomTerr1_2/results/top1pct_distances_BomTerr1_2.Rdata"
save(distance_1, file = output_file)
