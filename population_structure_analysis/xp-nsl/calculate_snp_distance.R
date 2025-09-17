#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Process data from a large file stream to improve memory efficiency.
# The `cut` command selects columns 3 onwards from the input file,
# and `read.table` reads this piped data directly into a data frame.
file_path <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/Population_genomics_results/xp-nsl/Bter1_0/results/combined_all_xpnsl.txt"
file_stream <- pipe(paste("cut -f3-", file_path))

# Read the data, assigning meaningful column names for clarity.
# The data corresponds to an extended haplotype homozygosity (EHH) analysis,
# specifically using the nSL (normalized S-L) statistic.
nsl_data <- read.table(
  file = file_stream,
  header = FALSE,
  col.names = c("chrom", "position", "freq1", "nsl1", "nsl0",
                "unstandard_nsl", "standard_nsl", "unusual_flag")
)

# Define the significance cutoff for identifying selection signatures.
# This value represents the top 1% of the standardized nSL distribution.
# This specific cutoff was pre-calculated for the Bter1_0 dataset.
top_cutoff <- 1.329084

# Isolate the top 1% of SNPs based on the absolute value of the standardized nSL score.
# These SNPs are candidates for recent positive selection.
top_one_percent <- subset(nsl_data, abs(standard_nsl) >= top_cutoff)

# Initialize a vector to store the calculated distances.
distance_1 <- vector()

# Iterate through each of the top 1% SNPs to find the nearest
# SNP with a standardized nSL value below the top 1% cutoff.
for (i in 1:nrow(top_one_percent)) {
  current_snp <- top_one_percent[i, ]
  row_index <- as.numeric(row.names(current_snp))

  # Define a search window of ±100 SNPs around the current SNP of interest.
  downstream_window_end <- max(row_index - 100, 1)
  upstream_window_end <- min(row_index + 100, nrow(nsl_data))

  # Search for the nearest SNP with a standard_nsl value below the cutoff
  # within the downstream window.
  downstream_subset <- nsl_data[downstream_window_end:row_index, ]
  downstream_subset$diff <- current_snp$position - downstream_subset$position
  downstream_min_dist <- min(
    subset(
      downstream_subset,
      chrom == current_snp$chrom & diff != 0 & abs(standard_nsl) <= 1
    )$diff,
    na.rm = TRUE
  )

  # Search for the nearest SNP with a standard_nsl value below the cutoff
  # within the upstream window.
  upstream_subset <- nsl_data[row_index:upstream_window_end, ]
  upstream_subset$diff <- upstream_subset$position - current_snp$position
  upstream_min_dist <- min(
    subset(
      upstream_subset,
      chrom == current_snp$chrom & diff != 0 & abs(standard_nsl) <= 1
    )$diff,
    na.rm = TRUE
  )

  # Store the smaller of the two minimum distances (upstream and downstream).
  distance_1 <- c(distance_1, min(c(downstream_min_dist, upstream_min_dist), na.rm = TRUE))
}

# Clean the result vector by removing any infinite values,
# which may result from a lack of nearby SNPs below the cutoff.
distance_1[is.infinite(distance_1)] <- NA
distance_1 <- na.omit(distance_1)

# Print a summary of the distances to understand the distribution.
summary_stats <- summary(distance_1)
print(summary_stats)

# Save the resulting distances for future analysis.
# The .Rdata file format is efficient for storing R objects.
output_file <- "/lustre/miifs01/project/m2_jgu-salmosex/nina/Population_genomics_results/xp-nsl/Bter1_0/results/top1pct_distances_Bter1_0.Rdata"
save(distance_1, file = output_file)
