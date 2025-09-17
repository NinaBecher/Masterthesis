#!/usr/bin/env Rscript
##############################################################################
## Purpose: Convert VCF -> haplotype/genotype matrix + SNP map for selscan nSL
## Compatible with haploid + diploid samples. This version is corrected
## to handle all possible genotype values and missing data.
##############################################################################

# Load required libraries
libraries <- c("gdsfmt", "SNPRelate", "MASS")
for (lib in libraries) {
  if (!require(lib, character.only = TRUE)) {
    stop(paste("Please install the R package:", lib))
  }
}

# Read command-line arguments
args <- commandArgs(TRUE)

# Check if the correct number of arguments were provided
if (length(args) < 4) {
  stop("Usage: Rscript vcf_to_genomatrix_hap.R <input_vcf> <input_gds> <output_hap> <output_map>", call.=FALSE)
}

input_vcf   <- args[1]  # Input VCF (can be gzipped)
input_gds   <- args[2]  # Temporary GDS file
output_hap  <- args[3]  # Haplotype/genotype matrix output
output_map  <- args[4]  # SNP map output

message(paste("Processing VCF:", input_vcf))
message(paste("Output haplotype matrix:", output_hap))
message(paste("Output SNP map:", output_map))

# Function to convert VCF -> hap + map
vcf_to_genomatrix_hap <- function(vcf_file, gds_file, hap_file, map_file) {
  # Convert VCF -> GDS, keeping only biallelic SNPs.
  # This is the first filtering step.
  message("Step 1: Converting VCF to GDS format...")
  snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only", verbose=FALSE)
  
  # Open GDS
  message("Step 2: Opening GDS file and reading genotype data...")
  genofile <- snpgdsOpen(gds_file)
  
  # Read the full genotype matrix from the GDS file.
  # The values are 0 (ref), 1 (het), 2 (alt), or NA (missing).
  genotype_matrix <- read.gdsn(index.gdsn(genofile, "genotype"))
  
  # IMPORTANT: Correctly convert the genotype matrix to a haploid matrix (0s and 1s).
  # We will also explicitly handle missing data (NA).
  message("Step 3: Converting genotype matrix to haploid matrix (0s and 1s)...")
  
  # Create a new matrix of the same dimensions, initialized with a placeholder
  # for missing data. Selscan uses '9' by default for missing alleles.
  haplotype_matrix <- matrix(9, nrow=nrow(genotype_matrix), ncol=ncol(genotype_matrix))
  
  # Set reference homozygotes (0) to 0.
  haplotype_matrix[genotype_matrix == 0] <- 0
  
  # Set heterozygotes (1) and alternate homozygotes (2) to 1.
  # This creates a presence/absence matrix of the alternate allele, which
  # is what selscan expects.
  haplotype_matrix[genotype_matrix == 1] <- 1
  haplotype_matrix[genotype_matrix == 2] <- 1
  
  # Read SNP positions and chromosomes.
  snp_chrom <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
  snp_pos   <- read.gdsn(index.gdsn(genofile, "snp.position"))
  
  # Build the map dataframe for selscan.
  output_df <- data.frame(
    chr = snp_chrom,
    snp_name = ".",
    pos1 = snp_pos,
    pos2 = snp_pos
  )
  
  # Write the outputs. Using write.table with " " as separator for consistency.
  message("Step 4: Writing output files...")
  write.table(haplotype_matrix, file = hap_file, sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(output_df, file = map_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  message("Done.")
  
  # Close the GDS file to release the resource.
  snpgdsClose(genofile)
}

# Run the function with provided arguments
vcf_to_genomatrix_hap(input_vcf, input_gds, output_hap, output_map)

