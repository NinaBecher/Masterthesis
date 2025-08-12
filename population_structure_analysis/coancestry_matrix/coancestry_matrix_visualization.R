#Code to visualize coancestry matrix 


library(pheatmap)

# Set the working directory
setwd("C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/ancestorymatrix")

# Read and prepare sample information with a flexible separator
sample_info <- read.table("sample_type_information.txt", header = FALSE, sep = "", na.strings = "", stringsAsFactors = FALSE)
colnames(sample_info) <- c("Sample_ID", "Supplier", "Country")
rownames(sample_info) <- sample_info$Sample_ID

# ----------------------------------------------------------------------
# ➡️ New Section: Re-classify suppliers as "Wild" or "Commercial"
# ----------------------------------------------------------------------

# Create a new column to store the simplified classification
sample_info$Supplier_Class <- sample_info$Supplier

# Re-classify all non-"wild" suppliers as "commercial"
sample_info$Supplier_Class[!grepl("wild", sample_info$Supplier_Class, ignore.case = TRUE)] <- "commercial"

# ----------------------------------------------------------------------
# ➡️ Read and subset the genotype matrix (same as before)
# ----------------------------------------------------------------------

# Read the genotype matrix
geno <- read.table("genotype_matrix_normalized.txt", header = TRUE, sep = "\t", na.strings = "NA")

# Define the samples you want to EXCLUDE.
samples_to_exclude <- c("SRR11742810", "SRR11744704", "SRR11780512", "SRR11780662", 
                        "SRR11880654", "ERR6055009", "ERR6055007", "ERR6055006", 
                        "ERR6055008", "ERR6558189", "SRR3928722", "SRR3928782", 
                        "SRR3928789", "SRR3928794")

# Remove CHROM and POS columns and convert to a numeric matrix
geno_mat <- as.matrix(sapply(geno[, -c(1, 2)], as.numeric))

# Subset the genotype matrix
geno_mat_subset <- geno_mat[, !colnames(geno_mat) %in% samples_to_exclude]

# ----------------------------------------------------------------------
# ➡️ Prepare annotations using only the simplified supplier class
# ----------------------------------------------------------------------

# Filter the sample info to include only the samples in the subsetted matrix
sample_info_filtered <- sample_info[colnames(geno_mat_subset), ]
annotation_df <- data.frame(
  Supplier_Type = sample_info_filtered$Supplier_Class
)
rownames(annotation_df) <- rownames(sample_info_filtered)

# ----------------------------------------------------------------------
# ➡️ Re-calculate co-ancestry and plot with new annotations
# ----------------------------------------------------------------------

# Calculate co-ancestry matrix
co_ancestry_subset <- cor(geno_mat_subset, use = "pairwise.complete.obs", method = "pearson")

# Visualize the new co-ancestry matrix as a heatmap with annotations
pheatmap(co_ancestry_subset,
         clustering_method = "ward.D2",
         main = "Co-Ancestry Matrix Heatmap with Supplier Type\n(Bter1.0 Reference Genome)",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_row = annotation_df,
         annotation_col = annotation_df)

