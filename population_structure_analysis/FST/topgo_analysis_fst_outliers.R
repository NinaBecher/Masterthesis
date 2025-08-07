# --------------------------
# Author: Nina
# Purpose: This script performs GO enrichment analysis using the topGO package.
# I use a list of genes of interest and a GO term annotation file from BioMart.
# The goal is to identify enriched Gene Ontology (GO) terms in the categories:
# Biological Process (BP), Cellular Component (CC), and Molecular Function (MF).
# --------------------------

# Load the required library for GO enrichment analysis
library(topGO)

# --- Define File Paths ---
# I define the paths to the GO mapping file and the input gene list
go_mapping_file <- "mart_export.txt"
gene_list_path <- "geneids.txt"

# --- Main Script ---

# Load GO term annotations obtained from BioMart
# I use tab-delimited format and disable string factors for easier handling
go_data <- read.table(go_mapping_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define which columns contain gene IDs and GO term accessions
bombus_col <- "Gene.stable.ID"
go_term_col <- "GO.term.accession"

# Filter out any entries with missing GO terms
# I only want to include genes that are annotated with at least one GO term
filtered_data <- go_data[go_data[[go_term_col]] != "", ]

# Construct the gene2GO list required by topGO
# This maps each gene ID to its list of associated GO terms
gene2GO <- by(
  filtered_data[[go_term_col]],
  filtered_data[[bombus_col]],
  function(x) unique(as.character(x))
)

# Load the list of genes I'm interested in (e.g., differentially expressed genes)
original_gene_list <- read.table(gene_list_path, header = FALSE, stringsAsFactors = FALSE)$V1

# Define the gene universe using all genes in the GO annotation
# I create a binary vector marking which genes are in my list of interest
all_genes_in_annotation <- names(gene2GO)
gene_universe <- factor(as.integer(all_genes_in_annotation %in% original_gene_list), levels = c(0, 1))
names(gene_universe) <- all_genes_in_annotation

# Check how many of my genes of interest are in the gene universe
overlap_count <- sum(gene_universe == 1)

# Stop the script if there's no overlap to avoid running topGO on an empty input
if (overlap_count == 0) {
  stop("No overlap between your gene list and GO mapping. Check IDs or column names.")
}

# --- Loop Through All Three Ontologies ---
# I run the GO enrichment separately for each ontology: BP, CC, MF
ontologies <- c("BP", "CC", "MF")

for (ontology_type in ontologies) {
  
  # Create the topGOdata object for this ontology
  # This object includes the gene universe and GO annotations
  go_data_object <- new("topGOdata",
                        ontology = ontology_type,
                        allGenes = gene_universe,
                        annot = annFUN.gene2GO,
                        gene2GO = gene2GO)
  
  # Run the classic Fisher's exact test and the elim KS test
  fisher_test <- runTest(go_data_object, algorithm = "classic", statistic = "fisher")
  ks_elim_test <- runTest(go_data_object, algorithm = "elim", statistic = "ks")
  
  # --- Save Uncorrected Top 20 Results ---
  # I generate a table of the top 20 GO terms ranked by elim KS p-value
  uncorrected_results <- GenTable(
    go_data_object,
    classicFisher = fisher_test,
    elimKS = ks_elim_test,
    orderBy = "elimKS",
    topNodes = 20
  )
  uncorrected_output_path <- paste0("topGO_", ontology_type, "_uncorrected_top20.tsv")
  write.table(uncorrected_results, file = uncorrected_output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # --- Save BH-Corrected Significant Results ---
  # I apply Benjamini-Hochberg correction to the p-values from the classic Fisher test
  bh_results_table <- GenTable(
    go_data_object,
    classicFisher = fisher_test,
    orderBy = "classicFisher",
    topNodes = 50
  )
  
  bh_results_table$bh.adj.p.value <- p.adjust(bh_results_table$classicFisher, method = "BH")
  
  # I filter for GO terms with an adjusted p-value less than 0.05
  significant_results <- bh_results_table[bh_results_table$bh.adj.p.value < 0.05, ]
  
  # Save the significant results to a file
  significant_output_path <- paste0("topGO_", ontology_type, "_significant_results.tsv")
  write.table(significant_results, file = significant_output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Print the table of significant GO terms
  print(significant_results)
}
