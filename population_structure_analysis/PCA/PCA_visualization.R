# ------------------------------------------------------------------------------
# Author: Nina Becher 
# Description: Visualizes PCA results split by BioProject and by Origin
# Reference Genome: Bter_1.0
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Load Required Libraries
# ---------------------------
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)

# ---------------------------
# Set Working Directory
# ---------------------------
setwd("C:/Users/User/Documents/Masterthesis/Results_ORG/PCA")

# ---------------------------
# Load PCA Results
# ---------------------------
load("pca_results_filtered_variants_pca4_excluded.RData")

# ---------------------------
# Set Plot Title and Extract PC Variance
# ---------------------------
JOB_RUN_NAME_FOR_PLOT <- "PCA (reference genome assembly Bter_1.0)"
pc_percent <- round(pca$varprop * 100, 2)
pc1_var <- pc_percent[1]
pc2_var <- pc_percent[2]

# ================================================================
# SECTION 1: PCA Split by BioProject
# ================================================================

# ---------------------------
# Load Sample Metadata (BioProject)
# ---------------------------
pop <- read.table("population_information.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(pop) <- c("sample", "bioproject")

# ---------------------------
# Match PCA Samples to BioProject Info
# ---------------------------
bioproject_info_matched <- pop$bioproject[match(pca$sample.id, pop$sample)]

# ---------------------------
# Create Data Frame for Plotting
# ---------------------------
tab <- data.frame(
  sample.id = pca$sample.id,
  bioproject = factor(bioproject_info_matched),
  EV1 = pca$eigenvect[, 1],
  EV2 = pca$eigenvect[, 2],
  stringsAsFactors = FALSE
)

# ---------------------------
# Quick Check of Sample Assignments
# ---------------------------
print(head(tab))
print(summary(tab$bioproject))

# Print samples per BioProject
samples_by_bioproject <- split(tab$sample.id, tab$bioproject)
for (bp in names(samples_by_bioproject)) {
  cat("\nBioproject:", bp, "\n")
  print(samples_by_bioproject[[bp]])
}

# Identify samples with missing BioProject info
na_samples <- tab$sample.id[is.na(tab$bioproject)]
cat("\nSamples with NA bioproject:\n")
print(na_samples)

# ---------------------------
# Plot: PCA Split by BioProject
# ---------------------------
# ---------------------------
# Plot: PCA Split by BioProject (without jitter)
# ---------------------------
pca_plot <- ggplot(tab, aes(x = EV1, y = EV2, color = bioproject)) +
  geom_point(alpha = 0.85, size = 3) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "PCA\nSplit by BioProject",
    subtitle = "Reference Genome: Bter_1.0",
    x = paste0("Principal Component 1 (", pc1_var, "% Variance Explained)"),
    y = paste0("Principal Component 2 (", pc2_var, "% Variance Explained)"),
    color = "BioProject"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = NA, color = "gray70"),
    text = element_text(family = "Arial")
  )

print(pca_plot)

# ================================================================
# SECTION 2: PCA Split by Sample Origin (Wild vs. Commercial)
# ================================================================

# ---------------------------
# Load Sample Origin Metadata (2 columns: sample ID, origin)
# ---------------------------
sample_origin <- read.table("sample_type_information.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(sample_origin) <- c("sample", "origin")

# ---------------------------
# Reclassify Origins into Two Groups: 'wild' and 'commercial'
# ---------------------------
sample_origin$origin_grouped <- ifelse(tolower(sample_origin$origin) == "wild", "wild", "commercial")

# ---------------------------
# Match Grouped Origins to PCA Sample IDs
# ---------------------------
origin_info_matched <- sample_origin$origin_grouped[match(pca$sample.id, sample_origin$sample)]

# ---------------------------
# Create Data Frame for Plotting
# ---------------------------
tab <- data.frame(
  sample.id = pca$sample.id,
  origin = factor(origin_info_matched),
  EV1 = pca$eigenvect[, 1],
  EV2 = pca$eigenvect[, 2],
  stringsAsFactors = FALSE
)

# ---------------------------
# Check Sample Counts per Group
# ---------------------------
print(summary(tab$origin))

# ---------------------------
# Plot: PCA Split by Origin
# ---------------------------
pca_plot <- ggplot(tab, aes(x = EV1, y = EV2, color = origin)) +
  geom_point(alpha = 0.85, size = 3) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "PCA\nSplit by Origin",
    subtitle = "Reference Genome: Bter_1.0",
    x = paste0("Principal Component 1 (", pc1_var, "% Variance Explained)"),
    y = paste0("Principal Component 2 (", pc2_var, "% Variance Explained)"),
    color = "Origin"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = NA, color = "gray70"),
    text = element_text(family = "Arial")
  )

print(pca_plot)


# ================================================================
# SECTION 3: PCA Split by Sample Origin (Detailed Origin Labels)
# ================================================================

# ---------------------------
# Load Sample Origin Metadata (2 columns: sample ID, origin)
# ---------------------------
sample_origin <- read.table("sample_type_information.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(sample_origin) <- c("sample", "origin")

# ---------------------------
# Match Origin Labels Directly to PCA Sample IDs
# ---------------------------
origin_info_matched <- sample_origin$origin[match(pca$sample.id, sample_origin$sample)]

# ---------------------------
# Format Origin Labels
# ---------------------------
origin_info_matched <- tolower(origin_info_matched)  # normalize case
origin_info_matched[origin_info_matched == "chinese_company"] <- "Fengdeng Biotech"
origin_info_matched <- ifelse(
  origin_info_matched == "Fengdeng Biotech",
  origin_info_matched,
  tools::toTitleCase(origin_info_matched)
)

# ---------------------------
# Create Data Frame for Plotting
# ---------------------------
tab <- data.frame(
  sample.id = pca$sample.id,
  origin = factor(origin_info_matched),
  EV1 = pca$eigenvect[, 1],
  EV2 = pca$eigenvect[, 2],
  stringsAsFactors = FALSE
)

# ---------------------------
# Check Sample Counts per Origin Label
# ---------------------------
print(summary(tab$origin))

# ---------------------------
# Plot: PCA Split by Specific Origin Labels
# ---------------------------
pca_plot <- ggplot(tab, aes(x = EV1, y = EV2, color = origin)) +
  geom_point(alpha = 0.85, size = 3) +
  scale_color_brewer(palette = "Dark2") +  # Or scale_color_manual() if preferred
  labs(
    title = "PCA\nSplit by Supplier",
    subtitle = "Reference Genome: Bter_1.0",
    x = paste0("Principal Component 1 (", pc1_var, "% Variance Explained)"),
    y = paste0("Principal Component 2 (", pc2_var, "% Variance Explained)"),
    color = "Origin"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = NA, color = "gray70"),
    text = element_text(family = "Arial")
  )

print(pca_plot)

