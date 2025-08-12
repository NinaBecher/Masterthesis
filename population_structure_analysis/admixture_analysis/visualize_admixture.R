################################################################################
# Admixture Analysis and Plotting Script
# Author: Nina
# Description: 
#   Loads admixture Q files and metadata, processes ancestry proportions,
#   visualizes admixture results with ggplot2, assigns clusters, and exports results.
################################################################################

# ----------------------------
# Load required libraries
# ----------------------------
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)
library(tidyr)

# ----------------------------
# Define file paths (update as needed)
# ----------------------------
q_file <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/all_bioprojects_snps_filtered_Bter1_0_for_admixture.5.q"
fam_file <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/famfile.txt"
pop_file <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/sample_type_information.txt"
cluster_assignments_file <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/cluster_assignments.csv"
output_pdf <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/admixture_plot_K5_combined_demarcated.pdf"
output_assignments_csv <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/cluster_assignments.csv"

# ----------------------------
# Load and prepare admixture data
# ----------------------------

# Read Q file with ancestry proportions
q_data <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
num_k <- ncol(q_data)  # Number of ancestry components (K)
colnames(q_data) <- paste0("V", 1:num_k)

# Read FAM file for individual IDs
fam_data <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")

# Add IID column to Q data
q_data$IID <- fam_data$IID

# Read population info
pop_data <- read.table(pop_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pop_data) <- c("IID", "Population")

# Standardize population names ("Wild" vs "Commercial")
pop_data$Population <- ifelse(toupper(trimws(pop_data$Population)) == "WILD", "Wild", "Commercial")

# Merge Q data with population info by IID
admix_data <- merge(q_data, pop_data, by = "IID")

# Reorder columns and sort individuals by population and first ancestry component
admix_data <- admix_data[, c("IID", paste0("V", 1:num_k), "Population")]
admix_data <- admix_data %>% arrange(Population, desc(V1))

# ----------------------------
# Load and merge cluster assignments
# ----------------------------

# Read raw cluster assignments CSV (semicolon-separated values in one column)
cluster_assignments_raw <- read.csv(cluster_assignments_file, 
                                    header = FALSE, stringsAsFactors = FALSE)

# Split the first column by semicolon into separate columns
cluster_assignments <- separate(cluster_assignments_raw, 
                                col = V1, 
                                into = c("BioProject", "IID", "V1", "V2", "V3", "V4", "V5", "V6", "Population", "Assigned_Cluster"), 
                                sep = ";", convert = TRUE)

# Remove header row if present
if (cluster_assignments$IID[1] == "IID") {
  cluster_assignments <- cluster_assignments[-1, ]
}

# Merge BioProject info into admixture data
admix_data <- merge(admix_data, cluster_assignments[, c("IID", "BioProject")], by = "IID")

# ----------------------------
# Reshape data for plotting by BioProject
# ----------------------------

admix_data_long <- admix_data %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Ancestry_Component",
    values_to = "Proportion"
  )

# Make BioProject a factor to preserve order on x-axis
admix_data_long$BioProject <- factor(admix_data_long$BioProject, levels = unique(admix_data_long$BioProject))

# Sort by BioProject if needed
admix_data_long <- admix_data_long[order(admix_data_long$BioProject), ]

# ----------------------------
# Plot admixture proportions per BioProject
# ----------------------------

colors_K <- viridis(n = num_k, option = "D")

p_bioproject <- ggplot(admix_data_long, aes(x = BioProject, y = Proportion, fill = Ancestry_Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = colors_K) +
  theme_minimal() +
  labs(title = paste0("Admixture Model (K=", num_k, ")"),
       x = "BioProject",
       y = "Ancestry Proportion",
       fill = "Ancestry Component") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_bioproject)

# ----------------------------
# Prepare data for individual-level plotting with population demarcation
# ----------------------------

# Reshape data again for plotting by individual and population
admix_data_long_ind <- melt(admix_data,
                           id.vars = c("IID", "Population"),
                           variable.name = "Ancestry_Component",
                           value.name = "Proportion")

# Set IID factor levels to preserve order
admix_data_long_ind$IID <- factor(admix_data_long_ind$IID, levels = admix_data$IID)

# Calculate population segments for background shading and labels
ordered_iids <- levels(admix_data_long_ind$IID)
population_segments <- admix_data %>%
  group_by(Population) %>%
  summarize(
    xmin = min(as.numeric(factor(IID, levels = ordered_iids))) - 0.5,
    xmax = max(as.numeric(factor(IID, levels = ordered_iids))) + 0.5,
    label_x = (min(as.numeric(factor(IID, levels = ordered_iids))) +
                 max(as.numeric(factor(IID, levels = ordered_iids)))) / 2
  )

# ----------------------------
# Plot admixture proportions per individual with population demarcation
# ----------------------------

p_individual <- ggplot(admix_data_long_ind, aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
  geom_rect(data = population_segments,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = rep(c("grey90", "grey80"), length.out = nrow(population_segments)),
            alpha = 0.3) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = colors_K) +
  geom_text(data = population_segments,
            aes(x = label_x, y = 1.05, label = paste0(Population, " Bees")),
            inherit.aes = FALSE, vjust = 0, fontface = "bold", size = 4.5) +
  theme_minimal() +
  labs(
    title = paste0("Admixture Model (K=", num_k, ")"),
    subtitle = "Genetic Ancestry Proportions",
    x = NULL,
    y = "Ancestry Proportion",
    fill = "Ancestry Component"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, margin = margin(t = 5)),
    axis.ticks.x = element_line(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = unit(c(1, 1, 3, 0.5), "cm")
  ) +
  scale_x_discrete(expand = expansion(add = 1)) +  # extra space on sides
  coord_cartesian(ylim = c(0, 1.07), clip = "off")

print(p_individual)

# ----------------------------
# Save plots to PDF
# ----------------------------

ggsave(output_pdf, plot = p_bioproject, width = 24, height = 6)
ggsave(output_pdf, plot = p_individual, width = 24, height = 6)

# ----------------------------
# Assign individuals to clusters based on highest ancestry proportion
# ----------------------------

admix_assignments <- admix_data %>%
  mutate(
    Assigned_Cluster = apply(select(., starts_with("V")), 1, function(x) {
      paste0("V", which.max(x))
    })
  ) %>%
  arrange(Assigned_Cluster)

# Save assignments to CSV
write.csv(admix_assignments, output_assignments_csv, row.names = FALSE)

# View first few assignments
head(admix_assignments)

