# ----------------------------
# Admixture Analysis and Plotting Script (K=3)
# Author: Adapted for messy population files
# ----------------------------

# ----------------------------
# Load required libraries
# ----------------------------
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)
library(tidyr)
library(scales) # for hue_pal()

# ----------------------------
# Set working directory and file paths
# ----------------------------
setwd("C:/Users/User/Documents/Masterthesis/Population_genomics/Bter1_0/admixture")

q_file <- "all_bioprojects_snps_strict_filtered_Bter1_0_for_admixture.3.q"
fam_file <- "all_bioprojects_snps_strict_filtered_Bter1_0_for_admixture.fam"
pop_file <- "sample_type_information.txt"
bioproject_file <- "population_information.txt" # File containing BioProject info
output_bioproject_pdf <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/admixture_plot_K3_by_bioproject.pdf"
output_individual_pdf <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/admixture_plot_K3_combined_demarcated.pdf"
output_assignments_csv <- "C:/Users/User/Documents/Masterthesis/results_bomterr/results_org/admixture/data/cluster_assignments.csv"

# ----------------------------
# Load and prepare admixture data
# ----------------------------
q_data <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
num_k <- ncol(q_data) # Number of ancestry components (K)
colnames(q_data) <- paste0("V", 1:num_k)

fam_data <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")

q_data$IID <- fam_data$IID

# Read population info
pop_data <- read.table(pop_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
pop_data <- pop_data[, 1:2] # take only IID and Population
colnames(pop_data) <- c("IID", "Population")

# Standardize population names ("Wild" vs "Commercial")
pop_data$Population <- ifelse(toupper(trimws(pop_data$Population)) == "WILD", "Wild", "Commercial")

# Read BioProject info from the new file
bioproject_data <- read.table(bioproject_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bioproject_data) <- c("IID", "BioProject")
bioproject_data$IID <- trimws(bioproject_data$IID)
bioproject_data$BioProject <- trimws(bioproject_data$BioProject)

# Merge Q data with population and bioproject info
admix_data <- merge(q_data, pop_data, by = "IID")
admix_data <- merge(admix_data, bioproject_data, by = "IID")

# Keep only K columns and the new info
admix_data <- admix_data[, c("IID", paste0("V", 1:num_k), "Population", "BioProject")]

# Sort individuals by Population and V1 proportion
admix_data <- admix_data %>% arrange(Population, desc(V1))

# ----------------------------
# Prepare long-format data for plotting
# ----------------------------
admix_data_long <- admix_data %>%
  pivot_longer(cols = starts_with("V"), names_to = "Ancestry_Component", values_to = "Proportion")

admix_data_long$BioProject <- factor(admix_data_long$BioProject, levels = unique(admix_data_long$BioProject))
admix_data_long <- admix_data_long[order(admix_data_long$BioProject), ]

# ----------------------------
# Colors
# ----------------------------
colors_K <- viridis(n = num_k, option = "D")
bioproject_colors <- scales::hue_pal()(length(unique(admix_data$BioProject)))
names(bioproject_colors) <- unique(admix_data$BioProject)

# ----------------------------
# Plot admixture proportions per BioProject
# ----------------------------
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
# Individual-level plot with population demarcation
# ----------------------------
admix_data_long_ind <- reshape2::melt(
  admix_data,
  id.vars = c("IID", "Population", "BioProject"),
  measure.vars = paste0("V", 1:num_k),
  variable.name = "Ancestry_Component",
  value.name = "Proportion"
)
admix_data_long_ind$IID <- factor(admix_data_long_ind$IID, levels = admix_data$IID)
admix_data_long_ind$Ancestry_Component <- factor(admix_data_long_ind$Ancestry_Component, levels = paste0("V", 1:num_k))
admix_data_long_ind <- droplevels(admix_data_long_ind)

# Population segments for top labels and background shading
ordered_iids <- levels(admix_data_long_ind$IID)
population_segments <- admix_data %>%
  group_by(Population) %>%
  summarize(
    xmin = min(as.numeric(factor(IID, levels = ordered_iids))) - 0.5,
    xmax = max(as.numeric(factor(IID, levels = ordered_iids))) + 0.5,
    label_x = (min(as.numeric(factor(IID, levels = ordered_iids))) +
                 max(as.numeric(factor(IID, levels = ordered_iids)))) / 2,
    .groups = "drop"
  )

# BioProject labels and background for the x-axis
bioproject_labels <- admix_data %>%
  group_by(BioProject) %>%
  summarize(
    start = min(as.numeric(factor(IID, levels = ordered_iids))),
    end = max(as.numeric(factor(IID, levels = ordered_iids))),
    midpoint = (start + end) / 2,
    .groups = "drop"
  )

# Build the plot incrementally to avoid the '+.gg' error
p_individual <- ggplot(admix_data_long_ind, aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
  
  # Add colored background for BioProjects
  geom_rect(data = bioproject_labels,
            aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = bioproject_colors[bioproject_labels$BioProject],
            alpha = 0.2) +
  
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = colors_K) +
  theme_minimal() +
  labs(title = "Bter 1.0",
       subtitle = "Genetic Ancestry Proportions",
       x = NULL,
       y = "Ancestry Proportion",
       fill = "Ancestry Component") +
  theme(axis.text.x = element_blank(), # Hide individual IID labels
        axis.ticks.x = element_blank()) +
  coord_cartesian(ylim = c(0, 1.07), clip = "off")

# Conditionally add the population demarcation layers
if (nrow(population_segments) > 0) {
  p_individual <- p_individual +
    geom_rect(data = population_segments,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE,
              fill = rep(c("grey90", "grey80"), length.out = nrow(population_segments)),
              alpha = 0.3) +
    geom_text(data = population_segments,
              aes(x = label_x, y = 1.05, label = paste0(Population, " Bees")),
              inherit.aes = FALSE, vjust = 0, fontface = "bold", size = 4.5)
}

# Add the BioProject labels at the bottom (without color)
p_individual <- p_individual +
  geom_text(data = bioproject_labels,
            aes(x = midpoint, y = -0.05, label = BioProject),
            inherit.aes = FALSE,
            angle = 45,
            hjust = 1,
            size = 3) +
  theme(plot.margin = unit(c(1, 1, 3, 1), "lines")) # Add margin for labels

print(p_individual)

