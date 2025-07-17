# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

# Define file paths
q_file <- "C:/Users/User/Documents/Masterthesis/Results_ORG/admixture/K6Q.txt"
fam_file <- "C:/Users/User/Documents/Masterthesis/Results_ORG/admixture/famfile.txt"
pop_file <- "C:/Users/User/Documents/Masterthesis/Results_ORG/admixture/sample_type_information.txt"
output_pdf <- "C:/Users/User/Documents/Masterthesis/Results_ORG/admixture/admixture_plot_K6_combined_demarcated.pdf"

# Read input files
q_data <- read.table(q_file)
fam_data <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
pop_data <- read.table(pop_file, header = FALSE, stringsAsFactors = FALSE)

# Rename columns for clarity
colnames(fam_data) <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")
colnames(pop_data) <- c("IID", "Population")

# Clean and group population names
pop_data$Population <- ifelse(toupper(trimws(pop_data$Population)) == "WILD", "Wild", "Commercial")

# Merge Q data with FAM and population info
q_data$IID <- fam_data$IID
admix_data <- merge(q_data, pop_data, by = "IID")
admix_data$Population <- factor(admix_data$Population, levels = c("Wild", "Commercial"))

# Order individuals by population and component 1 (just for consistency in bars)
admix_data <- admix_data %>% arrange(Population, desc(V1))

# Convert to long format for ggplot2
admix_data_long <- melt(admix_data, id.vars = c("IID", "Population"),
                        variable.name = "Ancestry_Component", value.name = "Proportion")
admix_data_long$IID <- factor(admix_data_long$IID, levels = admix_data$IID)

# Define population segments for plotting
ordered_iids <- levels(admix_data_long$IID)
population_segments <- admix_data %>%
  group_by(Population) %>%
  summarize(
    xmin = min(as.numeric(factor(IID, levels = ordered_iids))) - 0.5,
    xmax = max(as.numeric(factor(IID, levels = ordered_iids))) + 0.5,
    label_x = (min(as.numeric(factor(IID, levels = ordered_iids))) +
                 max(as.numeric(factor(IID, levels = ordered_iids)))) / 2
  )

# Generate plot
colors_K6 <- viridis(n = 6, option = "D")
p <- ggplot(admix_data_long, aes(x = IID, y = Proportion, fill = Ancestry_Component)) +
  geom_rect(data = population_segments,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = c("grey90", "grey80"), alpha = 0.3) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colors_K6) +
  geom_text(data = population_segments,
            aes(x = label_x, y = 1.05, label = paste0(Population, " Bees")),
            inherit.aes = FALSE, vjust = 0, fontface = "bold", size = 4.5) +
  theme_minimal() +
  labs(
    title = "Admixture Model",
    subtitle = "Genetic Ancestry Proportions (K=6)",
    x = NULL,
    y = "Ancestry Proportion",
    fill = "Populations K"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = unit(c(1, 1, 0.5, 0.5), "cm")
  ) +
  coord_cartesian(ylim = c(0, 1.07), clip = "off")
p
