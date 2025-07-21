library(ggplot2)
library(dplyr)
library(readxl)
setwd("C:/Users/User/Documents/Masterthesis/Results_ORG/admixture")
calculating_K_value <- read_excel("calculating_K_value.xlsx")
View(calculating_K_value)
cv_data <- calculating_K_value
cv_data$K_Value <- as.numeric(cv_data$K_Value)
cv_data <- cv_data %>% arrange(K_Value)

min_k_value <- cv_data$K_Value[which.min(cv_data$CV_Error)]

p <- ggplot(cv_data, aes(x = K_Value, y = CV_Error)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "darkblue", size = 3) +
  geom_vline(xintercept = min_k_value, linetype = "dashed", color = "red", size = 1) +
  geom_point(data = filter(cv_data, K_Value == min_k_value),
             aes(x = K_Value, y = CV_Error), color = "red", size = 5, shape = 21, fill = "red") +
  labs(
    title = "ADMIXTURE Model Selection Across K Values",
    subtitle = "Reference Genome: Bter_1.0",
    x = "Number of Ancestral Populations (K)",
    y = "Cross-Validation Error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),                     # Remove all grid lines
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  ) +
  scale_x_continuous(breaks = unique(cv_data$K_Value))
p
