library(readxl)
library("ggplot2")
library("dplyr")
library("tidyr")
library(tidyverse)
library(data.table)
library("writexl")
library(reshape2)
library(factoextra)
library(corrplot)
library(Hmisc)

setwd("/home/alexis/Desktop/github_script_finales_alexis/Objective_2_SCFAs_metabolic_interactions_invitro/Invitro_cultures_growth_&_SCFAs")

df <- read_excel("screening_only_BPBs.xlsx", sheet = 1, col_names = TRUE, na = "")
#names(df)

# Delete rows with "M38" in Experiment column
df <- subset(df, !grepl("M38", Experiment))

df

numeric_data <- df[, c("Growth [OD620]", "Butyrate [mmol/L]", "Acetate [mmol/L]", 
                       "Lactate [mmol/L]", "Propionate [mmol/L]", "Succinate [mmol/L]")]

numeric_data

# Calculate the correlation matrix and p-values.
corr_results <- rcorr(as.matrix(numeric_data), type = "spearman")

# Convert correlation matrices and p-values to long format
melted_cor <- melt(corr_results$r)
melted_p <- melt(corr_results$P)

# Combine the melted matrices
cor_data <- cbind(melted_cor, P_Value = melted_p$value)

# Add a new column for significance levels
cor_data$Significance <- ifelse(cor_data$P_Value < 0.001, '***', 
                                ifelse(cor_data$P_Value < 0.01, '**', 
                                       ifelse(cor_data$P_Value < 0.05, '*', '')))

# Creating the heat map with ggplot2
p1 <- ggplot(data = cor_data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limit = c(-1, 1), name="Spearman\nCorrelation") +
  geom_text(aes(label = ifelse(value == 1, '', paste0(sprintf("%.2f", value), Significance))), size = 4) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12), # Ajustar el tamaño aquí
        axis.text.y = element_text(size = 12), # Ajustar el tamaño aquí
        axis.title = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

p1

# Save the heatmap to a file
ggsave("heatmap_correlations_pvalues_invitro.png", p1, width = 10, height = 8, dpi = 600)


