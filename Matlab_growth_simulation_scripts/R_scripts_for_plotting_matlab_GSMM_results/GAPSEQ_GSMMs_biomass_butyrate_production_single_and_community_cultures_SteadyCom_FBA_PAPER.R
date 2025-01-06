library(readxl)
library("ggplot2")
library("dplyr")
library("tidyr")
library("readr")
library("purrr")
library("tibble")
library("stringr")
library("forcats")
library(tidyverse)
library(stringr)
library("pheatmap")
library("grid")
library(data.table)
library("writexl")
library(ggplot2)
library(matrixStats)
library(reshape2)
library(grid)
library(gridExtra)
#install.packages("devtools")
library(devtools)
#install.packages("patchwork")
library(patchwork)
#install.packages("cowplot")
library(cowplot)

# procesando sampling
setwd("/home/alexis/Desktop/github_script_finales_alexis/Objective_1_Metabolic_interactions_insilico/Simulation_of_growth/R_scripts_for_plotting_matlab_GSMM_results")

diseño_consorcios_data <- read_excel("design_consortium_modelseed_vFINAL.xlsx", sheet = 1, col_names = TRUE)
str(diseño_consorcios_data)
row.names(diseño_consorcios_data)
colnames(diseño_consorcios_data)

diseño_consorcios_df <- as.data.frame(diseño_consorcios_data)
str(diseño_consorcios_df)
summary(diseño_consorcios_df)

##----------------------------------------------------------------------------------------------------------------------------
# Subset the data to only include the relevant consortium and bacteria
mono_subset <- diseño_consorcios_df %>% 
  filter(Consortium_type %in% grep("Mono cultures", diseño_consorcios_df$Consortium_type, value = TRUE))
col_order_1 <- mono_subset$Consortium
col_order_1 <- unique(col_order_1)
mono_subset$Consortium <- factor(mono_subset$Consortium, levels = col_order_1)

#----------------------------------------------------------------------------------------------------------------------------------
innocuum_subset <- diseño_consorcios_df %>% 
  filter(Consortium %in% grep("HGF2", diseño_consorcios_df$Consortium, value = TRUE)) 
col_order_2 <- innocuum_subset$Consortium
col_order_2 <- unique(col_order_2)
innocuum_subset$Consortium <- factor(innocuum_subset$Consortium, levels = col_order_2)
# remove rows with "Mono cultures" in Consortium_type column (modify original dataframe)
innocuum_subset <- innocuum_subset[innocuum_subset$Consortium_type != "Mono cultures", ]

#----------------------------------------------------------------------------------------------------------------------------------
symbiosum_subset <- diseño_consorcios_df %>% 
  filter(Consortium %in% grep("WAL", diseño_consorcios_df$Consortium, value = TRUE)) 
col_order_3 <- symbiosum_subset$Consortium
col_order_3 <- unique(col_order_3)
symbiosum_subset$Consortium <- factor(symbiosum_subset$Consortium, levels = col_order_3)
# remove rows with "Mono cultures" in Consortium_type column (modify original dataframe)
symbiosum_subset <- symbiosum_subset[symbiosum_subset$Consortium_type != "Mono cultures", ]

#----------------------------------------------------------------------------------------------------------------------------------
saccharolyticum_subset <- diseño_consorcios_df %>% 
  filter(Consortium %in% grep("M62", diseño_consorcios_df$Consortium, value = TRUE)) 
col_order_4 <- saccharolyticum_subset$Consortium
col_order_4 <- unique(col_order_4)
saccharolyticum_subset$Consortium <- factor(saccharolyticum_subset$Consortium, levels = col_order_4)
# remove rows with "Mono cultures" in Consortium_type column (modify original dataframe)
saccharolyticum_subset <- saccharolyticum_subset[saccharolyticum_subset$Consortium_type != "Mono cultures", ]

#----------------------------------------------------------------------------------------------------------------------------------
mono_subset$Butyrate <- as.numeric(mono_subset$Butyrate)
innocuum_subset$Butyrate <- as.numeric(innocuum_subset$Butyrate)
symbiosum_subset$Butyrate <- as.numeric(symbiosum_subset$Butyrate)
saccharolyticum_subset$Butyrate <- as.numeric(saccharolyticum_subset$Butyrate)

#----------------------------------------------------------------------------------------------------------------------------------
mono_subset$Consortium <- gsub("\\d+-Mono-", "", mono_subset$Consortium)
innocuum_subset$Consortium <- gsub("\\d+-Mono-|\\d+-Co_|\\d+-Tri_|\\d+-Tetra_", "", innocuum_subset$Consortium)
symbiosum_subset$Consortium <- gsub("\\d+-Mono-|\\d+-Co_|\\d+-Tri_|\\d+-Tetra_", "", symbiosum_subset$Consortium)
saccharolyticum_subset$Consortium <- gsub("\\d+-Mono-|\\d+-Co_|\\d+-Tri_|\\d+-Tetra_", "", saccharolyticum_subset$Consortium)

mono_subset$Butyrate <- as.numeric(mono_subset$Butyrate)
innocuum_subset$Butyrate <- as.numeric(innocuum_subset$Butyrate)
symbiosum_subset$Butyrate <- as.numeric(symbiosum_subset$Butyrate)
saccharolyticum_subset$Butyrate <- as.numeric(saccharolyticum_subset$Butyrate)
#----------------------------------------------------------------------------------------------------------------------------------
# Define a new color palette
my_palette <- c("Clostridium_sp_HGF2" = "#FF6600", 
                "Lc_symbiosum_WAL_14673" = "#008080", 
                "Clostridium_sp_M62/1" = "#A2CBEF",
                "B_animalis_lactis_PT33" = "#3498DB", 
                "B_thetaiotaomicron_VPI_5482" = "#AA6AC8", 
                "L_paracasei_M38" = "#2ECC71")

# Define una función personalizada para cambiar los nombres de las facetas
custom_labeller <- labeller(
  Consortium_type = c(
    "Mono cultures" = "Mono-cultures",
    "Consortium Clostridium innocuum" = "Clostridium sp. HGF2",
    "Consortium Clostridium symbiosum" = "Lc. symbiosum WAL14673",
    "Consortium Clostridium saccharolyticum" = "Clostridium sp. M62"
  )
)

# Create separate plots for Biomass_bacterium and Butyrate
plot1 <- ggplot(mono_subset, aes(x = Consortium, y = Biomass_bacterium, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  labs(x = "", y = "Growth rate [1/h]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 0.35) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

plot1

plot2 <- ggplot(mono_subset, aes(x = Consortium, y = Butyrate, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(data = subset(mono_subset, Butyrate > 0), aes(label = round(Butyrate, 2)), vjust = -0.5, size = 3, color = "black") +
  labs(x = "", y = "Butyrate [mmol/(gDW·h)]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 3.65) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

# Arrange the plots side by side
grid.arrange(plot1, plot2, ncol = 2)

# Create separate plots for Biomass_bacterium and Butyrate
plot3 <- ggplot(innocuum_subset, aes(x = Consortium, y = Biomass_bacterium, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  labs(x = "", y = "Growth rate [1/h]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 0.35) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "#FF6600"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

plot4 <- ggplot(innocuum_subset, aes(x = Consortium, y = Butyrate, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(data = subset(innocuum_subset, Butyrate > 0), aes(label = round(Butyrate, 2)), vjust = -0.5, size = 3, color = "black") +
  labs(x = "", y = "Butyrate [mmol/(gDW·h)]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 3.65) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "#FF6600"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

# Arrange the plots side by side
grid.arrange(plot3, plot4, ncol = 2)

# Create separate plots for Biomass_bacterium and Butyrate
plot5 <- ggplot(symbiosum_subset, aes(x = Consortium, y = Biomass_bacterium, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  labs(x = "", y = "Growth rate [1/h]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 0.35) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "#008080"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

plot6 <- ggplot(symbiosum_subset, aes(x = Consortium, y = Butyrate, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(data = subset(symbiosum_subset, Butyrate > 0), aes(label = round(Butyrate, 2)), vjust = -0.5, size = 3, color = "black") +
  labs(x = "", y = "Butyrate [mmol/(gDW·h)]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 3.65) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "#008080"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

# Arrange the plots side by side
grid.arrange(plot5, plot6, ncol = 2)


# Create separate plots for Biomass_bacterium and Butyrate
plot9 <- ggplot(saccharolyticum_subset, aes(x = Consortium, y = Biomass_bacterium, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  labs(x = "", y = "Growth rate [1/h]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 0.35) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "#A2CBEF"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

plot10 <- ggplot(saccharolyticum_subset, aes(x = Consortium, y = Butyrate, fill = Bacteria)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(data = subset(saccharolyticum_subset, Butyrate > 1), aes(label = round(Butyrate, 2)),
            vjust = -0.5, size = 3, color = "black") +
  labs(x = "", y = "Butyrate [mmol/(gDW·h)]", fill = "Species") +
  scale_fill_manual(values = my_palette) +
  ylim(0, 3.65) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold", size = 14),
    axis.title.y = element_text(size = 14),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold", colour = "#A2CBEF"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 20)
  ) +
  facet_wrap(~ Consortium_type, labeller = custom_labeller) +  # Usa la función labeller directamente
  guides(fill = FALSE)

plot10

# Arrange the plots side by side
grid.arrange(plot9, plot10, ncol = 2)

library(cowplot)

# plots biomasas 1,3,5,7,9
# plots butiarto: 2,4,6,8,10
# Combine the plots into a single plot without a shared legend
# Combine the plots into a single plot without a shared legend
# Combine the plots into a single plot without a shared legend
# Ajusta el tamaño de fuente de los títulos de los gráficos
# Ajusta el tamaño de fuente de los títulos de los gráficos y la leyenda
font_size_titles <- 24
font_size_legend <- 12

# Ajusta el ancho relativo de los gráficos en combined_plot
combined_plot <- plot_grid(plot1, plot3, plot5, plot9,
                           plot2, plot4, plot6, plot10, ncol = 4,
                           align = "h", axis = "tb", 
                           labels = c("A1", "B1", "C1", "D1",
                                      "A2", "B2", "C2", "D2"),
                           rel_widths = c(1, 1, 1, 1, 1, 1, 1, 1))

# Crea la leyenda compartida usando el primer plot
shared_legend <- get_legend(
  ggplot(mono_subset, aes(x = Consortium, y = Biomass_bacterium, fill = Bacteria)) +
    geom_bar(stat = "identity") +
    labs(x = "Mono cultures", y = "Growth rate [1/h]", fill = "Species") +
    scale_fill_manual(values = my_palette) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 80, vjust = 0.95, hjust = 0.95, face = "bold"),
      panel.spacing = unit(0.1, "cm"),
      strip.text.x = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = font_size_legend),  # Ajusta el tamaño de la fuente de la leyenda
      legend.title = element_text(size = font_size_legend),  # Ajusta el tamaño de la fuente del título de la leyenda
      legend.margin = margin(r = -140, unit = "pt")  # Reduce el margen derecho de la leyenda
    ))

# Combina los gráficos y la leyenda en un único plot
final_plot <- plot_grid(
  combined_plot,
  shared_legend,
  NULL,  # Espacio en blanco
  ncol = 3,
  rel_widths = c(1.2, 0.05, 0.15)
)

# Aumenta el tamaño de fuente de los títulos de los gráficos
final_plot <- final_plot +
  theme(plot.title = element_text(size = font_size_titles))  # Ajusta el tamaño de la fuente del título de los gráficos

# Ajusta la orientación de la leyenda
final_plot <- final_plot +
  theme(
    legend.box = "horizontal",  # Cambia la orientación de la caja de la leyenda a horizontal
    legend.margin = margin(t = 0, unit = "pt"),  # Ajusta el margen superior de la leyenda
    legend.spacing = unit(0, "pt")  # Elimina el espacio entre la leyenda y el gráfico
  )

# Muestra el gráfico final
final_plot

# Save figure
ggsave("final_plot_gapseq_2_PAPER.pdf", plot = final_plot, width = 20, height = 12, units = "in", dpi = 320, bg = "white")
ggsave("final_plot_gapseq_2_PAPER.png", plot = final_plot, width = 20, height = 12, units = "in", dpi = 320, bg = "white")

# Biomass [h-1]
# Butyrate [mmol gDW-1 h-1]


