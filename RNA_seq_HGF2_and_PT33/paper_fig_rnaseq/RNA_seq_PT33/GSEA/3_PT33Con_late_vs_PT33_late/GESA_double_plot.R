library(tidyverse)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(grid)  # Asegurarse de incluir esta librería

setwd("/home/alexis/Desktop/RESULTADOS_FINALES_ERGO/kallisto_TRABAJADO/PT33/descargados_ERGO_deseq2/GSEA_gage/3_PT33Con_late_vs_PT33_late")

# Carga tus datos
upregulated_ergo <- read_csv("gage (ergo): Consortium_Late - PT8_Late_upregulated.csv")
downregulated_ergo <- read_csv("gage (ergo): Consortium_Late - PT8_Late_downregulated.csv")

str(upregulated_ergo)
str(downregulated_ergo)

# Reemplazar '-' por '_' en la columna 'Name'
upregulated_ergo <- upregulated_ergo %>% mutate(Name = str_replace_all(Name, "-", "_"))
downregulated_ergo <- downregulated_ergo %>% mutate(Name = str_replace_all(Name, "-", "_"))
downregulated_ergo$Name
# Filtrar y remover duplicados basándose en la columna 'Name'
upregulated_ergo_unique <- upregulated_ergo %>% distinct(Name, .keep_all = TRUE)
downregulated_ergo_unique <- downregulated_ergo %>% distinct(Name, .keep_all = TRUE)

# Filtrar los datos
umbral_p_valor <- 0.05
umbral_q_valor <- 0.05
upregulated_filtered <- upregulated_ergo_unique %>% filter(`P Value` <= umbral_p_valor, `Q Value` <= umbral_q_valor)
downregulated_filtered <- downregulated_ergo_unique %>% filter(`P Value` <= umbral_p_valor, `Q Value` <= umbral_q_valor)

# Categorías a excluir
categorias_excluir <- c("PRPPGLNIMP.ANA", "TAT_pathway_secreted_proteins", "TAT_pathway",
                        "Bacterial_Endogenous_Protein_Quality_Control_and_Degradation",
                        "Secretion", "Pyrimidine_salvage_pathways", "Sec_Dependent_Secretion_Stages",
                        "Photosynthesis_in_cyanobacteria", "Heteropolysaccharide_biosynthesis",
                        "Cytosolic_reactions_of_nitrogen_metabolism",
                        "Electron_transport_in_cyanobacterial_tylakoid_membrane",
                        "Sugar_precursors_of_streptococcal_capsular_and_exopolysaccharides",
                        "Purine_salvage_pathways", "Bacterial_Translation_Factors",
                        "Photosynthesis_and_respiration_in_cyanobacteria", "Bacterial_Photosynthesis",
                        "Reserve_carbohydrate_biosynthesis", "Bacterial_Protein_Fate",
                        "Bacterial_aerobic_respiration_proton_transport_and_oxidative_phosphorylation",
                        "Bacterial_respiration_proton_transport_and_oxidative_phosphorylation",
                        "Glycolysis_versions_and_bypasses", "Bacterial_Ribosomal_Complex",
                        "Proton_motive_cycle", "Bacterial_LSU", "Bacterial_SSU", "Anaerobic_respiratory_chain_proton_transport_and_oxidative_phosphorylation")


# Filtrar y eliminar las categorías seleccionadas
upregulated_filtered_clean <- upregulated_filtered %>% 
  filter(!Name %in% categorias_excluir) %>%
  arrange(`Q Value`)  # Orden ascendente para upregulated

downregulated_filtered_clean <- downregulated_filtered %>% 
  filter(!Name %in% categorias_excluir) %>%
  arrange(`Q Value`)  # Orden descendente para downregulated

library(scales)  # Para usar la función trans_new para transformaciones logarítmicas

# Colores para la escala
color_low <- "#E41A1C"  # Un tono pastel claro de rojo
color_high <- "#377EB8"  # Un tono pastel claro de azul

# Función para crear un gráfico con barras ordenadas por Q Value de forma ascendente
create_plot <- function(data, title) {
  # Reordena los datos en función del Q Value de forma ascendente
  data <- data %>%
    mutate(Name = reorder(Name, `Q Value`, FUN = function(x) -max(abs(x))))
  
  ggplot(data, aes(x = Name, y = `Stat Mean`, fill = `Q Value`)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_gradientn(colors = c(color_low, color_high),
                         trans = 'log10',  # Aplicar transformación logarítmica a la escala de colores
                         na.value = "grey50",  # Color para los valores NA
                         guide = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                                 label.theme = element_text(size = 18))) + # Aumenta el tamaño del texto de la leyenda
    labs(title = title, x = "", y = "Stat Mean", fill = "Log10(Q-Value)") +
    theme_minimal() +
    theme(
      text = element_text(face = "bold", size = 20), # Aumenta el tamaño general del texto
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 18), # Aumenta el tamaño del texto en el eje X
      axis.text.y = element_text(size = 18, face = "bold"), # Aumenta el tamaño del texto en el eje Y
      axis.title = element_text(face = "bold", size = 18), # Aumenta el tamaño de los títulos de los ejes
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # Aumenta el tamaño del título del gráfico
      legend.title = element_text(face = "bold", size = 13), # Aumenta el tamaño del título de la leyenda
      legend.text = element_text(face = "bold", size = 18), # Aumenta el tamaño del texto de la leyenda
      plot.margin = margin(t = 1, r = 2, b = 1, l = 4, unit = "cm") # Ajusta los márgenes del gráfico
    )
}

# Aplicar la función create_plot a tus datos filtrados y limpios
upregulated_plot <- create_plot(upregulated_filtered_clean, "Upregulated genes")
downregulated_plot <- create_plot(downregulated_filtered_clean, "Downregulated genes")

# Organizar en dos columnas con el mismo ancho de barras y tamaño de texto y agregar título principal
titulo_principal <- "PT33_Consortium_Late_vs_PT33_Single_Late"
combined_plot <- grid.arrange(upregulated_plot, downregulated_plot, ncol = 2,
                              top = textGrob(titulo_principal, gp = gpar(fontsize = 24, fontface = "bold"))); # Aumenta el tamaño del título principal


# Guardar la imagen
ggsave("combined_ergo_plot_1.png", combined_plot, width = 24, height = 15, dpi = 600)

# Combinar las tablas añadiendo una columna para identificar la fuente
combined_table <- bind_rows(
  upregulated_filtered_clean %>% mutate(Source = "Upregulated"),
  downregulated_filtered_clean %>% mutate(Source = "Downregulated")
)

# Especifica el nombre y ruta del archivo donde deseas guardar la tabla
path_to_file <- "/home/alexis/Desktop/github_script_finales_alexis/Objective_3_cocultures_vs_monocultures_RNAseq/RNA_seq/RNA_seq_Bi_animalis_lactis_PT33/GSEA/3_PT33Con_late_vs_PT33_late/combined_table_double.csv"

# Exportar la tabla a un archivo CSV
write.csv(combined_table, path_to_file, row.names = FALSE)





