# Cargar bibliotecas necesarias
library(readODS)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Configurar el directorio de trabajo
setwd("/home/alexis/Desktop/github_script_finales_alexis/Objective_3_cocultures_vs_monocultures_RNAseq/RNA_seq/paper_fig_rnaseq")
dir_path <- "/home/alexis/Desktop/github_script_finales_alexis/Objective_3_cocultures_vs_monocultures_RNAseq/RNA_seq/paper_fig_rnaseq"

# Lista de archivos ODS
files <- list.files(dir_path, pattern = "\\.ods$", full.names = TRUE)

# Cargar todas las tablas en una lista
tables <- lapply(files, function(file) read_ods(file, sheet = 1))
names(tables) <- gsub("\\.ods$", "", basename(files))

# Ajustar el rango del eje X para consistencia
x_limits <- c(-5, 5)

# Limpiar las tablas para eliminar valores fuera del rango o NA
clean_table <- function(data, x_limits) {
  data <- data %>%
    filter(!is.na(Log2FC)) %>%           # Eliminar filas con NA en Log2FC
    filter(Log2FC >= x_limits[1] & Log2FC <= x_limits[2])  # Filtrar valores dentro del rango
  return(data)
}

# Aplicar limpieza a todas las tablas
tables <- lapply(tables, clean_table, x_limits = x_limits)

# Función para gráficos personalizados sin tags
plot_function <- function(data, title, show_x_title = FALSE) {
  ggplot(data, aes(x = Log2FC, y = Gene_name, fill = Log2FC > 0)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("TRUE" = "salmon", "FALSE" = "lightblue")) +
    scale_x_continuous(limits = c(-5, 5), expand = c(0, 0)) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = if (show_x_title) "Log2FC" else NULL, y = NULL) +
    theme(
      axis.text.y = element_text(size = 12, hjust = 1),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.line.x = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(20, 15, 15, 15),  # Márgenes amplios
      legend.position = "none",
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank()
    ) +
    coord_fixed(ratio = 1.2)  # Proporción fija
}

# Crear gráficos con el ajuste fino
plots <- list(
  A = plot_function(tables[["HGF2_Co_8_vs_HGF2_Mono_8"]], "HGF2_Co_8_vs_HGF2_Mono_8"),
  D = plot_function(tables[["PT33_Co_8_vs_PT33_Mono_8"]], "PT33_Co_8_vs_PT33_Mono_8"),
  B = plot_function(tables[["HGF2_Co_45_vs_HGF2_Mono_45"]], "HGF2_Co_45_vs_HGF2_Mono_45"),
  E = plot_function(tables[["PT33_Co_45_vs_PT33_Mono_45"]], "PT33_Co_45_vs_PT33_Mono_45"),
  C = plot_function(tables[["HGF2_Mono_45_vs_HGF2_Mono_8"]], "HGF2_Mono_45_vs_HGF2_Mono_8"),
  F = plot_function(tables[["PT33_Mono_45_vs_PT33_Mono_8"]], "PT33_Mono_45_vs_PT33_Mono_8")
)

# Combinar gráficos en formato 2x3 con espaciado adecuado
final_plot_2x3 <- wrap_plots(
  list(
    plots$A, plots$D,  # Primera fila
    plots$B, plots$E,  # Segunda fila
    plots$C, plots$F   # Tercera fila
  ),
  ncol = 2,
  heights = c(1, 1, 1),
  guides = "collect"
)

# Mostrar el gráfico ajustado
print(final_plot_2x3)

# Guardar el gráfico ajustado sin tags
ggsave(
  "final_plot_2x3_no_tags.png",
  plot = final_plot_2x3,
  width = 15,
  height = 18,
  dpi = 600,
  bg = "white"
)


