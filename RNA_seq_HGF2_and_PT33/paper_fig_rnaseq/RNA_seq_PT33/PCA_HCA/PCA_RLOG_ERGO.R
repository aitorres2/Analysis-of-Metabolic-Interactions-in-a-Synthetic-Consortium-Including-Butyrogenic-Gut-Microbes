library(tidyverse)
library(readxl)
library("ggplot2")
library("dplyr")
library("tidyr")
library("readr")
library("purrr")
library("tibble")
library("stringr")
library("forcats")
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
library(stringr)
library(DESeq2)
library(apeglm)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(PoiClaClu)
library(apeglm)
library(purrr)
library(readr)
library(pheatmap)
library(clusterProfiler)
library(gridExtra)
library(cowplot)
library(ggalt)
library(ggforce)

setwd("/home/alexis/Desktop/RESULTADOS_FINALES_ERGO/kallisto_TRABAJADO/PT33/descargados_ERGO_deseq2/PCA")

counts <- read.csv2("counts_PT33_only_kallisto.csv",  row.names = 1, header = TRUE, sep=",")
# Verificar la estructura de la tabla importada
str(counts)

min_count <- 10 # Puedes ajustar este número según tus necesidades.
transformed_counts <- counts %>%
  replace(is.na(.), 0) %>%
  filter(rowSums(.) >= min_count)

# Crear un vector que represente las muestras en el orden correcto
PT33_samples <- c("Consortium_Early_R1", "Consortium_Early_R2", "Consortium_Early_R3", 
                  "Consortium_Late_R1", "Consortium_Late_R2", "Consortium_Late_R3", 
                  "PT33_Early_R1", "PT33_Early_R2", "PT33_Early_R3", 
                  "PT33_Late_R1", "PT33_Late_R2", "PT33_Late_R3")

# Crear un vector que represente los grupos correspondientes a cada muestra
PT33_groups <- c("Consortium_Early", "Consortium_Early", "Consortium_Early", 
                 "Consortium_Late", "Consortium_Late", "Consortium_Late", 
                 "PT33_Early", "PT33_Early", "PT33_Early", 
                 "PT33_Late", "PT33_Late", "PT33_Late")


# Crear el objeto colData
PT33_colData <- data.frame(exp = as.factor(PT33_groups),
                           row.names = PT33_samples)
PT33_colData

# Counts transformados de PT33
PT33_all_counts <- transformed_counts %>%
  filter(rowSums(.) >= 10)
# Verificar la estructura de la tabla filtrada
str(PT33_all_counts)

# Verificar los nombres de fila y columna en PT33_all_counts
rownames(PT33_all_counts)
colnames(PT33_all_counts)

# Verificar los nombres de fila y columna en PT33_colData
rownames(PT33_colData)
colnames(PT33_colData)

# Obtener los nuevos nombres de las columnas desde las filas de PT33_colData
new_column_names <- rownames(PT33_colData)

# Asignar los nuevos nombres de las columnas a PT33_all_counts
colnames(PT33_all_counts) <- new_column_names

# Crear el objeto DESeqDataSet utilizando el countData y el colData
design_PT33_all_counts <- DESeqDataSetFromMatrix(countData = PT33_all_counts,
                                                 colData = PT33_colData,
                                                 design = ~ exp)

# Normalizar counts, librearia en cantidad de reads (siempre se deben normalizar los datos de rna-seq); esta normalizacion es unica de deseq
dds_PT33_all_counts <- estimateSizeFactors(design_PT33_all_counts)
# View
sizeFactors(dds_PT33_all_counts)
# Extract normalized counts
counts_dds_PT33_all_counts <- counts(dds_PT33_all_counts, normalized = TRUE)

# Realizar la transformación rlog considerando las condiciones experimentales
rld_dds_PT33_all_counts <- rlog(dds_PT33_all_counts, blind = FALSE)
PT33_all_pca <- plotPCA(rld_dds_PT33_all_counts, intgroup = "exp")
PT33_all_pca
# Realizar la transformación rlog considerando las condiciones experimentales
#rld_dds_PT33_all_counts <- rlog(dds_PT33_all_counts, blind = TRUE)

rld_mat_dds_PT33_all_counts <- assay(rld_dds_PT33_all_counts) 
rld_cor_mat_dds_PT33_all_counts <- cor(rld_mat_dds_PT33_all_counts, method = "spearman")
#rld_cor_mat_dds_PT33_all_counts <- cor(rld_mat_dds_PT33_all_counts) # pearson correlation

p1 <- pheatmap(rld_cor_mat_dds_PT33_all_counts,
               clustering_distance_cols = "euclidean",
               clustering_distance_rows = "euclidean",
               clustering_method = "complete",
               legend = TRUE) # Añadir título a la leyenda
p1

#png("PT33_pheatmap_pre_filt.png", bg = "white", width = 20, height = 15, units = 'cm', res = 500)
#pheatmap(rld_cor_mat_dds_PT33_all_counts)
#dev.off()

# CUales exactamente la replica que esta molestando:
# Generar el gráfico PCA y utilizar "condition" como intgroup
PT33_all_pcaData <- plotPCA(rld_dds_PT33_all_counts, intgroup = "exp", returnData = TRUE)

# Agregar columna "replicate" a partir de los nombres de las muestras
PT33_all_pcaData$replicate <- PT33_all_pcaData$name

# Crear columna de replicado solo con valores 'replicate1', 'replicate2', 'replicate3'
PT33_all_pcaData$replicate <- paste0('Replicate', substr(PT33_all_pcaData$replicate, nchar(PT33_all_pcaData$replicate), nchar(PT33_all_pcaData$replicate)))

# Calcular el porcentaje de varianza explicada por cada componente
percentVar <- round(100 * attr(PT33_all_pcaData, "percentVar"), 2)

# Asignar formas automáticamente a las categorías únicas
shape_palette <- c(15, 16, 17)  # ahora tenemos 4 formas

# Definir unique_exp_categories
unique_exp_categories <- unique(PT33_all_pcaData$exp)

# Asignar colores automáticamente a las categorías únicas
color_palette <- colorRampPalette(c("#1f78b4", "#33a02c", "#e31a1c", "#9370DB"))(length(unique_exp_categories))

#HGF2_all_pcaData$group <- HGF2_all_pcaData$exp # Crea un nuevo factor que agrupe los datos

# Asignar colores automáticamente a las categorías únicas
color_palette <- colorRampPalette(c("#1f78b4", "#33a02c", "#e31a1c", "#9370DB"))(length(unique_exp_categories))

p2 <- ggplot(PT33_all_pcaData, aes(x = PC1, y = PC2, color = exp)) +
  geom_point(aes(shape = replicate), size = 4) +
  geom_mark_ellipse(aes(fill = exp), 
                    asp_ratio = 1, 
                    expand = 0.02,  # Valor más pequeño para hacer la elipse más ajustada
                    alpha = 0.2,
                    show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  ) +
  scale_color_manual("Experiment", values = color_palette) +
  scale_shape_manual("Replicates", values = shape_palette) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  scale_fill_manual("Experiment", values = color_palette)

# Guardar gráfico PCA previo a la remoción de outliers
#if (!dir.exists("results")) {
#  dir.create("results")
#}

p2

#ggsave("PT33_PCA_2_pre_filt.png", p2, bg = "white", dpi = 500, 
#       width = 30, height = 20, units = 'cm')

#########################################################################################################################################################################3
# Eliminar muestra
colnames(transformed_counts)
min_count <- 10 # Puedes ajustar este número según tus necesidades.
transformed_counts <- counts %>%
  select(-PT33_Early_R3) %>%
  replace(is.na(.), 0) %>%
  filter(rowSums(.) >= min_count)

# Crear un vector que represente las muestras en el orden correcto
PT33_samples <- c("Consortium_Early_R1", "Consortium_Early_R2", "Consortium_Early_R3",
                  "Consortium_Late_R1", "Consortium_Late_R2", "Consortium_Late_R3", 
                  "PT33_Early_R1", "PT33_Early_R2", 
                  "PT33_Late_R1", "PT33_Late_R2", "PT33_Late_R3")

# Crear un vector que represente los grupos correspondientes a cada muestra
PT33_groups <- c("Consortium_Early", "Consortium_Early", "Consortium_Early",
                 "Consortium_Late", "Consortium_Late", "Consortium_Late",
                 "PT33_Early", "PT33_Early", 
                 "PT33_Late", "PT33_Late", "PT33_Late")

# Crear el objeto colData
PT33_colData <- data.frame(exp = as.factor(PT33_groups),
                           row.names = PT33_samples)
PT33_colData

# Counts transformados de PT33
PT33_all_counts <- transformed_counts %>%
  filter(rowSums(.) >= 10)
# Verificar la estructura de la tabla filtrada
str(PT33_all_counts)

# Verificar los nombres de fila y columna en PT33_all_counts
rownames(PT33_all_counts)
colnames(PT33_all_counts)

# Verificar los nombres de fila y columna en PT33_colData
rownames(PT33_colData)
colnames(PT33_colData)

# Obtener los nuevos nombres de las columnas desde las filas de PT33_colData
new_column_names <- rownames(PT33_colData)

# Asignar los nuevos nombres de las columnas a PT33_all_counts
colnames(PT33_all_counts) <- new_column_names

# Crear el objeto DESeqDataSet utilizando el countData y el colData
design_PT33_all_counts <- DESeqDataSetFromMatrix(countData = PT33_all_counts,
                                                 colData = PT33_colData,
                                                 design = ~ exp)

# Normalizar counts, librearia en cantidad de reads (siempre se deben normalizar los datos de rna-seq); esta normalizacion es unica de deseq
dds_PT33_all_counts <- estimateSizeFactors(design_PT33_all_counts)
# View
sizeFactors(dds_PT33_all_counts)
# Extract normalized counts
counts_dds_PT33_all_counts <- counts(dds_PT33_all_counts, normalized = TRUE)

# Realizar la transformación rlog considerando las condiciones experimentales
rld_dds_PT33_all_counts <- rlog(dds_PT33_all_counts, blind = FALSE)
PT33_all_pca <- plotPCA(rld_dds_PT33_all_counts, intgroup = "exp")
PT33_all_pca
# Realizar la transformación rlog considerando las condiciones experimentales
#rld_dds_PT33_all_counts <- rlog(dds_PT33_all_counts, blind = TRUE)

rld_mat_dds_PT33_all_counts <- assay(rld_dds_PT33_all_counts) 
rld_cor_mat_dds_PT33_all_counts <- cor(rld_mat_dds_PT33_all_counts, method = "spearman")
p3 <- pheatmap(rld_cor_mat_dds_PT33_all_counts,
               clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean",
               clustering_method = "complete")
p3

#png("PT33_pheatmap_post_filt.png", bg = "white", width = 30, height = 20, units = 'cm', res = 500)
#pheatmap(rld_cor_mat_dds_PT33_all_counts)
#dev.off()

##########################################
# CUales exactamente la replica que esta molestando:
# Generar el gráfico PCA y utilizar "condition" como intgroup
PT33_all_pcaData <- plotPCA(rld_dds_PT33_all_counts, intgroup = "exp", returnData = TRUE)

# Agregar columna "replicate" a partir de los nombres de las muestras
PT33_all_pcaData$replicate <- PT33_all_pcaData$name

# Crear columna de replicado solo con valores 'replicate1', 'replicate2', 'replicate3'
PT33_all_pcaData$replicate <- paste0('Replicate', substr(PT33_all_pcaData$replicate, nchar(PT33_all_pcaData$replicate), nchar(PT33_all_pcaData$replicate)))

# Calcular el porcentaje de varianza explicada por cada componente
percentVar <- round(100 * attr(PT33_all_pcaData, "percentVar"), 2)

# Asignar formas automáticamente a las categorías únicas
shape_palette <- c(15, 16, 17)  # ahora tenemos 4 formas

# Definir unique_exp_categories
unique_exp_categories <- unique(PT33_all_pcaData$exp)

# Asignar colores automáticamente a las categorías únicas
color_palette <- colorRampPalette(c("#1f78b4", "#33a02c", "#e31a1c", "#9370DB"))(length(unique_exp_categories))

#HGF2_all_pcaData$group <- HGF2_all_pcaData$exp # Crea un nuevo factor que agrupe los datos

# Asignar colores automáticamente a las categorías únicas
color_palette <- colorRampPalette(c("#1f78b4", "#33a02c", "#e31a1c", "#9370DB"))(length(unique_exp_categories))

p4 <- ggplot(PT33_all_pcaData, aes(x = PC1, y = PC2, color = exp)) +
  geom_point(aes(shape = replicate), size = 4) +
  geom_mark_ellipse(aes(fill = exp), 
                    asp_ratio = 1, 
                    expand = 0.02, 
                    alpha = 0.2,
                    show.legend = F) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  ) +
  scale_color_manual("Experiment", values = color_palette) +
  scale_shape_manual("Replicates", values = shape_palette) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  scale_fill_manual("Experiment", values = color_palette)
p4

# Guardar gráfico PCA previo a la remoción de outliers
#if (!dir.exists("results")) {
#  dir.create("results")
#}

p4 

#ggsave("PT33_PCA_2_post_filt.png", ggsave_pca, bg = "white", dpi = 500, 
#       width = 30, height = 20, units = 'cm')

# Asegúrate de que las figuras p1, p2, p3, y p4 están correctamente generadas

# Guardar cada heatmap como una imagen
pheatmap_save <- function(p, filename) {
  png(filename, width = 7, height = 7, units = 'in', res = 600)
  print(p)
  dev.off()
}

# Reemplaza 'p1', 'p2', 'p3', 'p4' con tus objetos de gráficos reales
pheatmap_save(p1, "heatmap1.png")
pheatmap_save(p2, "heatmap2.png")
pheatmap_save(p3, "heatmap3.png")
pheatmap_save(p4, "heatmap4.png")

# Instalar y cargar las bibliotecas si aún no están instaladas
if (!require("grid")) install.packages("grid")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("png")) install.packages("png")

library(grid)
library(gridExtra)
library(png)

# Leer las imágenes
img1 <- readPNG("heatmap1.png")
img2 <- readPNG("heatmap2.png")
img3 <- readPNG("heatmap3.png")
img4 <- readPNG("heatmap4.png")

# Convertir imágenes en grobs
grob1 <- rasterGrob(img1, interpolate = TRUE)
grob2 <- rasterGrob(img2, interpolate = TRUE)
grob3 <- rasterGrob(img3, interpolate = TRUE)
grob4 <- rasterGrob(img4, interpolate = TRUE)

# Combinar las imágenes
combined_plot <- grid.arrange(
  grob1, grob2,
  grob3, grob4,
  ncol = 2,
  heights = c(1, 1) # Ajustar si es necesario
)

# Guardar la imagen combinada
ggsave("combined_heatmaps.png", combined_plot, width = 14, height = 14, units = 'in', dpi = 600)


