library(readxl)
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

# loading flux sampling output from matlab
setwd("/media/alexis/hdd2/objetivo_1_tesis_doctoral/matlab_scripts_genomas")
#sampling <- fread("samples_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_200000_modelseed_ids.csv")
sampling <- fread("samples_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_200000_modelseed_Asegurando_que_PT33_cresca_en_monocultivo.csv")

sampling_rxn_names <- fread("samples_EcComRXNS_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_modelseed_Asegurando_que_PT33_cresca_en_monocultivo.csv", header = FALSE)
#sampling_rxn_names <- fread("samples_EcComRXNS_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_modelseed_ids.csv", header = FALSE)

modelseed_names <- fread("modelseed_compounds.tsv", header = TRUE)
colnames(modelseed_names)
str(modelseed_names)

# Exploring tables to understand merge
colnames(sampling)[1:2]
colnames(sampling_rxn_names)[1:2]
rownames(sampling)[1:2]
rownames(sampling_rxn_names)[1:2]

# Set the row names of sampling using the first column of sampling_rxn_names
row.names(sampling) <- sampling_rxn_names[[1]]
rownames(sampling)[1:2]

# Identify rows containing "EX_" in their name
ex_rows <- grepl("EX_", row.names(sampling))

# Identify rows containing "biomass" in their name
biomass_rows <- grepl("biomass", row.names(sampling))

# Combine both conditions using the logical operator OR (|)
combined_rows <- ex_rows | biomass_rows

# Filter the original "sampling" table to obtain the rows that meet the following conditions
sampling_filtered <- sampling[combined_rows, ]

# Set the row names of the filtered dataset to the original row names that meet the condition
row.names(sampling_filtered) <- row.names(sampling)[combined_rows]

# Print selected rows
rownames(sampling_filtered)

# Convert row names to character strings
rownames(sampling_filtered) <- as.character(rownames(sampling_filtered))
rownames(sampling_filtered)[1:2]
rownames(sampling_filtered)

rm(sampling) # Use Restart R from Session menu to recover some RAM
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Transform modelseed id names to metabolite names

library(stringr)

# Step 1: Ensure that the row names are strings
rownames(sampling_filtered) <- as.character(rownames(sampling_filtered))
modelseed_names$id <- as.character(modelseed_names$id)

# Step 2: Create a function to look up and get the corresponding names in modelseed_names
get_metabolite_name <- function(metabolite_id) {
  index <- match(metabolite_id, modelseed_names$id)
  if (!is.na(index)) {
    return(modelseed_names$name[index])
  } else {
    return(NA) # Returns NA if the identifier is not found in modelseed_names
  }
}

# Step 3: Replace the row names in sampling_filtered with the corresponding names
unique_metabolite_names <- sapply(rownames(sampling_filtered), function(row_name) {
  bacterium_name <- str_extract(row_name, "^[^_]+") # Extract the name of the bacterium (part before the first underscore)
  metabolite_id <- str_extract(row_name, "(?<=cpd)\\d+") # # Extract the identifier only (cpd...)
  
  if (!is.na(metabolite_id)) {
    metabolite_name <- get_metabolite_name(paste0("cpd", metabolite_id))
    if (!is.na(metabolite_name)) {
      return(paste0(bacterium_name, "_", metabolite_name))
    }
  }
  
  return(row_name) # Return the original row name if the identifier is not found
})

# Resolve duplicates and make row names unique
unique_metabolite_names <- make.unique(unique_metabolite_names, sep = "_")

# Replace the row names in sampling_filtered with the corresponding unique names
sampling_filtered$metabolite_name <- unique_metabolite_names
rownames(sampling_filtered) <- sampling_filtered$metabolite_name
sampling_filtered$metabolite_name <- NULL # We remove the temporary column from metabolite_name

# Verification: Shows the first sampling_filtered rows after replacement
rownames(sampling_filtered)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Making a correlation matrix from external metabolite between bifidobacterium and clostridium from flux sampling

# Select only Clostridium and Bifidobacterium related fluxes
nombres_flujos <- rownames(sampling_filtered)
# Create an index based on whether the stream name contains "Clostridium" or "Bifidobacterium"
indice_flujos <- grep("Clostridium|Bifidobacterium", nombres_flujos)

# Filtering the data and keeping the original names of each flux
flujos_bacterias <- sampling_filtered[indice_flujos, ]
rownames(flujos_bacterias) <- nombres_flujos[indice_flujos]

rownames(flujos_bacterias)

# Transpose the data set
flujos_transpuestos <- t(flujos_bacterias)

# Separate Clostridium and Bifidobacterium streams
flujos_clostridium <- flujos_bacterias[, grepl("Clostridium_", colnames(flujos_bacterias))]
flujos_bifidobacterium <- flujos_bacterias[, grepl("Bifidobacterium_", colnames(flujos_bacterias))]

# debuggin
# Make sure that the transaction_flows have column names
colnames(flujos_transpuestos) <- rownames(flujos_bacterias)

# Separate data for Clostridium and Bifidobacterium using the correct names
flujos_clostridium <- flujos_transpuestos[, grepl("Clostridium_", colnames(flujos_transpuestos))]
flujos_bifidobacterium <- flujos_transpuestos[, grepl("Bifidobacterium_", colnames(flujos_transpuestos))]

# Verify that there is data before proceeding
if (ncol(flujos_clostridium) > 0 && ncol(flujos_bifidobacterium) > 0) {
  # Calculate standard deviations to verify the data
  sds_clostridium <- apply(flujos_clostridium, 2, sd)
  sds_bifidobacterium <- apply(flujos_bifidobacterium, 2, sd)
  print(sds_clostridium)
  print(sds_bifidobacterium)
  
  # If there are no standard deviations from zero, calculate the correlation
  if (!any(sds_clostridium == 0) && !any(sds_bifidobacterium == 0)) {
    matriz_correlacion_cruzada <- cor(flujos_clostridium, flujos_bifidobacterium, method = "spearman")
    print(matriz_correlacion_cruzada)
  } else {
    print("There are columns with zero standard deviation. Correlation cannot be calculated.")
  }
} else {
  print("Insufficient data in one or both data sets to calculate the correlation.")
}

# Calculates the standard deviation for the Clostridium and Bifidobacterium assemblies
sds_clostridium <- apply(flujos_clostridium, 2, sd)
sds_bifidobacterium <- apply(flujos_bifidobacterium, 2, sd)

# Filters columns where the standard deviation is non-zero
flujos_clostridium <- flujos_clostridium[, sds_clostridium > 0]
flujos_bifidobacterium <- flujos_bifidobacterium[, sds_bifidobacterium > 0]

# Check again if there is enough data to calculate the correlation
if (ncol(flujos_clostridium) > 0 && ncol(flujos_bifidobacterium) > 0) {
  matriz_correlacion_cruzada <- cor(flujos_clostridium, flujos_bifidobacterium, method = "spearman")
  print(matriz_correlacion_cruzada)
} else {
  print("There are not enough columns with variability to calculate the correlation.")
}

matriz_correlacion_cruzada

# Convert the correlation matrix to long format
matriz_larga <- melt(matriz_correlacion_cruzada)

# Plotting the heatmap
ggplot(data = matriz_larga, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + # Usar geom_tile para crear el heatmap
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 45)) +
  labs(fill = "Correlation") # Cambiar la etiqueta de la leyenda de color

# improve order for better heatmap visualization

# Extract the names of the metabolites of Bifidobacterium and Clostridium
nombres_metabolitos_bifido <- gsub("Bifidobacterium_", "", colnames(matriz_correlacion_cruzada))
nombres_metabolitos_clostridium <- gsub("Clostridium_", "", rownames(matriz_correlacion_cruzada))

# Identify common names
nombres_comunes <- intersect(nombres_metabolitos_bifido, nombres_metabolitos_clostridium)

# Create an order vector for rows and columns based on common names
orden_filas <- match(nombres_comunes, nombres_metabolitos_clostridium)
orden_columnas <- match(nombres_comunes, nombres_metabolitos_bifido)

# Create the new ordered correlation matrix
matriz_correlacion_ordenada <- matriz_correlacion_cruzada[orden_filas, orden_columnas]

# Make sure that row and column names are correct
rownames(matriz_correlacion_ordenada) <- paste("Clostridium", nombres_comunes, sep="_")
colnames(matriz_correlacion_ordenada) <- paste("Bifidobacterium", nombres_comunes, sep="_")

# Convert the correlation matrix to long format
matriz_larga_ordenada <- melt(matriz_correlacion_ordenada)

# Plotting the heatmap
p4 <- ggplot(data = matriz_larga_ordenada, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black", size = 0.1) +  # Define el color y el tamaño de las líneas entre los cuadros
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Spearman Correlation") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"), # Fondo blanco para el gráfico
    panel.background = element_rect(fill = "white", colour = "white"), # Fondo blanco para el panel
    axis.text.x = element_text(angle = 70, vjust = 1, hjust=1, size = 8, color = "black", face = "bold"), # Hacer texto del eje X más grande y negro
    axis.text.y = element_text(angle = 0, vjust = 1, hjust=1, size = 8, color = "black", face = "bold"), # Hacer texto del eje Y más grande y negro
    axis.title.x = element_blank(), # Quitar "Var1"
    axis.title.y = element_blank(), # Quitar "Var2"
    legend.position = "right" # Posición de la leyenda
  ) +
  labs(x = NULL, y = NULL) # Quitar etiquetas de los ejes

p4

ggsave("matriz_correlacion_sampling_spearman_200000_final.png", p4, width = 10, height = 8, dpi = 600)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# P-values
# Esta función toma dos matrices de datos y realiza pruebas de correlación por pares
calcular_correlaciones <- function(flujos1, flujos2, method = "spearman") {
  n <- ncol(flujos1)
  m <- ncol(flujos2)
  matriz_p_values <- matrix(NA, n, m)
  matriz_cor_values <- matrix(NA, n, m)
  
  for (i in 1:n) {
    for (j in 1:m) {
      test <- cor.test(flujos1[,i], flujos2[,j], method = method)
      matriz_p_values[i,j] <- test$p.value
      matriz_cor_values[i,j] <- test$estimate
    }
  }
  
  # Ajuste de los valores p para múltiples pruebas, si se desea
  # matriz_p_values_ajustados <- p.adjust(matriz_p_values, method = "bonferroni")
  
  list(correlation = matriz_cor_values, p.value = matriz_p_values)
}

# Llamar a la función con tus datos
resultados_correlacion <- calcular_correlaciones(flujos_clostridium, flujos_bifidobacterium)

# Obtener las matrices de correlación y valores p
matriz_correlacion <- resultados_correlacion$correlation
matriz_p_values <- resultados_correlacion$p.value

# Asignar los nombres a las filas y columnas de la matriz de correlación
rownames(matriz_correlacion) <- colnames(flujos_clostridium)
colnames(matriz_correlacion) <- colnames(flujos_bifidobacterium)

# Convertir la matriz de p-values en un vector
vector_p_values <- as.vector(t(matriz_p_values))  # Usamos la transposición para asegurar que se desenrolle por columnas adecuadamente

# Aplicar p.adjust al vector de p-values
vector_p_values_ajustados <- p.adjust(vector_p_values, method = "bonferroni")

# Reconvertir el vector ajustado en una matriz con las mismas dimensiones que la matriz de p-values original
matriz_p_values_ajustados <- matrix(vector_p_values_ajustados, nrow = nrow(matriz_p_values), ncol = ncol(matriz_p_values), byrow = TRUE)

# Verificar las dimensiones para asegurarse de que todo está correcto
print(dim(matriz_p_values_ajustados))

# Imprimir una pequeña parte para verificar
print(head(matriz_p_values_ajustados[1:5, 1:5]))


# Asegurar que los dataframes contienen los nombres correctos
df_correlacion <- as.data.frame(as.table(as.matrix(matriz_correlacion)))
names(df_correlacion) <- c("Metabolito_Clostridium", "Metabolito_Bifidobacterium", "Correlacion_Spearman")
colnames(df_correlacion)
df_correlacion

# Asignar nombres a las matrices p_values y p_values_ajustados antes de convertir a dataframes
rownames(matriz_p_values) <- colnames(flujos_clostridium)
colnames(matriz_p_values) <- colnames(flujos_bifidobacterium)
rownames(matriz_p_values_ajustados) <- colnames(flujos_clostridium)
colnames(matriz_p_values_ajustados) <- colnames(flujos_bifidobacterium)

# Convertir las matrices en dataframes
df_p_values <- as.data.frame(as.table(as.matrix(matriz_p_values)))
names(df_p_values) <- c("Metabolito_Clostridium", "Metabolito_Bifidobacterium", "P_Value")
df_p_values

df_p_values_ajustados <- as.data.frame(as.table(as.matrix(matriz_p_values_ajustados)))
names(df_p_values_ajustados) <- c("Metabolito_Clostridium", "Metabolito_Bifidobacterium", "P_Value_Ajustado")
df_p_values_ajustados

# Verificar que los nombres ahora son correctos
head(df_p_values)
head(df_p_values_ajustados)

# Fusionar los dataframes
df_final <- merge(df_correlacion, df_p_values, by = c("Metabolito_Clostridium", "Metabolito_Bifidobacterium"))
df_final <- merge(df_final, df_p_values_ajustados, by = c("Metabolito_Clostridium", "Metabolito_Bifidobacterium"))

head(df_final)
head(matriz_larga_ordenada)

colnames(matriz_larga_ordenada)
colnames(df_final)

# Ordenar y mostrar las primeras filas para verificar
#df_final <- df_final[order(df_final$P_Value_Ajustado), ]
#head(df_final)

# Añadir una nueva columna al dataframe donde las correlaciones no significativas serán establecidas a 0
#df_final$Correlacion_Filtrada <- ifelse(df_final$P_Value_Ajustado < 0.05, df_final$Correlacion_Spearman, 0)

# Primero, asegúrate de que los nombres de las columnas de 'matriz_larga_ordenada' coincidan con los de 'df_final' para poder mergearlas correctamente
colnames(matriz_larga_ordenada) <- c("Metabolito_Clostridium", "Metabolito_Bifidobacterium", "Correlacion_Spearman")

# Luego, realiza el merge utilizando las columnas de los nombres de metabolitos como llaves
df_combinado <- merge(matriz_larga_ordenada, df_final, by = c("Metabolito_Clostridium", "Metabolito_Bifidobacterium", "Correlacion_Spearman"))

# Ahora 'df_combinado' contendrá las correlaciones junto con los valores p y p ajustados correspondientes, en el orden de 'matriz_larga_ordenada'
head(df_combinado)
head(matriz_larga_ordenada)

# Exportar el dataframe combinado a un archivo CSV
write.csv(df_combinado, "Correlation_matrix_melted.csv", row.names = FALSE)
