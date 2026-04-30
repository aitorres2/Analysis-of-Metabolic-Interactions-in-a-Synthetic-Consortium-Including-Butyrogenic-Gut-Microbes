# =====================================================================
# HEATMAPS COMBINADOS HGF2 + PT33 - DESDE TPM (sin DESeq2)
# Version basica para principiantes
# =====================================================================
# Flujo:
#  1. Leer 12 TSVs de kallisto por organismo (uno por muestra)
#  2. Armar matriz TPM (filas = genes, columnas = muestras)
#  3. Filtrar TPM > 1 en al menos 2 muestras
#  4. Transformar log2(TPM + 1)
#  5. Leer anotaciones funcionales del genoma
#  6. Leer lista de genes de interes con sus categorias (.fasta)
#  7. Seleccionar genes por categoria (con reglas especiales)
#  8. Z-score por fila (gen) -> matriz centrada
#  9. Heatmap con pheatmap + anotaciones
# 10. Combinar HGF2 + PT33 con patchwork (1x2 y 2x1)
# =====================================================================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(ggplotify)


# ---------------------------------------------------------------------
# Carpetas
# ---------------------------------------------------------------------
BASE <- "/media/alexis/hdd2/objetivo_1_tesis_doctoral/matlab_scripts_genomas/0_ACTUALIZACION_paper_simulacion/scripts_rna_seq/0_actualizados/01_12_25_scripts_completos_rnaseq_correccion"
OUT  <- file.path(BASE, "deseq2_alex_kallisto/figuras_combinadas_HGF2_PT33")

setwd(OUT)


# =====================================================================
# ============================  HGF2  =================================
# =====================================================================

# ---------------------------------------------------------------------
# 1. Leer los 12 archivos TSV de kallisto (uno por muestra)
# ---------------------------------------------------------------------
dir_hgf2 <- file.path(BASE, "output_kallisto", "HGF2")

t1  <- read.table(file.path(dir_hgf2, "Coculture_PT33_HGF2_Early_Replicate_1_HGF2/Coculture_PT33_HGF2_Early_Replicate_1_HGF2.tsv"), header = TRUE, sep = "\t")
t2  <- read.table(file.path(dir_hgf2, "Coculture_PT33_HGF2_Early_Replicate_2_HGF2/Coculture_PT33_HGF2_Early_Replicate_2_HGF2.tsv"), header = TRUE, sep = "\t")
t3  <- read.table(file.path(dir_hgf2, "Coculture_PT33_HGF2_Early_Replicate_3_HGF2/Coculture_PT33_HGF2_Early_Replicate_3_HGF2.tsv"), header = TRUE, sep = "\t")
t4  <- read.table(file.path(dir_hgf2, "Coculture_PT33_HGF2_Late_Replicate_1_HGF2/Coculture_PT33_HGF2_Late_Replicate_1_HGF2.tsv"),   header = TRUE, sep = "\t")
t5  <- read.table(file.path(dir_hgf2, "Coculture_PT33_HGF2_Late_Replicate_2_HGF2/Coculture_PT33_HGF2_Late_Replicate_2_HGF2.tsv"),   header = TRUE, sep = "\t")
t6  <- read.table(file.path(dir_hgf2, "Coculture_PT33_HGF2_Late_Replicate_3_HGF2/Coculture_PT33_HGF2_Late_Replicate_3_HGF2.tsv"),   header = TRUE, sep = "\t")
t7  <- read.table(file.path(dir_hgf2, "Monoculture_Clostridium_HGF2_Early_Replicate_1_HGF2/Monoculture_Clostridium_HGF2_Early_Replicate_1_HGF2.tsv"), header = TRUE, sep = "\t")
t8  <- read.table(file.path(dir_hgf2, "Monoculture_Clostridium_HGF2_Early_Replicate_2_HGF2/Monoculture_Clostridium_HGF2_Early_Replicate_2_HGF2.tsv"), header = TRUE, sep = "\t")
t9  <- read.table(file.path(dir_hgf2, "Monoculture_Clostridium_HGF2_Early_Replicate_3_HGF2/Monoculture_Clostridium_HGF2_Early_Replicate_3_HGF2.tsv"), header = TRUE, sep = "\t")
t10 <- read.table(file.path(dir_hgf2, "Monoculture_Clostridium_HGF2_Late_Replicate_1_HGF2/Monoculture_Clostridium_HGF2_Late_Replicate_1_HGF2.tsv"),   header = TRUE, sep = "\t")
t11 <- read.table(file.path(dir_hgf2, "Monoculture_Clostridium_HGF2_Late_Replicate_2_HGF2/Monoculture_Clostridium_HGF2_Late_Replicate_2_HGF2.tsv"),   header = TRUE, sep = "\t")
t12 <- read.table(file.path(dir_hgf2, "Monoculture_Clostridium_HGF2_Late_Replicate_3_HGF2/Monoculture_Clostridium_HGF2_Late_Replicate_3_HGF2.tsv"),   header = TRUE, sep = "\t")


# ---------------------------------------------------------------------
# 2. Construir la matriz TPM (filas = genes, columnas = muestras)
# ---------------------------------------------------------------------
tpm_hgf2 <- data.frame(
  Consortium_Early_R1 = t1$tpm,
  Consortium_Early_R2 = t2$tpm,
  Consortium_Early_R3 = t3$tpm,
  Consortium_Late_R1  = t4$tpm,
  Consortium_Late_R2  = t5$tpm,
  Consortium_Late_R3  = t6$tpm,
  HGF2_Early_R1       = t7$tpm,
  HGF2_Early_R2       = t8$tpm,
  HGF2_Early_R3       = t9$tpm,
  HGF2_Late_R1        = t10$tpm,
  HGF2_Late_R2        = t11$tpm,
  HGF2_Late_R3        = t12$tpm
)
rownames(tpm_hgf2) <- t1[[1]]   # primera columna = target_id (el nombre real puede ser X.e.target_id por un bug del header)

str(tpm_hgf2)
head(rownames(tpm_hgf2))
colnames(tpm_hgf2)
head(tpm_hgf2)

# ---------------------------------------------------------------------
# 3. Filtrar genes con TPM > 1 en al menos 2 muestras
# ---------------------------------------------------------------------
tpm_hgf2 <- tpm_hgf2[rowSums(tpm_hgf2 > 1) >= 2, ]

head(tpm_hgf2)

# ---------------------------------------------------------------------
# 4. Transformacion log2(TPM + 1)  -> matriz para usar en heatmap
# ---------------------------------------------------------------------
log2_hgf2 <- as.matrix(log2(tpm_hgf2 + 1))

head(log2_hgf2)

# ---------------------------------------------------------------------
# 5. Leer anotaciones funcionales del genoma de HGF2
# ---------------------------------------------------------------------
annot_hgf2 <- read_csv(
  file.path(BASE, "deseq2_alex_kallisto/HGF2/Genome_annotation/Clostridium_HGF2_annotations.csv"),
  show_col_types = FALSE
)

head(annot_hgf2)
str(annot_hgf2)

# ---------------------------------------------------------------------
# 6. Leer lista de genes de interes (.fasta con categorias en ###)
# ---------------------------------------------------------------------
lineas_hgf2 <- readLines(
  file.path(BASE, "deseq2_alex_kallisto/HGF2/Genome_annotation/lista_genes_interes_RNAseq_paper_HGF2.fasta")
)

# Loop simple: cada linea que empieza con ### es una categoria nueva
# y cada linea que empieza con RCHGF es un gen de la categoria actual
categoria_actual <- NA
genes_id   <- c()
genes_cat  <- c()

for (linea in lineas_hgf2) {
  if (startsWith(linea, "###")) {
    categoria_actual <- trimws(sub("^#+", "", linea))
  } else if (startsWith(linea, "RCHGF") && !is.na(categoria_actual)) {
    gene_id <- strsplit(linea, "\\s+")[[1]][1]
    genes_id  <- c(genes_id, gene_id)
    genes_cat <- c(genes_cat, categoria_actual)
  }
}

genes_interes_hgf2 <- tibble(gene_id = genes_id, category = genes_cat)
categorias_hgf2    <- unique(genes_interes_hgf2$category)


# ---------------------------------------------------------------------
# 7. Seleccionar genes por categoria (reglas especiales para algunas)
# ---------------------------------------------------------------------

# Mapeo muestra -> condicion (sin el sufijo _R1/_R2/_R3)
muestras_hgf2   <- colnames(log2_hgf2)
condiciones_map_hgf2 <- tibble(
  Replicate = muestras_hgf2,
  Condition = sub("_R[0-9]+$", "", muestras_hgf2)
)
condiciones_unicas_hgf2 <- unique(condiciones_map_hgf2$Condition)

head(muestras_hgf2)
str(muestras_hgf2)

head(condiciones_unicas_hgf2)
str(condiciones_unicas_hgf2)

# Lista donde vamos guardando los genes seleccionados de cada categoria
genes_por_cat_hgf2 <- list()


# =====================================================================
# CATEGORIA ESPECIAL 1: Glycolysis_genes
# Regla: 1 gen por funcion eggNOG (el de mayor expresion media)
# =====================================================================
glyco_disponibles <- genes_interes_hgf2 %>%
  filter(category == "Glycolysis_genes",
         gene_id %in% rownames(log2_hgf2)) %>%
  pull(gene_id)

genes_por_cat_hgf2[["Glycolysis_genes"]] <- tibble(
    gene_id    = glyco_disponibles,
    media_expr = rowMeans(log2_hgf2[glyco_disponibles, , drop = FALSE])
  ) %>%
  left_join(annot_hgf2 %>% select(id, eggnog_Description),
            by = c("gene_id" = "id")) %>%
  mutate(funcion = if_else(
    is.na(eggnog_Description) | eggnog_Description %in% c("", "NA"),
    gene_id,
    tolower(eggnog_Description)
  )) %>%
  group_by(funcion) %>%
  slice_max(media_expr, n = 1, with_ties = FALSE) %>%
  pull(gene_id)


# =====================================================================
# CATEGORIA ESPECIAL 2: Sugar transport genes
# Regla: top 5 z-score promedio por cada condicion experimental
# =====================================================================
sugar_disponibles <- genes_interes_hgf2 %>%
  filter(category == "Sugar transport genes",
         gene_id %in% rownames(log2_hgf2)) %>%
  pull(gene_id)

# Z-score por fila (gen) sobre los genes disponibles
mat_z_sugar <- t(scale(t(log2_hgf2[sugar_disponibles, , drop = FALSE])))

genes_por_cat_hgf2[["Sugar transport genes"]] <- as_tibble(mat_z_sugar, rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "Replicate", values_to = "Z") %>%
  left_join(condiciones_map_hgf2, by = "Replicate") %>%
  group_by(gene_id, Condition) %>%
  summarise(avg_z = mean(Z, na.rm = TRUE), .groups = "drop") %>%
  group_by(Condition) %>%
  slice_max(avg_z, n = 5) %>%
  pull(gene_id) %>%
  unique()


# =====================================================================
# CATEGORIA ESPECIAL 3: Sporulation genes
# Regla: top 10 por expresion media (sin z-score)
# =====================================================================
spor_disponibles <- genes_interes_hgf2 %>%
  filter(category == "Sporulation genes",
         gene_id %in% rownames(log2_hgf2)) %>%
  pull(gene_id)

genes_por_cat_hgf2[["Sporulation genes"]] <- tibble(
    gene_id = spor_disponibles,
    media   = rowMeans(log2_hgf2[spor_disponibles, , drop = FALSE])
  ) %>%
  slice_max(media, n = 10) %>%
  pull(gene_id)


# =====================================================================
# RESTO DE CATEGORIAS: usar todos los genes disponibles (sin filtros)
# =====================================================================
categorias_normales_hgf2 <- setdiff(
  categorias_hgf2,
  c("Glycolysis_genes", "Sugar transport genes", "Sporulation genes")
)

for (cat in categorias_normales_hgf2) {
  genes_cat <- genes_interes_hgf2 %>%
    filter(category == cat, gene_id %in% rownames(log2_hgf2)) %>%
    pull(gene_id)

  if (length(genes_cat) > 0) {
    genes_por_cat_hgf2[[cat]] <- genes_cat
  }
}


# ---------------------------------------------------------------------
# 8. Z-score por gen -> matriz para el heatmap
# ---------------------------------------------------------------------
genes_todos_hgf2 <- unique(unlist(genes_por_cat_hgf2))
mat_hgf2         <- log2_hgf2[genes_todos_hgf2, , drop = FALSE]
mat_z_hgf2       <- t(scale(t(mat_hgf2)))   # z-score por fila (gen)

head(genes_todos_hgf2)
head(mat_hgf2)
head(mat_z_hgf2)

# Quitar filas con NA, NaN o Inf (genes con varianza 0 -> scale() divide entre 0)
filas_validas    <- apply(mat_z_hgf2, 1, function(x) all(is.finite(x)))
mat_z_hgf2       <- mat_z_hgf2[filas_validas, , drop = FALSE]
genes_todos_hgf2 <- genes_todos_hgf2[filas_validas]

head(filas_validas)
head(mat_z_hgf2)
head(genes_todos_hgf2)

cat("HGF2 - genes despues del z-score:", nrow(mat_z_hgf2), "\n")
cat("HGF2 - hay valores no finitos?:", any(!is.finite(mat_z_hgf2)), "\n")


# ---------------------------------------------------------------------
# 9. Etiquetas (gene_id - funcion) y anotaciones para el heatmap
# ---------------------------------------------------------------------

# --- Funcion: para cada gen, sacar la mejor anotacion funcional ---
# Cascada: GenBank.Function -> ERGO.Function -> eggnog_Description -> id
mapear_anotacion_hgf2 <- function(gene_ids) {
  sapply(gene_ids, function(gid) {
    fila <- annot_hgf2[annot_hgf2$id == gid, ]
    if (nrow(fila) == 0) return(gid)

    anot <- fila$GenBank.Function[1]
    if (is.na(anot) || anot == "" || anot == "NA") anot <- fila$ERGO.Function[1]
    if (is.na(anot) || anot == "" || anot == "NA") anot <- fila$eggnog_Description[1]
    if (is.na(anot) || anot == "" || anot == "NA") return(gid)

    return(paste0(gid, " - ", anot))
  })
}

# --- Funcion: anotacion especifica para genes de glicolisis (KEGG Orthology) ---
mapear_anotacion_glyco_hgf2 <- function(gene_ids) {
  sapply(gene_ids, function(gid) {
    fila <- annot_hgf2[annot_hgf2$id == gid, ]
    if (nrow(fila) == 0) return(gid)

    anot <- fila$Kegg.Orthology[1]
    if (is.na(anot) || anot == "" || anot == "NA") return(gid)

    # Quitar el prefijo "Kxxxxx - " del KEGG Orthology
    anot_limpia <- sub("^K[0-9]+ - ", "", anot)
    return(paste0(gid, " - ", anot_limpia))
  })
}


# Construir las etiquetas de cada gen (distinto si es glicolisis)
glyco_hgf2 <- if ("Glycolysis_genes" %in% names(genes_por_cat_hgf2)) {
  genes_por_cat_hgf2[["Glycolysis_genes"]]
} else {
  character(0)
}

# Construimos las etiquetas con un loop simple para garantizar
# que el resultado sea un character vector (sapply puede devolver lista)
etiquetas_hgf2 <- character(length(genes_todos_hgf2))
for (i in seq_along(genes_todos_hgf2)) {
  gid <- genes_todos_hgf2[i]
  if (gid %in% glyco_hgf2) {
    etiquetas_hgf2[i] <- as.character(mapear_anotacion_glyco_hgf2(gid))
  } else {
    etiquetas_hgf2[i] <- as.character(mapear_anotacion_hgf2(gid))
  }
}

# Hacer las etiquetas unicas (por si hay funciones repetidas)
etiquetas_hgf2 <- make.unique(etiquetas_hgf2, sep = "_")
rownames(mat_z_hgf2) <- etiquetas_hgf2


# Anotacion de filas: a que pathway pertenece cada gen
row_annot_hgf2 <- data.frame(
  Pathway = character(length(genes_todos_hgf2)),
  row.names = etiquetas_hgf2,
  stringsAsFactors = FALSE
)

for (cat in names(genes_por_cat_hgf2)) {
  genes_de_cat <- genes_por_cat_hgf2[[cat]]

  if (cat == "Glycolysis_genes") {
    etiquetas_cat <- mapear_anotacion_glyco_hgf2(genes_de_cat)
  } else {
    etiquetas_cat <- mapear_anotacion_hgf2(genes_de_cat)
  }
  etiquetas_cat <- make.unique(etiquetas_cat, sep = "_")

  filas_match <- rownames(row_annot_hgf2) %in% etiquetas_cat
  row_annot_hgf2$Pathway[filas_match] <- cat
}


# Anotacion de columnas: condicion experimental
col_annot_hgf2 <- data.frame(
  Condition = factor(condiciones_map_hgf2$Condition),
  row.names = condiciones_map_hgf2$Replicate
)


# ---------------------------------------------------------------------
# 10. Colores para anotaciones
# ---------------------------------------------------------------------
colores_condicion_hgf2 <- c(
  "Consortium_Early" = "#E41A1C",
  "Consortium_Late"  = "#377EB8",
  "HGF2_Early"       = "#4DAF4A",
  "HGF2_Late"        = "#984EA3"
)

n_cat_hgf2 <- length(genes_por_cat_hgf2)
if (n_cat_hgf2 <= 12) {
  colores_pathway_hgf2 <- brewer.pal(n_cat_hgf2, "Set3")
} else {
  colores_pathway_hgf2 <- colorRampPalette(brewer.pal(12, "Set3"))(n_cat_hgf2)
}
names(colores_pathway_hgf2) <- names(genes_por_cat_hgf2)

ann_colors_hgf2 <- list(
  Condition = colores_condicion_hgf2,
  Pathway   = colores_pathway_hgf2
)


# ---------------------------------------------------------------------
# 11. Generar heatmap de HGF2
# ---------------------------------------------------------------------
p_complete_HGF2 <- pheatmap(
  mat               = mat_z_hgf2,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  annotation_col    = col_annot_hgf2,
  annotation_row    = row_annot_hgf2,
  annotation_colors = ann_colors_hgf2,
  color             = colorRampPalette(c("green", "black", "red"))(100),
  fontsize          = 10,
  fontsize_row      = 6,
  fontsize_col      = 10,
  angle_col         = "90",
  border_color      = NA,
  scale             = "none",
  cellwidth         = 11,
  cellheight        = 7
)




# =====================================================================
# ============================  PT33  =================================
# =====================================================================

# ---------------------------------------------------------------------
# 1. Leer los 12 archivos TSV de kallisto (uno por muestra)
# ---------------------------------------------------------------------
dir_pt33 <- file.path(BASE, "output_kallisto", "PT33")

s1  <- read.table(file.path(dir_pt33, "Coculture_PT33_HGF2_Early_Replicate_1_PT33/Coculture_PT33_HGF2_Early_Replicate_1_PT33.tsv"), header = TRUE, sep = "\t")
s2  <- read.table(file.path(dir_pt33, "Coculture_PT33_HGF2_Early_Replicate_2_PT33/Coculture_PT33_HGF2_Early_Replicate_2_PT33.tsv"), header = TRUE, sep = "\t")
s3  <- read.table(file.path(dir_pt33, "Coculture_PT33_HGF2_Early_Replicate_3_PT33/Coculture_PT33_HGF2_Early_Replicate_3_PT33.tsv"), header = TRUE, sep = "\t")
s4  <- read.table(file.path(dir_pt33, "Coculture_PT33_HGF2_Late_Replicate_1_PT33/Coculture_PT33_HGF2_Late_Replicate_1_PT33.tsv"),   header = TRUE, sep = "\t")
s5  <- read.table(file.path(dir_pt33, "Coculture_PT33_HGF2_Late_Replicate_2_PT33/Coculture_PT33_HGF2_Late_Replicate_2_PT33.tsv"),   header = TRUE, sep = "\t")
s6  <- read.table(file.path(dir_pt33, "Coculture_PT33_HGF2_Late_Replicate_3_PT33/Coculture_PT33_HGF2_Late_Replicate_3_PT33.tsv"),   header = TRUE, sep = "\t")
s7  <- read.table(file.path(dir_pt33, "Monoculture_Bi_animalis_lactis_PT33_Early_Replicate_1_PT33/Monoculture_Bi_animalis_lactis_PT33_Early_Replicate_1_PT33.tsv"), header = TRUE, sep = "\t")
s8  <- read.table(file.path(dir_pt33, "Monoculture_Bi_animalis_lactis_PT33_Early_Replicate_2_PT33/Monoculture_Bi_animalis_lactis_PT33_Early_Replicate_2_PT33.tsv"), header = TRUE, sep = "\t")
s9  <- read.table(file.path(dir_pt33, "Monoculture_Bi_animalis_lactis_PT33_Early_Replicate_3_PT33/Monoculture_Bi_animalis_lactis_PT33_Early_Replicate_3_PT33.tsv"), header = TRUE, sep = "\t")
s10 <- read.table(file.path(dir_pt33, "Monoculture_Bi_animalis_lactis_PT33_Late_Replicate_1_PT33/Monoculture_Bi_animalis_lactis_PT33_Late_Replicate_1_PT33.tsv"),   header = TRUE, sep = "\t")
s11 <- read.table(file.path(dir_pt33, "Monoculture_Bi_animalis_lactis_PT33_Late_Replicate_2_PT33/Monoculture_Bi_animalis_lactis_PT33_Late_Replicate_2_PT33.tsv"),   header = TRUE, sep = "\t")
s12 <- read.table(file.path(dir_pt33, "Monoculture_Bi_animalis_lactis_PT33_Late_Replicate_3_PT33/Monoculture_Bi_animalis_lactis_PT33_Late_Replicate_3_PT33.tsv"),   header = TRUE, sep = "\t")


# ---------------------------------------------------------------------
# 2. Construir la matriz TPM
# ---------------------------------------------------------------------
tpm_pt33 <- data.frame(
  Consortium_Early_R1 = s1$tpm,
  Consortium_Early_R2 = s2$tpm,
  Consortium_Early_R3 = s3$tpm,
  Consortium_Late_R1  = s4$tpm,
  Consortium_Late_R2  = s5$tpm,
  Consortium_Late_R3  = s6$tpm,
  PT33_Early_R1       = s7$tpm,
  PT33_Early_R2       = s8$tpm,
  PT33_Early_R3       = s9$tpm,
  PT33_Late_R1        = s10$tpm,
  PT33_Late_R2        = s11$tpm,
  PT33_Late_R3        = s12$tpm
)
rownames(tpm_pt33) <- s1[[1]]   # primera columna = target_id (el nombre real puede ser X.e.target_id por un bug del header)


# ---------------------------------------------------------------------
# 3. Filtro + 4. log2(TPM + 1)
# ---------------------------------------------------------------------
tpm_pt33  <- tpm_pt33[rowSums(tpm_pt33 > 1) >= 2, ]
log2_pt33 <- as.matrix(log2(tpm_pt33 + 1))


# ---------------------------------------------------------------------
# 5. Leer anotaciones funcionales del genoma de PT33
# ---------------------------------------------------------------------
annot_pt33 <- read_csv(
  file.path(BASE, "deseq2_alex_kallisto/PT33/Genome_annotation/Bi_animalis_lactis_PT33_annotations.csv"),
  show_col_types = FALSE
)


# ---------------------------------------------------------------------
# 6. Leer lista de genes de interes (.fasta con categorias en ###)
# ---------------------------------------------------------------------
lineas_pt33 <- readLines(
  file.path(BASE, "deseq2_alex_kallisto/PT33/Genome_annotation/lista_genes_interes_RNAseq_paper_PT33.fasta")
)

categoria_actual <- NA
genes_id_pt  <- c()
genes_cat_pt <- c()

for (linea in lineas_pt33) {
  if (startsWith(linea, "###")) {
    categoria_actual <- trimws(sub("^#+", "", linea))
  } else if (startsWith(linea, "RVCYZ") && !is.na(categoria_actual)) {
    gene_id <- strsplit(linea, "\\s+")[[1]][1]
    genes_id_pt  <- c(genes_id_pt, gene_id)
    genes_cat_pt <- c(genes_cat_pt, categoria_actual)
  }
}

genes_interes_pt33 <- tibble(gene_id = genes_id_pt, category = genes_cat_pt)
categorias_pt33    <- unique(genes_interes_pt33$category)


# ---------------------------------------------------------------------
# 7. Seleccionar genes por categoria (todo escrito explicitamente)
# ---------------------------------------------------------------------
muestras_pt33 <- colnames(log2_pt33)
condiciones_map_pt33 <- tibble(
  Replicate = muestras_pt33,
  Condition = sub("_R[0-9]+$", "", muestras_pt33)
)
condiciones_unicas_pt33 <- unique(condiciones_map_pt33$Condition)


genes_por_cat_pt33 <- list()


# =====================================================================
# CATEGORIA ESPECIAL: Transport protein genes
# Regla: top 10 z-score promedio por cada condicion experimental
# =====================================================================
transport_disponibles <- genes_interes_pt33 %>%
  filter(category == "Transport protein genes",
         gene_id %in% rownames(log2_pt33)) %>%
  pull(gene_id)

mat_z_transport <- t(scale(t(log2_pt33[transport_disponibles, , drop = FALSE])))

genes_por_cat_pt33[["Transport protein genes"]] <- as_tibble(mat_z_transport, rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "Replicate", values_to = "Z") %>%
  left_join(condiciones_map_pt33, by = "Replicate") %>%
  group_by(gene_id, Condition) %>%
  summarise(avg_z = mean(Z, na.rm = TRUE), .groups = "drop") %>%
  group_by(Condition) %>%
  slice_max(avg_z, n = 10) %>%
  pull(gene_id) %>%
  unique()


# =====================================================================
# RESTO DE CATEGORIAS: usar todos los genes disponibles (sin filtros)
# =====================================================================
categorias_normales_pt33 <- setdiff(categorias_pt33, "Transport protein genes")

for (cat in categorias_normales_pt33) {
  genes_cat <- genes_interes_pt33 %>%
    filter(category == cat, gene_id %in% rownames(log2_pt33)) %>%
    pull(gene_id)

  if (length(genes_cat) > 0) {
    genes_por_cat_pt33[[cat]] <- genes_cat
  }
}


# ---------------------------------------------------------------------
# 8. Z-score por gen
# ---------------------------------------------------------------------
genes_todos_pt33 <- unique(unlist(genes_por_cat_pt33))
mat_pt33         <- log2_pt33[genes_todos_pt33, , drop = FALSE]
mat_z_pt33       <- t(scale(t(mat_pt33)))

# Quitar filas con NA, NaN o Inf (genes con varianza 0 -> scale() divide entre 0)
filas_validas    <- apply(mat_z_pt33, 1, function(x) all(is.finite(x)))
mat_z_pt33       <- mat_z_pt33[filas_validas, , drop = FALSE]
genes_todos_pt33 <- genes_todos_pt33[filas_validas]

cat("PT33 - genes despues del z-score:", nrow(mat_z_pt33), "\n")
cat("PT33 - hay valores no finitos?:", any(!is.finite(mat_z_pt33)), "\n")


# ---------------------------------------------------------------------
# 9. Etiquetas (gene_id - funcion) y anotaciones para el heatmap
# ---------------------------------------------------------------------

# Cascada: GenBank_Function -> ERGO_Function -> Description_eggnog -> id
# Ojo: en PT33 las columnas usan _ (no .) y eggnog se llama Description_eggnog
mapear_anotacion_pt33 <- function(gene_ids) {
  sapply(gene_ids, function(gid) {
    fila <- annot_pt33[annot_pt33$id == gid, ]
    if (nrow(fila) == 0) return(gid)

    anot <- fila$GenBank_Function[1]
    if (is.na(anot) || anot == "" || anot == "NA") anot <- fila$ERGO_Function[1]
    if (is.na(anot) || anot == "" || anot == "NA") anot <- fila$Description_eggnog[1]
    if (is.na(anot) || anot == "" || anot == "NA") return(gid)

    return(paste0(gid, " - ", anot))
  })
}

etiquetas_pt33 <- character(length(genes_todos_pt33))
for (i in seq_along(genes_todos_pt33)) {
  etiquetas_pt33[i] <- as.character(mapear_anotacion_pt33(genes_todos_pt33[i]))
}
etiquetas_pt33 <- make.unique(etiquetas_pt33, sep = "_")
rownames(mat_z_pt33) <- etiquetas_pt33


# Anotacion de filas: pathway por gen
row_annot_pt33 <- data.frame(
  Pathway = character(length(genes_todos_pt33)),
  row.names = etiquetas_pt33,
  stringsAsFactors = FALSE
)

for (cat in names(genes_por_cat_pt33)) {
  genes_de_cat  <- genes_por_cat_pt33[[cat]]
  etiquetas_cat <- mapear_anotacion_pt33(genes_de_cat)
  etiquetas_cat <- make.unique(etiquetas_cat, sep = "_")

  filas_match <- rownames(row_annot_pt33) %in% etiquetas_cat
  row_annot_pt33$Pathway[filas_match] <- cat
}


# Anotacion de columnas
col_annot_pt33 <- data.frame(
  Condition = factor(condiciones_map_pt33$Condition),
  row.names = condiciones_map_pt33$Replicate
)


# ---------------------------------------------------------------------
# 10. Colores
# ---------------------------------------------------------------------
colores_condicion_pt33 <- c(
  "Consortium_Early" = "#E41A1C",
  "Consortium_Late"  = "#377EB8",
  "PT33_Early"       = "#4DAF4A",
  "PT33_Late"        = "#984EA3"
)

n_cat_pt33 <- length(genes_por_cat_pt33)
if (n_cat_pt33 <= 12) {
  colores_pathway_pt33 <- brewer.pal(n_cat_pt33, "Set3")
} else {
  colores_pathway_pt33 <- colorRampPalette(brewer.pal(12, "Set3"))(n_cat_pt33)
}
names(colores_pathway_pt33) <- names(genes_por_cat_pt33)

ann_colors_pt33 <- list(
  Condition = colores_condicion_pt33,
  Pathway   = colores_pathway_pt33
)


# ---------------------------------------------------------------------
# 11. Generar heatmap de PT33
# ---------------------------------------------------------------------
p_complete_PT33 <- pheatmap(
  mat               = mat_z_pt33,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  annotation_col    = col_annot_pt33,
  annotation_row    = row_annot_pt33,
  annotation_colors = ann_colors_pt33,
  color             = colorRampPalette(c("green", "black", "red"))(100),
  fontsize          = 10,
  fontsize_row      = 6,
  fontsize_col      = 10,
  angle_col         = "90",
  border_color      = NA,
  scale             = "none",
  cellwidth         = 11,
  cellheight        = 7
)




# =====================================================================
# ===========  COMBINAR LOS DOS HEATMAPS CON PATCHWORK  ===============
# =====================================================================

setwd(OUT)

# Convertir cada pheatmap a ggplot para poder combinarlos
gg_HGF2 <- as.ggplot(p_complete_HGF2) +
  annotate("text", x = 0.02, y = 0.98, label = "A)",
           size = 35, fontface = "bold", hjust = 0, vjust = 1)

gg_PT33 <- as.ggplot(p_complete_PT33) +
  annotate("text", x = 0.02, y = 0.98, label = "B)",
           size = 35, fontface = "bold", hjust = 0, vjust = 1)


# ---------------------------------------------------------------------
# Layout 1x2 (lado a lado)
# ---------------------------------------------------------------------
combined_1x2 <- gg_HGF2 | gg_PT33

ggsave("HGF2_PT33_heatmaps_TPM_basic_1x2.png",
       plot = combined_1x2, width = 26, height = 13, units = "in", dpi = 600)

ggsave("HGF2_PT33_heatmaps_TPM_basic_1x2.pdf",
       plot = combined_1x2, width = 26, height = 13, units = "in")

ggsave("HGF2_PT33_heatmaps_TPM_basic_1x2.tiff",
       plot = combined_1x2, width = 26, height = 13, units = "in",
       dpi = 600, compression = "lzw")


# ---------------------------------------------------------------------
# Layout 2x1 (apilados)
# ---------------------------------------------------------------------
combined_2x1 <- gg_HGF2 / gg_PT33

ggsave("HGF2_PT33_heatmaps_TPM_basic_2x1.png",
       plot = combined_2x1, width = 14, height = 22, units = "in", dpi = 600)

ggsave("HGF2_PT33_heatmaps_TPM_basic_2x1.pdf",
       plot = combined_2x1, width = 14, height = 22, units = "in")

ggsave("HGF2_PT33_heatmaps_TPM_basic_2x1.tiff",
       plot = combined_2x1, width = 14, height = 22, units = "in",
       dpi = 600, compression = "lzw")

