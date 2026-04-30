# =====================================================================
# plot_flux_sampling_Late_AcConstrained_basic_simple.R
# Version basica del script de visualizacion del flux sampling
# Solo procesa el escenario Late + Inulina (AcConstrained, basic)
# =====================================================================
#
# Pasos:
#  1. Leer rxn_names CSV  -> nombres de las reacciones
#  2. Leer samples CSV    -> matriz de flujos (rxns x puntos)
#  3. Filtrar solo reacciones IEX (transferencias entre especies)
#  4. Parsear bacteria (HGF2/PT33) y metabolito (BiGG id)
#  5. Construir data frame en formato largo para ggplot
#  6. Generar histogramas por categoria (5 categorias)
#  7. Generar matriz de correlacion Spearman entre HGF2 y PT33
# =====================================================================

library(tidyverse)
library(data.table)
library(scales)


# ---------------------------------------------------------------------
# Directorio de trabajo y archivos
# ---------------------------------------------------------------------
# setwd("/media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme/1_contextualizacion_trasncriptomica/flux_sampling")

# file_samples <- "samples_Consortium_Late_Inulina_AcConstrained_basic_200k.csv"
# file_rxns    <- "rxn_names_Consortium_Late_Inulina_AcConstrained_basic.csv"

setwd("/media/alexis/hdd2/postdoc/0_Simulaciones_MATLAB_paper_HGF2_PT33_correcciones/00_nuevos_modelos_metabolicos_carveme_gapseq/1_output_modelos_carveme/1_contextualizacion_trasncriptomica/flux_sampling")

file_samples <- "samples_Consortium_Late_Inulina_AcConstrained_sin_proteccion_200k.csv"                                               
file_rxns    <- "rxn_names_Consortium_Late_Inulina_AcConstrained_sin_proteccion.csv"

# Sufijo para los archivos de salida
# sufijo_out <- "Late_Inulina_AcConstr_basic"
sufijo_out <- "Late_Inulina_AcConstr_sin_proteccion" 

# ---------------------------------------------------------------------
# Mapeo de IDs BiGG a nombres legibles
# ---------------------------------------------------------------------
bigg_to_name <- c(
  "but"      = "Butyrate",
  "ac"       = "Acetate",
  "lac__L"   = "L-Lactate",
  "lac__D"   = "D-Lactate",
  "succ"     = "Succinate",
  "for"      = "Formate",
  "etoh"     = "Ethanol",
  "co2"      = "CO2",
  "fru"      = "D-Fructose",
  "glc__D"   = "D-Glucose",
  "inulin"   = "Inulin",
  "kesto"    = "Kestose",
  "kestopt"  = "Kestopentaose",
  "kestottr" = "Kestotetraose",
  "glu__L"   = "L-Glutamate",
  "gln__L"   = "L-Glutamine",
  "ala__L"   = "L-Alanine",
  "asp__L"   = "L-Aspartate",
  "asn__L"   = "L-Asparagine",
  "val__L"   = "L-Valine",
  "phe__L"   = "L-Phenylalanine",
  "ile__L"   = "L-Isoleucine",
  "leu__L"   = "L-Leucine",
  "pro__L"   = "L-Proline",
  "ser__L"   = "L-Serine",
  "thr__L"   = "L-Threonine",
  "his__L"   = "L-Histidine",
  "gly"      = "Glycine",
  "lys__L"   = "L-Lysine",
  "trp__L"   = "L-Tryptophan",
  "tyr__L"   = "L-Tyrosine",
  "met__L"   = "L-Methionine",
  "arg__L"   = "L-Arginine",
  "cys__L"   = "L-Cysteine",
  "gthrd"    = "Glutathione",
  "cgly"     = "Cys-Gly",
  "inost"    = "Inositol",
  "ade"      = "Adenine",
  "hxan"     = "Hypoxanthine",
  "h2s"      = "H2S",
  "nh4"      = "Ammonium",
  "acald"    = "Acetaldehyde",
  "mal__L"   = "L-Malate",
  "ala__D"   = "D-Alanine",
  "glyc"     = "Glycerol",
  "hom__L"   = "L-Homoserine",
  "fum"      = "Fumarate",
  "rib__D"   = "D-Ribose"
)


# ---------------------------------------------------------------------
# Categorias de metabolitos para histogramas
# ---------------------------------------------------------------------
cat_fermentacion <- c("Butyrate", "Acetate", "L-Lactate", "D-Lactate",
                      "Succinate", "Formate", "Ethanol", "CO2")

cat_azucares <- c("D-Fructose", "D-Glucose", "Inulin",
                  "Kestopentaose", "Kestose", "Kestotetraose")

cat_aminoacidos_1 <- c("L-Glutamate", "L-Glutamine", "L-Alanine",
                       "L-Aspartate", "L-Asparagine", "L-Valine",
                       "L-Isoleucine", "L-Serine", "L-Threonine",
                       "Glycine", "L-Histidine", "L-Cysteine")

cat_aminoacidos_2 <- c("L-Leucine", "L-Proline", "L-Arginine",
                       "L-Lysine", "L-Phenylalanine", "L-Tyrosine",
                       "L-Tryptophan", "L-Methionine",
                       "D-Alanine", "L-Homoserine")

cat_otros <- c("Glutathione", "Cys-Gly", "Inositol", "Adenine",
               "Hypoxanthine", "H2S", "Ammonium", "Acetaldehyde",
               "L-Malate", "Glycerol", "Fumarate", "D-Ribose")

# Todos juntos (para la matriz de correlacion)
todos_metabolitos <- c(cat_fermentacion, cat_azucares,
                       cat_aminoacidos_1, cat_aminoacidos_2, cat_otros)


# =====================================================================
# 1. Leer rxn_names CSV
# =====================================================================
cat("\n[1] Leyendo rxn_names ...\n")
rxn_names <- fread(file_rxns, header = FALSE)$V1
cat("    Total de reacciones:", length(rxn_names), "\n")


# =====================================================================
# 2. Leer samples CSV  (matriz grande, fread es la opcion rapida)
# =====================================================================
cat("\n[2] Leyendo samples (puede tardar varios minutos por el tamano del archivo) ...\n")
sampling <- fread(file_samples)
cat("    Dimensiones:", nrow(sampling), "filas x", ncol(sampling), "columnas\n")

# Verificacion: nrow(sampling) debe coincidir con length(rxn_names)
if (nrow(sampling) != length(rxn_names)) {
  stop("ERROR: nrow(sampling)=", nrow(sampling),
       " no coincide con length(rxn_names)=", length(rxn_names))
}
rownames(sampling) <- rxn_names


# =====================================================================
# 3. Filtrar solo reacciones IEX
# =====================================================================
cat("\n[3] Filtrando reacciones IEX ...\n")
es_iex <- grepl("IEX_", rxn_names)
sampling_iex <- sampling[es_iex, ]
rxn_iex      <- rxn_names[es_iex]
cat("    Reacciones IEX:", length(rxn_iex), "\n")


# =====================================================================
# 4. Parsear bacteria + metabolito de cada nombre de reaccion IEX
#    Patron: HGF2IEX_xxx[u]tr  o  PT33IEX_xxx[u]tr
# =====================================================================
cat("\n[4] Parseando bacteria + metabolito ...\n")

# Vectores paralelos a rxn_iex con bacteria, met_id, met_name
bacteria_vec <- character(length(rxn_iex))
met_id_vec   <- character(length(rxn_iex))
met_name_vec <- character(length(rxn_iex))

for (i in seq_along(rxn_iex)) {
  rxn <- rxn_iex[i]
  m <- regmatches(rxn, regexec("^(HGF2|PT33)IEX_(.+)\\[u\\]tr$", rxn))[[1]]

  if (length(m) == 3) {
    bacteria_vec[i] <- m[2]
    met_id_vec[i]   <- m[3]
    # Si el met_id esta en el mapeo, usar el nombre legible; si no, usar el id
    if (m[3] %in% names(bigg_to_name)) {
      met_name_vec[i] <- bigg_to_name[[m[3]]]
    } else {
      met_name_vec[i] <- m[3]
    }
  } else {
    bacteria_vec[i] <- NA
    met_id_vec[i]   <- NA
    met_name_vec[i] <- NA
  }
}
cat("    Reacciones IEX parseadas:", sum(!is.na(bacteria_vec)), "\n")


# =====================================================================
# 5. Construir data frame en formato largo para ggplot
#    Solo metabolitos que estan en alguna categoria y que tienen flujo > 0
# =====================================================================
cat("\n[5] Construyendo data frame largo para histogramas ...\n")

all_data <- data.frame()

for (i in seq_along(rxn_iex)) {

  if (is.na(met_name_vec[i])) next
  if (!(met_name_vec[i] %in% todos_metabolitos)) next

  vals <- as.numeric(sampling_iex[i, ])

  # Saltar metabolitos con flujo practicamente cero en todos los puntos
  if (all(abs(vals) < 1e-9)) next

  df_temp <- data.frame(
    value      = vals,
    bacteria   = bacteria_vec[i],
    metabolite = met_name_vec[i]
  )
  all_data <- rbind(all_data, df_temp)
}

cat("    Filas en all_data:", nrow(all_data), "\n")


# =====================================================================
# 6. Histogramas por categoria (5 categorias)
#    Loop simple sobre una lista de (mets, titulo, archivo, ncol)
# =====================================================================
cat("\n[6] Generando histogramas por categoria ...\n")

categorias_para_plot <- list(
  list(mets = cat_fermentacion,  titulo = "Fermentation products",     archivo = "hist_fermentacion",  ncol = 4),
  list(mets = cat_azucares,      titulo = "Sugars & substrates",       archivo = "hist_azucares",      ncol = 3),
  list(mets = cat_aminoacidos_1, titulo = "Amino acids (group 1)",     archivo = "hist_aminoacidos1",  ncol = 4),
  list(mets = cat_aminoacidos_2, titulo = "Amino acids (group 2)",     archivo = "hist_aminoacidos2",  ncol = 4),
  list(mets = cat_otros,         titulo = "Other metabolites",         archivo = "hist_otros",         ncol = 4)
)

for (cat_info in categorias_para_plot) {

  # Filtrar al subset de la categoria
  datos_cat <- all_data %>%
    filter(metabolite %in% cat_info$mets) %>%
    mutate(metabolite = factor(metabolite, levels = cat_info$mets)) %>%
    filter(!is.na(metabolite))

  if (nrow(datos_cat) == 0) {
    cat("    Sin datos para", cat_info$titulo, "\n")
    next
  }

  # Calcular tamano de la figura segun cuantos metabolitos hay
  n_mets <- length(unique(datos_cat$metabolite))
  n_rows <- ceiling(n_mets / cat_info$ncol)
  alto   <- max(4, n_rows * 3.5 + 1)
  ancho  <- cat_info$ncol * 5.5 + 1

  # Plot
  p <- ggplot(datos_cat, aes(x = value, fill = bacteria)) +
    geom_histogram(position = "identity", bins = 120, alpha = 0.5) +
    labs(x = expression("Flux [mmol gDW"^{-1}~"h"^{-1}*"]"),
         y = "Frequency",
         title = paste0(cat_info$titulo, " - Late Inulina (Ac constrained)")) +
    scale_fill_manual(name   = "Species",
                      values = c("HGF2" = "#d95f02", "PT33" = "#1b9e77"),
                      labels = c("HGF2" = "Clostridium sp. HGF2",
                                 "PT33" = "Bi. animalis lactis PT33")) +
    facet_wrap(~metabolite, ncol = cat_info$ncol, scales = "free", strip.position = "top") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text       = element_text(size = 14, face = "bold"),
          panel.spacing    = unit(1, "lines"),
          axis.text        = element_text(size = 12, face = "bold"),
          axis.title       = element_text(size = 13, face = "bold"),
          legend.text      = element_text(size = 13),
          legend.title     = element_text(size = 13, face = "bold"),
          plot.title       = element_text(size = 16, face = "bold"),
          text             = element_text(size = 12))

  # Guardar PNG y PDF
  nombre_archivo <- paste0(cat_info$archivo, "_", sufijo_out)
  ggsave(paste0(nombre_archivo, ".png"), plot = p,
         width = ancho, height = alto, units = "in", dpi = 600, bg = "white")
  ggsave(paste0(nombre_archivo, ".pdf"), plot = p,
         width = ancho, height = alto, units = "in", dpi = 600, bg = "white")

  cat("    Guardado:", nombre_archivo, "\n")
}


# =====================================================================
# 7. Matriz de correlacion Spearman entre HGF2 y PT33 (cross-species)
# =====================================================================
cat("\n[7] Generando matriz de correlacion Spearman ...\n")

# Transponer la matriz: filas = puntos del sampling, columnas = reacciones
sampling_t <- as.data.frame(t(sampling_iex))
colnames(sampling_t) <- rxn_iex

# Separar columnas HGF2 y PT33
cols_hgf2 <- grep("^HGF2IEX_", colnames(sampling_t), value = TRUE)
cols_pt33 <- grep("^PT33IEX_", colnames(sampling_t), value = TRUE)

# Eliminar reacciones con sd = 0 (cor() falla con varianza cero)
sd_hgf2 <- apply(sampling_t[, cols_hgf2, drop = FALSE], 2, sd)
sd_pt33 <- apply(sampling_t[, cols_pt33, drop = FALSE], 2, sd)
cols_hgf2 <- cols_hgf2[sd_hgf2 > 1e-12]
cols_pt33 <- cols_pt33[sd_pt33 > 1e-12]

# Vector con el nombre legible de cada reaccion (paralelo a las columnas)
nombres_hgf2 <- character(length(cols_hgf2))
for (i in seq_along(cols_hgf2)) {
  m <- regmatches(cols_hgf2[i], regexec("IEX_(.+)\\[u\\]tr$", cols_hgf2[i]))[[1]]
  if (length(m) == 2 && m[2] %in% names(bigg_to_name)) {
    nombres_hgf2[i] <- bigg_to_name[[m[2]]]
  } else {
    nombres_hgf2[i] <- NA
  }
}

nombres_pt33 <- character(length(cols_pt33))
for (i in seq_along(cols_pt33)) {
  m <- regmatches(cols_pt33[i], regexec("IEX_(.+)\\[u\\]tr$", cols_pt33[i]))[[1]]
  if (length(m) == 2 && m[2] %in% names(bigg_to_name)) {
    nombres_pt33[i] <- bigg_to_name[[m[2]]]
  } else {
    nombres_pt33[i] <- NA
  }
}

# Filtrar solo metabolitos que estan en la lista de interes
keep_hgf2 <- nombres_hgf2 %in% todos_metabolitos
keep_pt33 <- nombres_pt33 %in% todos_metabolitos

cols_hgf2_filt    <- cols_hgf2[keep_hgf2]
cols_pt33_filt    <- cols_pt33[keep_pt33]
nombres_hgf2_filt <- nombres_hgf2[keep_hgf2]
nombres_pt33_filt <- nombres_pt33[keep_pt33]

if (length(cols_hgf2_filt) > 0 && length(cols_pt33_filt) > 0) {

  # Calcular matriz de correlacion Spearman
  mat_corr <- cor(sampling_t[, cols_hgf2_filt, drop = FALSE],
                  sampling_t[, cols_pt33_filt, drop = FALSE],
                  method = "spearman")

  rownames(mat_corr) <- paste0("HGF2_", nombres_hgf2_filt)
  colnames(mat_corr) <- paste0("PT33_", nombres_pt33_filt)

  # Pasar a formato largo con tidyr (en vez de reshape2::melt)
  mat_long <- as.data.frame(mat_corr) %>%
    rownames_to_column(var = "HGF2") %>%
    pivot_longer(cols = -HGF2, names_to = "PT33", values_to = "Correlation")

  n_h     <- length(unique(mat_long$HGF2))
  n_p     <- length(unique(mat_long$PT33))
  ancho_c <- max(10, n_p * 0.55 + 4)
  alto_c  <- max(8,  n_h * 0.55 + 3)

  # Plot
  p_corr <- ggplot(mat_long, aes(x = PT33, y = HGF2, fill = Correlation)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_gradientn(
      colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
      values  = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),
      limits  = c(-1, 1),
      name    = "Spearman\nCorrelation"
    ) +
    labs(title = "Cross-species Spearman correlation - Late Inulina (Ac constrained)") +
    theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", colour = "white"),
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.title       = element_text(size = 14, face = "bold"),
      axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1,
                                      size = 9, color = "black", face = "bold"),
      axis.text.y      = element_text(size = 9, color = "black", face = "bold"),
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      legend.title     = element_text(size = 12, face = "bold"),
      legend.text      = element_text(size = 10, face = "bold")
    ) +
    coord_fixed(ratio = 1)

  # Guardar
  nombre_corr <- paste0("correlation_spearman_", sufijo_out)
  ggsave(paste0(nombre_corr, ".png"), plot = p_corr,
         width = ancho_c, height = alto_c, dpi = 600, bg = "white")
  ggsave(paste0(nombre_corr, ".pdf"), plot = p_corr,
         width = ancho_c, height = alto_c, dpi = 600, bg = "white")

  cat("    Guardado:", nombre_corr, "\n")

} else {
  cat("    No hay reacciones validas para calcular correlacion\n")
}


# =====================================================================
# 8. HISTOGRAMA DE METABOLITOS SELECCIONADOS PARA EL PAPER
#    Una sola figura combinada con los 16 metabolitos clave
# =====================================================================
cat("\n[8] Generando histograma de metabolitos seleccionados (paper) ...\n")

# Lista de metabolitos seleccionados (en el orden deseado para la figura)
mets_selected <- c("Butyrate", "Acetate", "L-Lactate", "D-Lactate",
                   "D-Fructose", "D-Glucose", "Inulin", "Formate",
                   "Kestopentaose", "Inositol", "Ammonium",
                   "L-Alanine", "L-Aspartate", "L-Valine",
                   "L-Isoleucine", "L-Cysteine")

# Filtrar all_data al subset seleccionado y reordenar como factor
datos_sel <- all_data %>%
  filter(metabolite %in% mets_selected) %>%
  mutate(metabolite = factor(metabolite, levels = mets_selected)) %>%
  filter(!is.na(metabolite))

n_mets_sel <- length(unique(datos_sel$metabolite))
cat("    Metabolitos con datos:", n_mets_sel, "\n")
cat("    Metabolitos encontrados:", paste(unique(datos_sel$metabolite), collapse = ", "), "\n")

# Tamano de la figura
ncol_facet_sel <- 4
n_rows_sel <- ceiling(n_mets_sel / ncol_facet_sel)
alto_sel  <- max(4, n_rows_sel * 3.5 + 1.5)
ancho_sel <- ncol_facet_sel * 5.5 + 1

# Plot — estilo paper: fuentes mas grandes, nombres de especies en cursiva
p_sel <- ggplot(datos_sel, aes(x = value, fill = bacteria)) +
  geom_histogram(position = "identity", bins = 120, alpha = 0.5) +
  labs(x = expression(bold("Flux [mmol gDW"^{-1}~"h"^{-1}*"]")),
       y = "Frequency") +
  scale_fill_manual(name   = "Species",
                    values = c("HGF2" = "#d95f02", "PT33" = "#1b9e77"),
                    labels = c("HGF2" = expression(italic("Clostridium")~"sp. HGF2"),
                               "PT33" = expression(italic("B. animalis")~"subsp."~italic("lactis")~"PT33"))) +
  facet_wrap(~metabolite, ncol = ncol_facet_sel, scales = "free", strip.position = "top") +
  theme_bw() +
  theme(strip.background  = element_blank(),
        strip.text        = element_text(size = 16, face = "bold"),
        panel.spacing     = unit(1, "lines"),
        axis.text         = element_text(size = 14, face = "bold"),
        axis.title        = element_text(size = 16, face = "bold"),
        legend.text       = element_text(size = 15),
        legend.title      = element_text(size = 15, face = "bold"),
        legend.text.align = 0,
        legend.key.size   = unit(0.7, "cm"),
        text              = element_text(size = 14, face = "bold"))

# Guardar PNG y PDF
nombre_sel <- paste0("hist_selected_metabolites_", sufijo_out)
ggsave(paste0(nombre_sel, ".png"), plot = p_sel,
       width = ancho_sel, height = alto_sel, units = "in", dpi = 600, bg = "white")
ggsave(paste0(nombre_sel, ".pdf"), plot = p_sel,
       width = ancho_sel, height = alto_sel, units = "in", dpi = 600, bg = "white")

cat("    Guardado:", nombre_sel, "\n")


cat("\nTerminado.\n")

