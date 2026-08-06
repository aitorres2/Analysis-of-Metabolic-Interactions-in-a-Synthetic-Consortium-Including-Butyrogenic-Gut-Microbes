
# Change this next path directory to your own computer:
setwd("/home/alexis/Desktop/actualizacion_github_paper_HGF2_PT33/0_github_corregido/Enzmyes_in_genomes")

# Install libraries
#instal.packages("tidyverse")
#instal.packages("pheatmap")

# Load libraries
library(tidyverse)
library(pheatmap)

# The 8 genomes from this study and their abbreviations used here
genomas <- c(
  "VPI"     = "B. thetaiotaomicron VPI-5482",
  "PT33"    = "B. animalis lactis PT33",
  "M38"     = "L. paracasei M38",
  "7243FAA" = "Clostridium sp. 7_2_43FAA",
  "HGF2"    = "Clostridium sp. HGF2",
  "L250"    = "Clostridium sp. L2-50",
  "M62"     = "Clostridium sp. M62/1",
  "WAL"     = "C. symbiosum WAL-14673"
)

# The 17 aminoacids characterized by GapMind, 3 letters abbreviation are needed
aminoacidos <- c(
  "phe" = "Phenylalanine", "tyr" = "Tyrosine",   "trp" = "Tryptophan",
  "ile" = "Isoleucine",    "leu" = "Leucine",    "val" = "Valine",
  "ser" = "Serine",        "gly" = "Glycine",    "cys" = "Cysteine",
  "thr" = "Threonine",     "met" = "Methionine", "lys" = "Lysine",
  "asn" = "Asparagine",    "arg" = "Arginine",   "pro" = "Proline",
  "gln" = "Glutamine",     "his" = "Histidine"
)


# Read HTML from GapMind output

# =====================================================================
#  2. LEER LOS HTML DE GAPMIND
#
#     GapMind entrega una pagina web por genoma. Cada fila de su tabla es
#     un aminoacido, y el COLOR de cada enlace codifica la confianza:
#
#         verde  #007000  ->  encontrado con confianza alta
#         negro  #000000  ->  confianza media
#         rojo   #CC4444  ->  no encontrado (hueco en la ruta)
#
#     Hay dos tipos de enlace y la diferencia es importante:
#       - enlace con "&path=" solamente  -> es la RUTA entera
#       - enlace con "&path=" y "&step=" -> es UN PASO de la ruta
#
#     Para el panel A solo necesitamos el color de la RUTA, porque esa es
#     la conclusion de GapMind: "esta bacteria puede o no fabricarlo".
# =====================================================================

leer_gapmind <- function(archivo, nombre_genoma) {

  # Leer el HTML entero como un solo texto largo
  texto <- paste(readLines(archivo, warn = FALSE), collapse = "")

  # Buscar los enlaces de RUTA: llevan &path= pero NO llevan &step=
  # El patron captura dos cosas: (1) el color y (2) el nombre de la ruta
  patron_ruta <- 'style="color: (#[0-9A-Fa-f]{6})[^"]*"[^>]*>([a-z]{3})</a>'
  encontrados <- str_match_all(texto, patron_ruta)[[1]]

  tibble(
    genoma  = nombre_genoma,
    color   = encontrados[, 2],
    pathway = encontrados[, 3]
  )
}

# Aplicar la funcion a los 8 archivos y juntar todo en una sola tabla
rutas <- map_dfr(names(genomas), function(g) {
  leer_gapmind(paste0(g, ".html"), g)
})

# Quedarse solo con los 17 aminoacidos de interes.
# GapMind tambien reporta "chorismate", que es un intermediario y no un
# aminoacido, asi que se descarta aqui.
rutas <- rutas %>% filter(pathway %in% names(aminoacidos))

cat("Leidos", nrow(rutas), "resultados de GapMind\n")


# =====================================================================
#  3. PANEL A: MATRIZ DE CAPACIDAD
#
#     Convertir el color en un numero, porque pheatmap pinta numeros:
#         1 = rojo   (no puede producirlo)
#         2 = medio
#         3 = verde  (puede producirlo)
# =====================================================================

rutas <- rutas %>%
  mutate(valor = case_when(
    color == "#007000" ~ 3,   # verde
    color == "#000000" ~ 2,   # negro
    color == "#CC4444" ~ 1    # rojo
  ))

# pivot_wider convierte la tabla larga (una fila por genoma+aminoacido)
# en una tabla ancha (genomas en filas, aminoacidos en columnas).
tabla_A <- rutas %>%
  select(genoma, pathway, valor) %>%
  pivot_wider(names_from = pathway, values_from = valor)

# pheatmap necesita una MATRIZ con nombres de fila, no un data frame
matriz_A <- tabla_A %>%
  column_to_rownames("genoma") %>%
  as.matrix()

# Reordenar filas y columnas segun definimos arriba, y poner los nombres
# largos para que la figura sea legible
matriz_A <- matriz_A[names(genomas), names(aminoacidos)]
rownames(matriz_A) <- genomas
colnames(matriz_A) <- aminoacidos


# =====================================================================
#  4. DIBUJAR EL PANEL A
#
#     Sobre el tamano: se calcula a partir del numero de celdas en vez de
#     escribirlo a mano, para que la figura no se deforme si algun dia se
#     anaden mas genomas o mas aminoacidos.
# =====================================================================

ancho_A <- ncol(matriz_A) * 0.65 + 5    # 0.65 pulgadas por columna + margen
alto_A  <- nrow(matriz_A) * 0.42 + 2.5  # 0.42 pulgadas por fila + margen

figura_A <- pheatmap(
  matriz_A,
  color         = c("#D6604D", "#FDDBA6", "#1A9850"),   # rojo, crema, verde
  breaks        = c(0.5, 1.5, 2.5, 3.5),                # 3 tramos: 1, 2 y 3
  legend_breaks = c(1, 2, 3),
  legend_labels = c("gap (cannot produce)", "medium confidence",
                    "all steps found"),
  cluster_rows  = FALSE,     # mantener NUESTRO orden, no el del algoritmo
  cluster_cols  = FALSE,
  cellwidth     = 46,
  cellheight    = 30,
  angle_col     = 45,
  border_color  = "white",
  main          = ""         # sin titulo: los titulos van en el texto del paper
)

pdf("Fig3A_aa_capacity.pdf", width = ancho_A, height = alto_A)
grid::grid.draw(figura_A$gtable)
dev.off()

png("Fig3A_aa_capacity.png", width = ancho_A, height = alto_A,
    units = "in", res = 200)
grid::grid.draw(figura_A$gtable)
dev.off()

cat("Panel A guardado:", ncol(matriz_A), "aminoacidos x",
    nrow(matriz_A), "genomas\n")


# =====================================================================
#  5. PANEL B: CONTAR PROTEINAS DE LA TABLA CURADA
#
#     La tabla curada tiene una fila por proteina y una columna por
#     genoma. En cada celda hay una lista de locus separados por ";",
#     o un "0" si esa bacteria no tiene ninguno. Por ejemplo:
#
#        protein                 | M38
#        Isoleucine transporter  | pgaptmp_002848; pgaptmp_000425; ...
#        Butyrate_kinase         | 0
#
#     El panel B cuenta cuantos locus hay en cada celda.
# =====================================================================

tabla_curada <- read_tsv("tabla_heatmap_prots_crossfeeding.tsv",
                         col_types = cols(.default = col_character()))

# Los nombres de columna de la tabla curada no son los mismos que usamos
# nosotros, asi que hace falta una equivalencia.
equivalencia <- c("VPI-5482" = "VPI", "PT33" = "PT33", "M38" = "M38",
                  "C7243FAA" = "7243FAA", "HGF2" = "HGF2", "L250" = "L250",
                  "M62/1" = "M62", "WAL14673" = "WAL")

conteos <- tabla_curada %>%
  # pasar de ancho a largo: una fila por proteina + genoma
  pivot_longer(all_of(names(equivalencia)),
               names_to = "columna", values_to = "locus") %>%
  # contar los locus: si la celda es "0" o esta vacia el conteo es 0;
  # si no, es el numero de ";" mas uno
  mutate(
    genoma = equivalencia[columna],
    n = if_else(is.na(locus) | locus == "0",
                0L,
                str_count(locus, ";") + 1L)
  ) %>%
  select(proteina = protein, metabolito = metabolite_produced, genoma, n)


# =====================================================================
#  6. CLASIFICAR CADA PROTEINA EN UN BLOQUE FUNCIONAL
#
#     Cuatro bloques, que en la figura salen separados por un hueco y
#     con su propio color en la pista "Function":
#
#       Transporters            el nombre dice "transporter" o "permease"
#       Fermentation acids      butirato, acetato, lactato, succinato
#       Inulin / FOS / sugars   hidrolisis de fibra
#       Amino acid              el resto
#
#     SOBRE LA DIRECCION: la pista NO dice si una enzima de aminoacidos
#     produce o consume. Se probo y no era defendible: de 17 enzimas,
#     solo 4 tienen direccion inequivoca, 10 son reversibles (la fija la
#     fisiologia, no el gen) y 3 estaban mal clasificadas. Poner un color
#     a ese juicio le daria la misma autoridad visual que un conteo
#     medido. Donde SI es inequivoca por la clase de proteina se conserva:
#     un transportador transporta, la acetato quinasa produce acetato.
# =====================================================================

es_transportador <- function(x) str_detect(x, regex("transporter|permease", TRUE))
acidos_fermentacion <- c("Butyrate", "Acetate", "Lactate", "Succinate")

conteos <- conteos %>%
  mutate(
    bloque = case_when(
      es_transportador(proteina)                       ~ "Transporters",
      metabolito %in% acidos_fermentacion              ~ "Fermentation acids",
      str_detect(metabolito, "Glucose|fructose|FOS")   ~ "Inulin / FOS / sugars",
      TRUE                                             ~ "Amino acid"
    ),
    funcion = case_when(
      bloque == "Transporters"          ~ "transport",
      bloque == "Inulin / FOS / sugars" ~ "fiber hydrolysis",
      bloque == "Fermentation acids"    ~ "SCFA / fermentation",
      TRUE                              ~ "amino acid (not GapMind)"
    )
  )

# --- Filtrar el bloque de aminoacidos --------------------------------
# Se quitan las enzimas de aminoacidos que YA salen en el panel A por
# GapMind: repetirlas aqui alarga el eje sin anadir informacion.
# Se conservan dos grupos:
#   (a) Ala / Asp / Glu, que GapMind no cubre (se forman por
#       transaminacion desde intermediarios centrales)
#   (b) enzimas de DEGRADACION reconocibles por su clase de reaccion:
#       amonio-liasas, deshidratasas y desamidasas. Rompen el aminoacido
#       de forma irreversible, asi que indican capacidad de CONSUMIRLO
#       -- el lado receptor del cross-feeding, que GapMind no anota.
sin_gapmind <- "Alanine|Aspartate|Glutamate"
degradacion <- regex("ammonia-lyase|dehydratase|Glutaminase|Asparaginase", TRUE)

conteos <- conteos %>%
  filter(bloque != "Amino acid" |
         str_detect(metabolito, sin_gapmind) |
         str_detect(proteina, degradacion)) %>%
  # 'Threonine dehydratase' pasa el filtro, pero el nombre no distingue
  # ilvA (biosintesis de Ile) de tdcB (catabolica). En vez de excluirla o
  # de afirmar una direccion, se marca como ambigua en la etiqueta.
  mutate(proteina = if_else(proteina == "Threonine dehydratase",
                            "Threonine dehydratase (ilvA/tdcB)", proteina))

# Una tabla con una fila por proteina: es la metadata de las columnas.
# El orden de esta tabla ES el orden de las columnas de la figura.
orden_bloques <- c("Fermentation acids", "Inulin / FOS / sugars",
                   "Amino acid", "Transporters")

columnas <- conteos %>%
  distinct(proteina, bloque, funcion, metabolito) %>%
  group_by(proteina) %>%
  slice(1) %>%          # si una proteina aparece con dos metabolitos, uno
  ungroup() %>%
  mutate(bloque = factor(bloque, levels = orden_bloques)) %>%
  arrange(bloque, metabolito, proteina)

cat("Panel B:", nrow(columnas), "proteinas tras el filtro\n")

# Matriz de conteos, en el orden que acabamos de fijar
matriz_B <- conteos %>%
  distinct(proteina, genoma, n) %>%
  complete(proteina = columnas$proteina, genoma = names(genomas),
           fill = list(n = 0)) %>%
  pivot_wider(names_from = proteina, values_from = n) %>%
  column_to_rownames("genoma") %>%
  as.matrix()

matriz_B <- matriz_B[names(genomas), columnas$proteina]
rownames(matriz_B) <- genomas


# Using microshaedd for many labels "metadata"

# remotes::install_github("KarstensLab/microshades", dependencies = TRUE)
library(microshades)

familia_de <- function(x) case_when(
  str_detect(x, "Butyrate|Acetate|Lactate|Succinate")               ~ "Fermentation",
  str_detect(x, "Glucose|fructose|FOS|Inulin")                      ~ "Fiber & sugars",
  str_detect(x, "Phenylalanine|Tyrosine|Tryptophan")                ~ "Aromatic",
  str_detect(x, "Isoleucine|Leucine|Valine|Alanine")                ~ "BCAA",
  str_detect(x, "Serine|Glycine|Cysteine")                          ~ "Serine fam.",
  str_detect(x, "Threonine|Methionine|Lysine|Asparagine|Aspartate") ~ "Aspartate fam.",
  str_detect(x, "Arginine|Proline|Glutam")                          ~ "Glutamate fam.",
  str_detect(x, "Histidine")                                        ~ "Histidine",
  TRUE                                                              ~ "Other")

tonos <- c("Fermentation"   = "micro_cvd_gray",
           "Fiber & sugars" = "micro_cvd_turquoise",
           "Aromatic"       = "micro_cvd_orange",
           "BCAA"           = "micro_cvd_green",
           "Serine fam."    = "micro_cvd_purple",
           "Aspartate fam." = "micro_cvd_blue",
           "Glutamate fam." = "micro_cvd_turquoise",
           "Histidine"      = "micro_cvd_gray",
           "Other"          = "micro_cvd_gray")

# Para cada familia: coger su paleta de 5 tonos e interpolar tantas
# sombras como metabolitos tenga esa familia.
paleta_metabolitos <- function(niveles) {
  familias <- familia_de(niveles)
  colores <- unlist(lapply(names(tonos), function(f) {
    de_esta_familia <- niveles[familias == f]
    if (length(de_esta_familia) == 0) return(NULL)
    base <- as.character(microshades_palette(unname(tonos[f]), 5))
    setNames(colorRampPalette(base)(length(de_esta_familia)), de_esta_familia)
  }))
  colores <- colores[!duplicated(names(colores))]
  colores[niveles]
}


# =====================================================================
#  8. DIBUJAR EL PANEL B
#
#     Los conteos van de 0 a ~30, pero la mayoria son numeros pequenos.
#     Con una escala continua casi todo saldria del mismo color, asi que
#     se agrupan en tramos: 0, 1-2, 3-7, 8-15, 16+.
# =====================================================================

niveles_metabolito <- unique(columnas$metabolito)

# Las dos pistas de color que van encima de las columnas
anotacion_B <- data.frame(
  Metabolite = factor(columnas$metabolito, levels = niveles_metabolito),
  Function   = factor(columnas$funcion,
                      levels = c("SCFA / fermentation", "fiber hydrolysis",
                                 "amino acid (not GapMind)", "transport")),
  row.names  = columnas$proteina
)

colores_funcion <- c("SCFA / fermentation"      = "#E68310",
                     "fiber hydrolysis"         = "#2F9E44",
                     "amino acid (not GapMind)" = "#7A8B99",
                     "transport"                = "#7048E8")

# Etiqueta de columna: proteina + metabolito entre corchetes. Con ~28
# metabolitos el sombreado no siempre alcanza para identificarlos, asi
# que el corchete da una segunda via de lectura.
etiquetas_B <- paste0(columnas$proteina, "  [", columnas$metabolito, "]")

# Hueco entre bloques: donde cambia el bloque, se separa visualmente
huecos <- which(diff(as.integer(columnas$bloque)) != 0)

maximo <- max(matriz_B)

figura_B <- pheatmap(
  matriz_B,
  color  = c("white", "#D8F1B2", "#99C5E3", "#F5A040", "#EB4F48", "#8B1A1A"),
  breaks = c(-0.5, 0.5, 1.5, 3.5, 7.5, 15.5, maximo + 0.5),
  legend_breaks     = c(0, 1, 3, 5, 8, 16, maximo),
  labels_col        = etiquetas_B,
  annotation_col    = anotacion_B,
  annotation_colors = list(Metabolite = paleta_metabolitos(niveles_metabolito),
                           Function   = colores_funcion),
  gaps_col      = huecos,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  cellwidth     = 11,
  cellheight    = 10,
  fontsize      = 8,
  fontsize_row  = 8,
  angle_col     = 90,
  border_color  = NA,
  main          = ""
)

ancho_B <- ncol(matriz_B) * 14 / 72 + 4.5 + 1.5 + 3.5
alto_B  <- max(nrow(matriz_B) * 10 / 72 + max(nchar(etiquetas_B)) * 6.5 / 72 + 2,
               (length(niveles_metabolito) + 4 + 6 + 6) * 14 / 72 + 1.5)

pdf("Fig3B_scfa_fiber_transporters.pdf", width = ancho_B, height = alto_B)
grid::grid.draw(figura_B$gtable)
dev.off()

png("Fig3B_scfa_fiber_transporters.png", width = ancho_B, height = alto_B,
    units = "in", res = 200)
grid::grid.draw(figura_B$gtable)
dev.off()

cat("Panel B guardado:", ncol(matriz_B), "proteinas x",
    nrow(matriz_B), "genomas\n")


# Save data
write_tsv(as.data.frame(matriz_A) %>% rownames_to_column("bacteria"),
          "Fig3A_capacidad_aminoacidos.tsv")

write_tsv(as.data.frame(matriz_B) %>% rownames_to_column("bacteria"),
          "Fig3B_matriz_conteos.tsv")
