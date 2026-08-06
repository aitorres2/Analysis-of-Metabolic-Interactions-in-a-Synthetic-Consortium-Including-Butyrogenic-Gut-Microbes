# =====================================================================
#  fig7.R  --  Figure 7: gene expression (RNA-seq)
#
#  WHAT IT DOES
#    Panel A: HGF2 in monoculture, Early (8 h) vs Late (45 h)
#    Panel B: PT33, monoculture vs co-culture, both timepoints
#
#  HOW TO RUN
#    Put this file in the same folder as the data and run:
#        Rscript fig7.R
#
#  INPUT FILES (all in the same folder)
#    Raw_Counts_HGF2_only_kallisto_alex.tsv    <- raw kallisto counts
#    Raw_Counts_PT33_only_kallisto_alex.tsv
#    genes_interes_RNAseq_HGF2_PT33.tsv        <- which genes to plot
#
#  OUTPUT FILES
#    Fig7A_HGF2_monoculture_early_late.pdf / .png
#    Fig7B_PT33_complete.pdf / .png
#
#  DEPENDENCIES
#    install.packages(c("tidyverse", "pheatmap", "RColorBrewer"))
#    remotes::install_github("KarstensLab/microshades")
#
#  microshades is cited as:
#    Dahl EM, Neer E, Bowie KR, Leung ET, Karstens L. 2022. microshades:
#    an R package for improving color accessibility and organization of
#    microbiome data. Microbiol Resour Announc 11:e00795-22.
#    https://doi.org/10.1128/mra.00795-22
# =====================================================================

# No setwd() on purpose: the script reads and writes in the current working
# directory, so it runs anywhere the data folder is. Open the folder in
# RStudio, or `cd` into it before calling Rscript.

setwd("/home/alexis/Desktop/actualizacion_github_paper_HGF2_PT33/0_github_corregido/Figure_7")

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(microshades)


# =====================================================================
#  1. READ THE DATA
# =====================================================================

# The count tables carry gene IDs in the first column ("Feature").
# column_to_rownames moves them to row names, which is what pheatmap
# expects, and as.matrix turns the rest into a numeric matrix.
counts_HGF2 <- read_tsv("Raw_Counts_HGF2_only_kallisto_alex.tsv",
                        show_col_types = FALSE) %>%
  column_to_rownames("Feature") %>%
  as.matrix()

counts_PT33 <- read_tsv("Raw_Counts_PT33_only_kallisto_alex.tsv",
                        show_col_types = FALSE) %>%
  column_to_rownames("Feature") %>%
  as.matrix()

# One sample is named "HGF2_Late_R13" where it should be "_R3". Fixed here
# so every sample follows the same NAME_TIMEPOINT_REPLICATE pattern.
colnames(counts_HGF2) <- sub("_R13$", "_R3", colnames(counts_HGF2))

# The gene list says which genes to plot and how to group them.
# Columns used here: genome, locus_tag, gene, categoria, metab_lab, familia
genes <- read_tsv("genes_interes_RNAseq_HGF2_PT33.tsv", show_col_types = FALSE)

cat("HGF2 counts:", nrow(counts_HGF2), "genes x",
    ncol(counts_HGF2), "samples\n")
cat("PT33 counts:", nrow(counts_PT33), "genes x",
    ncol(counts_PT33), "samples\n")
cat("Genes of interest:", nrow(genes), "\n")


# =====================================================================
#  2. FROM RAW COUNTS TO Z-SCORE
#
#     Three steps, each solving a different problem:
#
#     (a) CPM: each sample was sequenced to a different depth. Without
#         correcting for that, a sample with twice the reads looks twice
#         as expressed in EVERY gene. CPM = counts per million, i.e. the
#         fraction of the library each gene represents.
#
#     (b) log2: expression spans several orders of magnitude (from 10 to
#         100,000). On a linear scale only the highest genes would be
#         visible. We add 1 before the log because log2(0) is -Inf.
#
#     (c) per-gene z-score: each gene is centred on its own mean and
#         divided by its own SD. The colour then means "high or low FOR
#         THIS GENE", which is the comparison of interest, instead of
#         "high or low compared to other genes".
#
#     IMPORTANT: the z-score uses ONLY the columns the panel compares.
#     Computing it over all samples at once folds variation from
#     conditions we are not comparing into the SD, which flattens the
#     contrast we actually want to see.
# =====================================================================

compute_zscore <- function(count_matrix, columns, chosen_genes) {

  # (a) CPM.
  #     WATCH THE DENOMINATOR: library size is the total over ALL genes in
  #     the genome, not just the panel genes. Dividing by the sum of the
  #     chosen genes makes "per million" mean something different in each
  #     panel and the values stop being comparable (measured: shifts the
  #     z-score by up to 3 units).
  library_size <- colSums(count_matrix[, columns])

  # Now subset to the panel rows
  sub_matrix <- count_matrix[chosen_genes, columns]

  # t() appears twice because R divides row-wise: transpose, divide by
  # each sample's library size, transpose back.
  cpm <- t(t(sub_matrix) / library_size * 1e6)

  # (b) log
  log_cpm <- log2(cpm + 1)

  # (c) per-gene z-score. scale() works column-wise and we want it
  #     per gene (row), hence the two t().
  z <- t(scale(t(log_cpm)))

  return(z)
}


# =====================================================================
#  3. COLOURS FOR THE ANNOTATION TRACKS
#
#     Each panel carries two colour tracks on the left:
#
#       Pathway                   the functional category of the gene
#       Amino acid / metabolite   the metabolite it maps to
#
#     KEY POINT: levels and colours are computed ONCE, over the union of
#     both genomes. If each panel built its own palette from its own
#     levels, the same metabolite would come out one colour in panel A
#     and another in panel B -- which rules out a single shared legend
#     and invites readers to compare colours that are not comparable.
# =====================================================================

# --- "Pathway" track -------------------------------------------------
# Set2 maxes out at 8 colours. With 10 categories brewer.pal recycles and
# two categories end up identical, which makes the track unreadable. We
# interpolate instead, guaranteeing one colour per category.
pathway_levels <- sort(unique(genes$categoria))
pathway_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(pathway_levels)),
  pathway_levels
)
stopifnot(!anyDuplicated(pathway_colors))   # a collision voids the track

# --- "Amino acid / metabolite" track ---------------------------------
# ~38 metabolites, far too many for a qualitative palette. We use a
# HIERARCHICAL palette: one hue per biological family, and within each
# family a different shade per metabolite. microshades supplies
# colour-vision-deficiency (CVD) safe palettes.
#
# Only 7 CVD hues are available, so "Fermentation" and "Fiber & sugars"
# are merged into a single non-amino-acid family. The distinction is NOT
# lost: the Pathway track still separates them.
hues <- c("Fermentation & fiber" = "micro_brown",
          "Aromatic"             = "micro_cvd_orange",
          "Pyruvate / BCAA"      = "micro_cvd_green",
          "Serine family"        = "micro_cvd_purple",
          "Aspartate family"     = "micro_cvd_blue",
          "Glutamate family"     = "micro_cvd_turquoise",
          "Unspecified"          = "micro_cvd_gray")
stopifnot(!anyDuplicated(hues))

# The gene list carries "Fermentation" and "Fiber & sugars" as separate
# families. They are merged here into the single non-amino-acid family so
# that the number of families (7) matches the number of hues (7).
genes <- genes %>%
  mutate(familia = if_else(familia %in% c("Fermentation", "Fiber & sugars"),
                           "Fermentation & fiber", familia))

# Metabolite order: grouped by category and family, so similar colours sit
# together in the legend
metabolite_levels <- genes %>%
  arrange(categoria, familia, metab_lab) %>%
  pull(metab_lab) %>%
  unique()

family_of_each <- genes$familia[match(metabolite_levels, genes$metab_lab)]

# Guard: a family with no hue leaves its metabolites colourless and
# pheatmap fails with an unhelpful message.
stopifnot(all(family_of_each %in% names(hues)))

# For each family: take its 5-shade palette and interpolate as many
# shades as that family has metabolites
metabolite_colors <- unlist(lapply(names(hues), function(f) {
  in_this_family <- metabolite_levels[family_of_each == f]
  if (length(in_this_family) == 0) return(NULL)
  base <- as.character(microshades_palette(unname(hues[f]), 5))
  setNames(colorRampPalette(base)(length(in_this_family)), in_this_family)
}))
metabolite_colors <- metabolite_colors[!duplicated(names(metabolite_colors))]
metabolite_colors <- metabolite_colors[metabolite_levels]

# --- Column track colours --------------------------------------------
timepoint_colors <- c("Early (8 h)" = "#4DAF4A", "Late (45 h)"  = "#984EA3")
culture_colors   <- c("Co-culture"  = "#D95F02", "Monoculture"  = "#1B9E77")

# --- Heatmap colour scale --------------------------------------------
# Green-black-red, clipped at +-2: without clipping a single extreme gene
# eats the whole scale and everything else looks flat.
z_scale  <- colorRampPalette(c("#00FF00", "#000000", "#FF0000"))(100)
z_breaks <- seq(-2, 2, length.out = 101)


# =====================================================================
#  4. ONE FUNCTION TO DRAW EITHER PANEL
#
#     Both panels are drawn the same way except for the column tracks,
#     so the drawing code is written once.
# =====================================================================

draw_panel <- function(z, panel_genes, column_annotation,
                       column_colors, filename) {

  # Order genes by category and family so the colour tracks come out in
  # solid blocks instead of scattered
  ordered <- panel_genes %>%
    arrange(categoria, familia, metab_lab, locus_tag)
  z <- z[ordered$locus_tag, ]

  # Row label: locus + gene name. make.unique keeps two genes with the
  # same name from clashing (pheatmap requires unique row names).
  labels <- make.unique(
    paste0(ordered$locus_tag, " - ", str_trunc(ordered$gene, 44)), sep = "_")
  rownames(z) <- labels

  # The two left-hand tracks
  row_annotation <- data.frame(
    `Amino acid / metabolite` = factor(ordered$metab_lab,
                                       levels = metabolite_levels),
    Pathway = factor(ordered$categoria, levels = pathway_levels),
    row.names = labels,
    check.names = FALSE       # keep the name with spaces and a slash
  )

  # Row-track colours plus column-track colours
  all_colors <- c(
    list(`Amino acid / metabolite` = metabolite_colors,
         Pathway = pathway_colors),
    column_colors
  )

  figure <- pheatmap(
    z,
    color             = z_scale,
    breaks            = z_breaks,
    scale             = "none",   # the z-score is already computed
    annotation_row    = row_annotation,
    annotation_col    = column_annotation,
    annotation_colors = all_colors,
    cluster_rows      = FALSE,    # keep OUR ordering by category
    cluster_cols      = TRUE,     # group similar samples
    # visually separate the category blocks
    gaps_row          = which(diff(as.integer(row_annotation$Pathway)) != 0),
    cellwidth         = 12,
    cellheight        = 8,
    fontsize          = 9,
    fontsize_row      = 6,
    fontsize_col      = 9,
    angle_col         = 90,
    border_color      = NA,
    main              = ""        # titles belong in the paper text
  )

  # Canvas size. Height is driven by whichever is TALLER: the heatmap body
  # or the legend. With ~38 metabolites the legend usually wins, and if
  # the canvas is sized from the body alone the legend gets clipped.
  n_legend_entries <- length(metabolite_levels) +
                      length(pathway_levels) +
                      sum(sapply(column_annotation, function(x)
                        length(unique(x))))
  height <- max(nrow(z) * 8 / 72 + 3.5, (n_legend_entries + 8) * 0.19 + 3)
  width  <- ncol(z) * 12 / 72 + 11

  pdf(paste0(filename, ".pdf"), width = width, height = height)
  grid::grid.newpage(); grid::grid.draw(figure$gtable); dev.off()

  png(paste0(filename, ".png"), width = width, height = height,
      units = "in", res = 180)
  grid::grid.newpage(); grid::grid.draw(figure$gtable); dev.off()

  cat(sprintf("  %s: %d genes x %d samples (%.1f x %.1f in)\n",
              filename, nrow(z), ncol(z), width, height))
  invisible(z)
}


# =====================================================================
#  5. PANEL A: HGF2 IN MONOCULTURE, EARLY vs LATE
# =====================================================================

# HGF2 monoculture samples start with "HGF2_"
# (co-culture samples start with "Consortium_")
columns_A <- grep("^HGF2_", colnames(counts_HGF2), value = TRUE)

genes_A <- genes %>% filter(genome == "HGF2")

z_A <- compute_zscore(counts_HGF2, columns_A, genes_A$locus_tag)

# Column track: which sample is Early and which is Late
annotation_A <- data.frame(
  Timepoint = if_else(str_detect(columns_A, "Late"),
                      "Late (45 h)", "Early (8 h)"),
  row.names = columns_A
)

z_A <- draw_panel(z_A, genes_A, annotation_A,
                  list(Timepoint = timepoint_colors),
                  "Fig7A_HGF2_monoculture_early_late")


# =====================================================================
#  6. PANEL B: PT33, ALL SAMPLES
#
#     Here all 12 samples ARE used, because PT33 was sequenced to a
#     comparable depth in monoculture and in co-culture.
#
#     (The same cannot be done for HGF2: in co-culture it has roughly
#     300x fewer reads than in monoculture, so most of its genes are not
#     quantifiable there.)
# =====================================================================

genes_B <- genes %>% filter(genome == "PT33")

z_B <- compute_zscore(counts_PT33, colnames(counts_PT33), genes_B$locus_tag)

# Two column tracks are needed here: timepoint and culture type
annotation_B <- data.frame(
  Timepoint = if_else(str_detect(colnames(z_B), "Late"),
                      "Late (45 h)", "Early (8 h)"),
  Culture   = if_else(str_starts(colnames(z_B), "Consortium"),
                      "Co-culture", "Monoculture"),
  row.names = colnames(z_B)
)

z_B <- draw_panel(z_B, genes_B, annotation_B,
                  list(Timepoint = timepoint_colors,
                       Culture   = culture_colors),
                  "Fig7B_PT33_complete")


# =====================================================================
#  7. SEQUENCING DEPTH CHECK
#
#     Not a figure, but worth looking at before interpreting anything:
#     if two conditions were sequenced to very different depths, part of
#     what the heatmap shows may be that imbalance rather than biology.
# =====================================================================

depth <- bind_rows(
  tibble(genome = "HGF2", sample = colnames(counts_HGF2),
         reads = colSums(counts_HGF2)),
  tibble(genome = "PT33", sample = colnames(counts_PT33),
         reads = colSums(counts_PT33))
) %>%
  mutate(culture = if_else(str_starts(sample, "Consortium"),
                           "Co-culture", "Monoculture")) %>%
  group_by(genome, culture) %>%
  summarise(median_reads = round(median(reads)), .groups = "drop")

cat("\nSequencing depth (median reads per sample):\n")
print(as.data.frame(depth), row.names = FALSE)
cat("\nHGF2 in co-culture has ~300x fewer reads than in monoculture.\n")
cat("That is why panel A uses monoculture samples only.\n")

cat("\nDone.\n")

