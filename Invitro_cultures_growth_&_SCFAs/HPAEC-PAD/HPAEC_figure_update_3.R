## Código Simplificado - Cromatogramas con 51h en 0 y resto stackeado
library(tidyverse)
library(patchwork)
library(ggplot2)

# ===== FUNCIÓN PARA PROCESAR CONTROLES =====
process_controls <- function(directory, offsets) {
  setwd(directory)
  files <- list.files(pattern="*.txt")
  
  data_list <- lapply(1:length(files), function(i) {
    x <- files[i]
    dat <- read.table(x, header=FALSE, sep="\t", skip=43)
    names(dat) <- c("time", "step", "value")
    dat$Sample <- gsub(".txt", "", x)
    
    # Crear nombres descriptivos y asignar offsets específicos para 8 controles
    if(grepl("Inulin", x)) { 
      dat$time_label <- "Inulin"
      dat$offset <- 0  # Inulin siempre arriba
    }
    else if(grepl("Fructose", x)) { 
      dat$time_label <- "Fructose"
      dat$offset <- -15
    }
    else if(grepl("1,1,1-Kestose pentaose", x)) { 
      dat$time_label <- "Pentaose"
      dat$offset <- -30
    }
    else if(grepl("1,1-Kestose tetraose", x)) { 
      dat$time_label <- "Tetraose"
      dat$offset <- -45
    }
    else if(grepl("1-Kestose", x)) { 
      dat$time_label <- "1-Kestose"
      dat$offset <- -60
    }
    else if(grepl("a-D-Glucose", x)) { 
      dat$time_label <- "Glucose"
      dat$offset <- -75
    }
    else if(grepl("Isomalto-oligosaccharides", x)) { 
      dat$time_label <- "Isomalto-oligo"
      dat$offset <- -90
    }
    else if(grepl("Malto-oligosaccharide", x)) { 
      dat$time_label <- "Malto-oligo"
      dat$offset <- -105  # Malto-oligo al fondo
    }
    else { 
      dat$time_label <- paste("Control", i)
      dat$offset <- offsets[i]
    }
    
    dat$plot <- dat$value + dat$offset
    return(dat)
  })
  
  return(do.call("rbind", data_list))
}

# ===== FUNCIÓN PARA PROCESAR ARCHIVOS DE MUESTRAS =====
process_samples <- function(directory, offsets) {
  setwd(directory)
  files <- list.files(pattern="*.txt")
  
  data_list <- lapply(files, function(x) {
    dat <- read.table(x, header=FALSE, sep="\t", skip=43)
    names(dat) <- c("time", "step", "value")
    dat$Sample <- gsub(".txt", "", x)
    
    # Mapear tiempos reales
    if(grepl("15-", x)) { dat$time_label <- "51 hours"; dat$offset <- offsets[1] }
    else if(grepl("12-", x)) { dat$time_label <- "45 hours"; dat$offset <- offsets[2] }
    else if(grepl("4-", x)) { dat$time_label <- "8 hours"; dat$offset <- offsets[3] }
    else if(grepl("1-", x)) { dat$time_label <- "2 hours"; dat$offset <- offsets[4] }
    else if(grepl("0-", x)) { dat$time_label <- "0 hours"; dat$offset <- offsets[5] }
    else { dat$time_label <- "unknown"; dat$offset <- 0 }
    
    dat$plot <- dat$value + dat$offset
    return(dat)
  })
  
  return(do.call("rbind", data_list))
}

# ===== PROCESAR TODOS LOS DATASETS =====
base_dir <- "/home/alexis/Desktop/GITHUB_REspaldos_linx_2025/Objective_3_cocultures_vs_monocultures_RNAseq/HPAECPAD_actualizacion_mayo_2025/HPAEC-PAD"

# Offsets: 51h en 0, resto hacia abajo
sample_offsets <- c(0, -20, -40, -60, -80)
control_offsets <- c(0, -15, -30, -45, -60, -75, -90, -105)  # 8 offsets para 8 controles

# Procesar datasets
control <- process_controls(paste0(base_dir, "/controls"), control_offsets)
hgf2 <- process_samples(paste0(base_dir, "/clostridium"), sample_offsets)
pt33 <- process_samples(paste0(base_dir, "/bifidobacterium"), sample_offsets)
cocult <- process_samples(paste0(base_dir, "/co_culture"), sample_offsets)

# ===== FUNCIÓN PARA CREAR GRÁFICOS DE CONTROLES =====
create_control_plot <- function(data, title) {
  # Colores para 8 controles
  control_colors <- c("dodgerblue3", "goldenrod2", "firebrick3", "slategray3", 
                      "darkgreen", "purple", "orange", "brown")
  
  p <- ggplot(data, aes(x=time, y=plot, color=time_label)) +
    geom_line(size=1) +
    scale_color_manual(values=control_colors[1:length(unique(data$time_label))], name=NULL) +
    labs(title=title, x="Time (min)", y="Value (nC)") +
    ylim(c(-120, 200)) +  # Rango ajustado para 8 controles
    scale_x_continuous(breaks=seq(0, 120, 20), limits=c(0, 120)) +  # Mismo rango que las muestras
    scale_y_continuous(breaks=seq(0, 200, 40), limits=c(-120, 200)) +
    theme_linedraw(base_size=16) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(face="bold"),
          axis.title = element_text(face="bold"),
          axis.text = element_text(face="bold"),
          legend.text = element_text(face="bold", size=12)) +  # Texto más pequeño para más leyendas
    coord_fixed(ratio=0.15)
  
  return(p)
}

# ===== FUNCIÓN PARA CREAR GRÁFICOS DE MUESTRAS =====
create_sample_plot <- function(data, title) {
  # Crear factor ordenado para leyenda
  data$time_label_factor <- factor(data$time_label, 
                                   levels = c("51 hours", "45 hours", "8 hours", "2 hours", "0 hours"))
  
  colors <- c("red", "darkblue", "darkorange", "plum3", "seagreen3")
  
  p <- ggplot(data, aes(x=time, y=plot, color=time_label_factor)) +
    geom_line(size=1) +
    scale_color_manual(values=colors, name=NULL) +
    labs(title=title, x="Time (min)", y="Value (nC)") +
    ylim(c(-100, 200)) +
    scale_x_continuous(breaks=seq(0, 180, 20)) +
    scale_y_continuous(breaks=seq(0, 200, 40), limits=c(-100, 200)) +
    theme_linedraw(base_size=16) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(face="bold"),
          axis.title = element_text(face="bold"),
          axis.text = element_text(face="bold"),
          legend.text = element_text(face="bold")) +
    coord_fixed(ratio=0.15)
  
  return(p)
}

# ===== CREAR GRÁFICOS =====
plot_control <- create_control_plot(control, "Controls")
plot_hgf2 <- create_sample_plot(hgf2, "Clostridium sp. HGF2 mono-culture")
plot_bifido <- create_sample_plot(pt33, "Bi. animalis subsp. lactis PT33 mono-culture")
plot_coculture <- create_sample_plot(cocult, "HGF2 and PT33 Co-culture")

# ===== COMBINAR Y GUARDAR =====
combined_plot <- plot_control / plot_hgf2 / plot_bifido / plot_coculture + 
  plot_annotation(tag_levels='A') & 
  theme(legend.position="right", plot.background=element_rect(fill="white", colour=NA))

combined_plot_nc <- plot_hgf2 / plot_bifido / plot_coculture + 
  plot_annotation(tag_levels='A') & 
  theme(legend.position="right", plot.background=element_rect(fill="white", colour=NA))

# Mostrar gráficos
print(combined_plot)
print(combined_plot_nc)

# Guardar
setwd(base_dir)
ggsave("figure_final_51h_at_zero.png", combined_plot, width=16, height=12, dpi=600, bg="white")
ggsave("figure_final_51h_at_zero_nc.png", combined_plot_nc, width=16, height=12, dpi=600, bg="white")

cat("Gráficos guardados exitosamente!\n")

