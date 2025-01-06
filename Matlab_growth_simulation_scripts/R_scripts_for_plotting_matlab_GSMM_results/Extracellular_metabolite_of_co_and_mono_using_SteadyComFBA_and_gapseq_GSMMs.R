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
library(dplyr)

#setwd("/home/alexis/Desktop/github_script_finales_alexis/Objective_1_Metabolic_interactions_insilico/Simulation_of_growth/R_scripts_plotting_SteadyCom_results/Extracellular_metabolite_profiles_plots_SteadyCom")
setwd("/media/alexis/hdd2/objetivo_1_tesis_doctoral/matlab_scripts_genomas")

# Mono-cultures
bt_mono <- read_excel("Bacteroides_thetaiotaomicron_VPI_5482_monocultivo_modelseed.xlsx")
anim_mono <- read_excel("Bifidobacterium_animalis_lactis_PT33_monocultivo_modelseed.xlsx")
ino_mono <- read_excel("Clostridium_innocuum_HFG2_monocultivo_modelseed.xlsx")
m38_mono <- read_excel("Lacticaseibacillus_paracasei_M38_monocultivo_modelseed.xlsx")

# Co-cultures
bt_ino <- read_excel("Bacteroides_thetaiotaomicron_VPI_5482_consorcio_con_Clostridium_innocuum_HFG2_modelseed.xlsx")
ino_bt <- read_excel("Clostridium_innocuum_HFG2_consorcio_con_Bacteroides_thetaiotaomicron_VPI_5482_modelseed.xlsx")
anim_ino <- read_excel("Bifidobacterium_animalis_lactis_PT33_consorcio_con_Clostridium_innocuum_HFG2_modelseed.xlsx")
ino_anim <- read_excel("Clostridium_innocuum_HFG2_consorcio_con_Bifidobacterium_animalis_lactis_PT33_modelseed.xlsx")
bt_anim <- read_excel("Bacteroides_thetaiotaomicron_VPI_5482_consorcio_con_Bifidobacterium_animalis_lactis_PT33_modelseed.xlsx")
anim_bt <- read_excel("Bifidobacterium_animalis_lactis_PT33_consorcio_con_Bacteroides_thetaiotaomicron_VPI_5482_modelseed.xlsx")
m38_ino <- read_excel("Lacticaseibacillus_paracasei_M38_consorcio_con_Clostridium_innocuum_HFG2_modelseed.xlsx")
ino_m38 <- read_excel("Clostridium_innocuum_HFG2_consorcio_con_Lacticaseibacillus_paracasei_M38_modelseed.xlsx")

# Filter Mono-cultures to exchange metabolites; 
bt_mono_EX <- grep("EX_", bt_mono$rxn)
bt_mono_EXmet <- bt_mono[c(bt_mono_EX), ]

anim_mono_EX <- grep("EX_", anim_mono$rxn)
anim_mono_EXmet <- anim_mono[c(anim_mono_EX), ]

ino_mono_EX <- grep("EX_", ino_mono$rxn)
ino_mono_EXmet <- ino_mono[c(ino_mono_EX), ]

m38_mono_EX <- grep("EX_", m38_mono$rxn)
m38_mono_EXmet <- m38_mono[c(m38_mono_EX), ]

bt_ino <- bt_ino %>%
  mutate_at("rxn", str_replace, "Bacteroides_thetaiotaomicron_VPI_5482I", "")
bt_ino$rxn <- gsub("\\[u]tr","(e)", bt_ino$rxn)
bt_ino$rxn <- gsub("\\[c]tr","(e)", bt_ino$rxn)

ino_bt <- ino_bt %>%
  mutate_at("rxn", str_replace, "Clostridium_innocuum_HFG2I", "")
ino_bt$rxn <- gsub("\\[u]tr","(e)", ino_bt$rxn)
ino_bt$rxn <- gsub("\\[c]tr","(e)", ino_bt$rxn)


ino_m38 <- ino_m38 %>%
  mutate_at("rxn", str_replace, "Clostridium_innocuum_HFG2I", "")
ino_m38$rxn <- gsub("\\[u]tr","(e)", ino_m38$rxn)
ino_m38$rxn <- gsub("\\[c]tr","(e)", ino_m38$rxn)

m38_ino <- m38_ino %>%
  mutate_at("rxn", str_replace, "Lacticaseibacillus_paracasei_M38I", "")
m38_ino$rxn <- gsub("\\[u]tr","(e)", m38_ino$rxn)
m38_ino$rxn <- gsub("\\[c]tr","(e)", m38_ino$rxn)


anim_ino <- anim_ino %>%
  mutate_at("rxn", str_replace, "Bifidobacterium_animalis_lactis_PT33I", "")
anim_ino$rxn <- gsub("\\[u]tr","(e)", anim_ino$rxn)
anim_ino$rxn <- gsub("\\[c]tr","(e)", anim_ino$rxn)

ino_anim <- ino_anim %>%
  mutate_at("rxn", str_replace, "Clostridium_innocuum_HFG2I", "")
ino_anim$rxn <- gsub("\\[u]tr","(e)", ino_anim$rxn)
ino_anim$rxn <- gsub("\\[c]tr","(e)", ino_anim$rxn)

bt_anim <- bt_anim %>%
  mutate_at("rxn", str_replace, "Bacteroides_thetaiotaomicron_VPI_5482I", "")
bt_anim$rxn <- gsub("\\[u]tr","(e)", bt_anim$rxn)
bt_anim$rxn <- gsub("\\[c]tr","(e)", bt_anim$rxn)

anim_bt <- anim_bt %>%
  mutate_at("rxn", str_replace, "Bifidobacterium_animalis_lactis_PT33I", "")
anim_bt$rxn <- gsub("\\[u]tr","(e)", anim_bt$rxn)
anim_bt$rxn <- gsub("\\[c]tr","(e)", anim_bt$rxn)

#---------------------------------------------------------------------------------------------------------------------------
# Modify modelseed reaction names with metabolite names
modelseed_names <- fread("modelseed_compounds.tsv", header = TRUE)
#head(modelseed_names)

#bt_mono_EXmet
#anim_mono_EXmet
#ino_mono_EXmet
#bt_ino
#ino_bt
#anim_ino
#ino_anim
#bt_anim
#anim_bt
#ino_bt_anim
#anim_bt_ino
#bt_anim_ino
#1.1
m38_mono_EXmet$rxn <- gsub("^EX_|_e0$", "", m38_mono_EXmet$rxn)
m38_mono_EXmet_ljoin <- left_join(m38_mono_EXmet, modelseed_names, by = c("rxn" = "id"))
m38_mono_EXmet_ljoin_select <- select(m38_mono_EXmet_ljoin, v, name)
m38_mono_EXmet_ljoin_select <- na.omit(m38_mono_EXmet_ljoin_select)
m38_mono_EXmet_ljoin_select <- select(m38_mono_EXmet_ljoin, name, v)
m38_mono_EXmet_ljoin_select
#1.2
m38_ino$rxn <- gsub("^EX_|_e0$", "", m38_ino$rxn)
m38_ino$rxn <- gsub("\\(e\\)", "", m38_ino$rxn)
m38_ino_ljoin <- left_join(m38_ino, modelseed_names, by = c("rxn" = "id"))
m38_ino_ljoin_select <- select(m38_ino_ljoin, v, name)
m38_ino_ljoin_select <- na.omit(m38_ino_ljoin_select)
m38_ino_ljoin_select <- select(m38_ino_ljoin, name, v)
m38_ino_ljoin_select
#1.3
ino_m38$rxn <- gsub("^EX_|_e0$", "", ino_m38$rxn)
ino_m38$rxn <- gsub("\\(e\\)", "", ino_m38$rxn)
ino_m38_ljoin <- left_join(ino_m38, modelseed_names, by = c("rxn" = "id"))
ino_m38_ljoin_select <- select(ino_m38_ljoin, v, name)
ino_m38_ljoin_select <- na.omit(ino_m38_ljoin_select)
ino_m38_ljoin_select <- select(ino_m38_ljoin, name, v)
ino_m38_ljoin_select
#1
bt_mono_EXmet$rxn <- gsub("^EX_|_e0$", "", bt_mono_EXmet$rxn)
bt_mono_EXmet_ljoin <- left_join(bt_mono_EXmet, modelseed_names, by = c("rxn" = "id"))
bt_mono_EXmet_ljoin_select <- select(bt_mono_EXmet_ljoin, v, name)
bt_mono_EXmet_ljoin_select <- na.omit(bt_mono_EXmet_ljoin_select)
bt_mono_EXmet_ljoin_select <- select(bt_mono_EXmet_ljoin, name, v)
bt_mono_EXmet_ljoin_select
#2
anim_mono_EXmet$rxn <- gsub("^EX_|_e0$", "", anim_mono_EXmet$rxn)
anim_mono_EXmet_ljoin <- left_join(anim_mono_EXmet, modelseed_names, by = c("rxn" = "id"))
anim_mono_EXmet_ljoin_select <- select(anim_mono_EXmet_ljoin, v, name)
anim_mono_EXmet_ljoin_select <- na.omit(anim_mono_EXmet_ljoin_select)
anim_mono_EXmet_ljoin_select <- select(anim_mono_EXmet_ljoin, name, v)
anim_mono_EXmet_ljoin_select
#3
ino_mono_EXmet$rxn <- gsub("^EX_|_e0$", "", ino_mono_EXmet$rxn)
ino_mono_EXmet_ljoin <- left_join(ino_mono_EXmet, modelseed_names, by = c("rxn" = "id"))
ino_mono_EXmet_ljoin_select <- select(ino_mono_EXmet_ljoin, v, name)
ino_mono_EXmet_ljoin_select <- na.omit(ino_mono_EXmet_ljoin_select)
ino_mono_EXmet_ljoin_select <- select(ino_mono_EXmet_ljoin, name, v)
ino_mono_EXmet_ljoin_select
#4
bt_ino$rxn <- gsub("^EX_|_e0$", "", bt_ino$rxn)
bt_ino$rxn <- gsub("\\(e\\)", "", bt_ino$rxn)
bt_ino_ljoin <- left_join(bt_ino, modelseed_names, by = c("rxn" = "id"))
bt_ino_ljoin_select <- select(bt_ino_ljoin, v, name)
bt_ino_ljoin_select <- na.omit(bt_ino_ljoin_select)
bt_ino_ljoin_select <- select(bt_ino_ljoin, name, v)
bt_ino_ljoin_select
#5
ino_bt$rxn <- gsub("^EX_|_e0$", "", ino_bt$rxn)
ino_bt$rxn <- gsub("\\(e\\)", "", ino_bt$rxn)
ino_bt_ljoin <- left_join(ino_bt, modelseed_names, by = c("rxn" = "id"))
ino_bt_ljoin_select <- select(ino_bt_ljoin, v, name)
ino_bt_ljoin_select <- na.omit(ino_bt_ljoin_select)
ino_bt_ljoin_select <- select(ino_bt_ljoin, name, v)
ino_bt_ljoin_select
#6
anim_ino$rxn <- gsub("^EX_|_e0$", "", anim_ino$rxn)
anim_ino$rxn <- gsub("\\(e\\)", "", anim_ino$rxn)
anim_ino_ljoin <- left_join(anim_ino, modelseed_names, by = c("rxn" = "id"))
anim_ino_ljoin_select <- select(anim_ino_ljoin, v, name)
anim_ino_ljoin_select <- na.omit(anim_ino_ljoin_select)
anim_ino_ljoin_select <- select(anim_ino_ljoin, name, v)
anim_ino_ljoin_select
#7
ino_anim$rxn <- gsub("^EX_|_e0$", "", ino_anim$rxn)
ino_anim$rxn <- gsub("\\(e\\)", "", ino_anim$rxn)
ino_anim_ljoin <- left_join(ino_anim, modelseed_names, by = c("rxn" = "id"))
ino_anim_ljoin_select <- select(ino_anim_ljoin, v, name)
ino_anim_ljoin_select <- na.omit(ino_anim_ljoin_select)
ino_anim_ljoin_select <- select(ino_anim_ljoin, name, v)
ino_anim_ljoin_select
#8
bt_anim$rxn <- gsub("^EX_|_e0$", "", bt_anim$rxn)
bt_anim$rxn <- gsub("\\(e\\)", "", bt_anim$rxn)
bt_anim_ljoin <- left_join(bt_anim, modelseed_names, by = c("rxn" = "id"))
bt_anim_ljoin_select <- select(bt_anim_ljoin, v, name)
bt_anim_ljoin_select <- na.omit(bt_anim_ljoin_select)
bt_anim_ljoin_select <- select(bt_anim_ljoin, name, v)
bt_anim_ljoin_select
#9
anim_bt$rxn <- gsub("^EX_|_e0$", "", anim_bt$rxn)
anim_bt$rxn <- gsub("\\(e\\)", "", anim_bt$rxn)
anim_bt_ljoin <- left_join(anim_bt, modelseed_names, by = c("rxn" = "id"))
anim_bt_ljoin_select <- select(anim_bt_ljoin, v, name)
anim_bt_ljoin_select <- na.omit(anim_bt_ljoin_select)
anim_bt_ljoin_select <- select(anim_bt_ljoin, name, v)
anim_bt_ljoin_select

#---------------------------------------------------------------------------------------------------------------------------
## Merge bacteria tables
m38_ino_data <- full_join(m38_ino_ljoin_select, ino_m38_ljoin_select, by = 'name')
colnames(m38_ino_data)[2] = "Co-culture L. paracasei M38"
colnames(m38_ino_data)[3] = "Co-culture Clostridium inn. HGF2"

m38_ino_data <- full_join(m38_ino_data, m38_mono_EXmet_ljoin_select, by = 'name')
colnames(m38_ino_data)[4] = "Mono-culture L. paracasei M38"
m38_ino_data <- full_join(m38_ino_data, ino_mono_EXmet_ljoin_select, by = 'name')
colnames(m38_ino_data)[5] = "Mono-culture Clostridium inn. HGF2"

m38_ino_data <- m38_ino_data[complete.cases(m38_ino_data$name), ]
# replace NAs with 0
m38_ino_data[is.na(m38_ino_data)] <- 0

# first column into rownames
m38_ino_data2 <- m38_ino_data %>% remove_rownames %>% column_to_rownames(var="name")

# delete rows by row names
m38_ino_data3 <- m38_ino_data2[!(row.names(m38_ino_data2) %in% c("EX_biomass(e)")),]

#transpose data
m38_ino_data4 <- t(m38_ino_data3)
m38_ino_data41 <- as.data.frame(m38_ino_data4)

str(m38_ino_data41)
colnames(m38_ino_data41)
rownames(m38_ino_data41)

min(m38_ino_data41)
max(m38_ino_data41)

# Define a new color palette with white in the middle
nueva_paleta <- colorRampPalette(c("navy", "white", "firebrick3"))(41)

# Set the cut-off points on the color scale
breaks <- c(seq(min(m38_ino_data41), -0.00001, length.out = 20), 0, seq(0.00001, max(m38_ino_data41), length.out = 20))

# Plot the heatmap with the new color palette and the cut points
p3 <- pheatmap(m38_ino_data41, 
               main = "SteadyCom-FBA: Co-culture vs Mono-culture",
               clustering_distance_cols = "euclidean", 
               clustering_distance_rows = "euclidean", 
               clustering_method = "average",
               #cluster_rows = FALSE,
               #cluster_cols = FALSE,
               cellwidth = 13, 
               cellheight = 20,
               col = nueva_paleta,  # Use the new color palette with white in the center.
               breaks = breaks,      # Establish cut-off points
               fontsize = 13)

######################################################################################################################################################
# Graph with data normalized by the scale function without centering
m38_ino_data41
str(m38_ino_data41)

m38_ino_data_scale = scale(m38_ino_data4, center = FALSE, scale = TRUE)

m38_ino_data_scale
write.csv(m38_ino_data_scale, "m38_ino_data_scale_2.csv", row.names = TRUE)

nueva_paleta_2 <- colorRampPalette(c("navy", "white", "firebrick3"))(219)

p4 <- pheatmap(m38_ino_data_scale, 
               main = "SteadyCom-FBA: Co-culture vs Mono-culture",
               clustering_distance_cols = "euclidean", 
               clustering_distance_rows = "euclidean", 
               clustering_method = "average",
               cluster_rows = FALSE,
               #cluster_cols = FALSE,
               cellwidth = 13, 
               cellheight = 20,
               col = nueva_paleta_2,  # Use the new color palette with white in the center
               fontsize = 13)

p4

ggsave(('Steadycom_FBA_gapseq_m38_ino_2.pdf'), plot=p4, dpi = 600, width = 23, height = 10)

####################################################################################################################################################################
## Merge bacteria tables
bt_ino_data <- full_join(bt_ino_ljoin_select, ino_bt_ljoin_select, by = 'name')
colnames(bt_ino_data)[2] = "Co-culture Bacteroides the. VPI-5482"
colnames(bt_ino_data)[3] = "Co-culture Clostridium inn. HGF2"

bt_ino_data <- full_join(bt_ino_data, bt_mono_EXmet_ljoin_select, by = 'name')
colnames(bt_ino_data)[4] = "Mono-culture Bacteroides the. VPI-5482"

bt_ino_data <- full_join(bt_ino_data, ino_mono_EXmet_ljoin_select, by = 'name')
colnames(bt_ino_data)[5] = "Mono-culture Clostridium inn. HGF2"

bt_ino_data <- bt_ino_data[complete.cases(bt_ino_data$name), ]
# replace NAs with 0
bt_ino_data[is.na(bt_ino_data)] <- 0

# first column into rownames
bt_ino_data2 <- bt_ino_data %>% remove_rownames %>% column_to_rownames(var="name")

# delete rows by row names
bt_ino_data3 <- bt_ino_data2[!(row.names(bt_ino_data2) %in% c("EX_biomass(e)")),]

#transpose data
bt_ino_data4 <- t(bt_ino_data3)
bt_ino_data41 <- as.data.frame(bt_ino_data4)

str(bt_ino_data41)
colnames(bt_ino_data41)
rownames(bt_ino_data41)

min(bt_ino_data41)
max(bt_ino_data41)

# Define a new color palette with white in the middle
nueva_paleta <- colorRampPalette(c("navy", "white", "firebrick3"))(41)

# Set the cut-off points on the color scale
breaks <- c(seq(min(bt_ino_data41), -0.00001, length.out = 20), 0, seq(0.00001, max(bt_ino_data41), length.out = 20))

# Plot the heatmap with the new color palette and the cut points
p1 <- pheatmap(bt_ino_data41, 
               main = "SteadyCom-FBA: Co-culture vs Mono-culture",
               clustering_distance_cols = "euclidean", 
               clustering_distance_rows = "euclidean", 
               clustering_method = "average",
               #cluster_rows = FALSE,
               #cluster_cols = FALSE,
               cellwidth = 13, 
               cellheight = 20,
               col = nueva_paleta,  # Use the new color palette with white in the center
               breaks = breaks,      # Establish cut-off points
               fontsize = 13)

######################################################################################################################################################
# Graph with data normalized by the scale function without centering
bt_ino_data41
str(bt_ino_data41)

bt_ino_data_scale = scale(bt_ino_data4, center = FALSE, scale = TRUE)

bt_ino_data_scale

nueva_paleta_2 <- colorRampPalette(c("navy", "white", "firebrick3"))(219)

p2 <- pheatmap(bt_ino_data_scale, 
               main = "SteadyCom-FBA: Co-culture vs Mono-culture",
               clustering_distance_cols = "euclidean", 
               clustering_distance_rows = "euclidean", 
               clustering_method = "average",
               cluster_rows = FALSE,
               #cluster_cols = FALSE,
               cellwidth = 13, 
               cellheight = 20,
               col = nueva_paleta_2,  # Use the new color palette with white in the center
               fontsize = 13)

p2

ggsave(('Steadycom_FBA_gapseq_bt_ino_2.pdf'), plot=p2, dpi = 600, width = 23, height = 10)

#---------------------------------------------------------------------------------------------------------------------------

anim_ino_data <- full_join(anim_ino_ljoin_select, ino_anim_ljoin_select, by = 'name')
colnames(anim_ino_data)[2] = "Co-culture B. animalis lactis PT33"
# change column names
colnames(anim_ino_data)[3] = "Co-culture Clostridium inn. HGF2"

anim_ino_data <- full_join(anim_ino_data, anim_mono_EXmet_ljoin_select, by = 'name')
colnames(anim_ino_data)[4] = "Mono-culture B. animalis lactis PT33"

anim_ino_data <- full_join(anim_ino_data, ino_mono_EXmet_ljoin_select, by = 'name')
colnames(anim_ino_data)[5] = "Mono-culture Clostridium inn. HGF2"

anim_ino_data <- anim_ino_data[complete.cases(anim_ino_data$name), ]
anim_ino_data[is.na(anim_ino_data)] <- 0


# first column into rownames
anim_ino_data2 <- anim_ino_data %>% remove_rownames %>% column_to_rownames(var="name")

# delete rows by row names:
anim_ino_data3 <- anim_ino_data2[!(row.names(anim_ino_data2) %in% c("EX_biomass(e)")),]

#transpose data
anim_ino_data4 <- t(anim_ino_data3)
anim_ino_data41 <- as.data.frame(anim_ino_data4)

str(anim_ino_data4)

#ruta_archivo_csv <- "/home/alexis/Desktop/figuras_tesis_doctoral_importantes/REVISAR_DATA_STEADYCOMFBA.csv"
#write.csv(anim_ino_data41, file = ruta_archivo_csv, row.names = TRUE)

# scale function:
anim_ino_data_scale = scale(anim_ino_data4, center = FALSE, scale = TRUE)

nueva_paleta_2 <- colorRampPalette(c("navy", "white", "firebrick3"))(219)

p5 <- pheatmap(anim_ino_data_scale, 
               main = "SteadyCom-FBA: Co-culture vs Mono-culture",
               clustering_distance_cols = "euclidean", 
               clustering_distance_rows = "euclidean", 
               clustering_method = "average",
               cluster_rows = FALSE,
               #cluster_cols = FALSE,
               cellwidth = 13, 
               cellheight = 20,
               col = nueva_paleta_2,  # Use the new color palette with white in the center
               fontsize = 13)

p5

ggsave(('Steadycom_FBA_gapseq_anim_ino_2.pdf'), plot=p5, dpi = 600, width = 23, height = 10)
ggsave(('Steadycom_FBA_gapseq_anim_ino_2.png'), plot=p5, dpi = 600, width = 23, height = 10)

#--------------------------------------------------------------------------------------------------------------------------

bt_anim_data <- full_join(bt_anim, anim_bt, by = 'name')
colnames(bt_anim_data)[2] = "Co-culture Bacteroides the. VPI-5482"
# change column names:
colnames(bt_anim_data)[3] = "Co-culture B. animalis lactis PT33"

bt_anim_data <- full_join(bt_anim_data, bt_mono_EXmet_ljoin_select, by = 'name')
colnames(bt_anim_data)[4] = "Mono-culture Bacteroides the. VPI-5482"

bt_anim_data <- full_join(bt_anim_data, anim_mono_EXmet_ljoin_select, by = 'name')
colnames(bt_anim_data)[5] = "Mono-culture B. animalis lactis PT33"

# replace NAs with 0:
bt_anim_data[is.na(bt_anim_data)] <- 0

# first column into rownames:
bt_anim_data2 <- bt_anim_data %>% remove_rownames %>% column_to_rownames(var="rxn")


# delete rows by row names:
bt_anim_data3 <- bt_anim_data2[!(row.names(bt_anim_data2) %in% c("EX_biomass(e)")),]


#transpose data
bt_anim_data4 <- t(bt_anim_data3)
bt_anim_data41 <- as.data.frame(bt_anim_data4)

# scale function:
bt_anim_data_scale = scale(bt_anim_data4, center = FALSE, scale = TRUE)

pheatmap(bt_anim_data_scale, cutree_rows = 4, main = "SteadyCom-FBA: Co-culture vs Mono-culture",
         clustering_distance_cols = "maximum", clustering_distance_rows = "maximum", 
         clustering_method = "single",
         cluster_rows=FALSE,
         cellwidth = 9, cellheight = 20,
         treeheight_row = 0,
         treeheight_column=0,
         col = colorRampPalette(c("navy", "white", "firebrick3"))(2019),
         legend_breaks = c(1.5, 0.85, 0, -0.85, -1.5),
         legend_labels = c("1.5", "(Produccion)", 0, "(Consumo)", "-1.5"))

