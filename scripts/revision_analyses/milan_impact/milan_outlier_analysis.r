##############################################
########## Setup and Configuration ###########
##############################################

base_path <- '/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering'
output_dir <- "/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/output"
setwd(base_path)

library(sf)
library(tidyverse)
library(readxl)
library(cluster)
library(factoextra) 
library(NbClust)
library(ClustGeo)
library(viridis)
library(ggpubr)
library(latex2exp)

# Define Analysis Parameters
TARGET_REGION <- "Lombardia"
MILAN_LAU_CODE <- "015146"
EXCLUDE_MILAN <- FALSE  # Set to FALSE for the robustness check
ALPHA_SPATIAL <- 0.1
K_CLUSTERS <- 8

##############################################
########## Data Loading and Merging ##########
##############################################

# Load raw DPCM data
dpcm_raw <- read_excel("./data/DPCM/DPCM_cleaned.xlsx") %>%
  dplyr::rename(LAU_Code = ISTAT_Code_LAU) %>%
  dplyr::select(-Name_LAU)

# Load and format municipal registry for code matching
lau_matches <- read_excel("./data/anagrafica_comuni/Manually_matched_LAUs_Oct2023.xlsx", sheet = "Modified") %>%
  mutate(across(c(`Denominazione Comune NEW`), toupper),
         across(c(`Denominazione Comune NEW`), stringr::str_to_title)) %>%
  dplyr::select(LAU_Code_old = `Codice Istat Comune`,
                LAU_Code_new = `Codice Istat Comune NEW`) %>%
  distinct()

# Update outdated LAU codes
old_codes_padded <- stringr::str_pad(string = lau_matches$LAU_Code_old, width = 6, pad = "0")
new_codes_padded <- stringr::str_pad(string = lau_matches$LAU_Code_new, width = 6, pad = "0")

for (i in seq_len(nrow(dpcm_raw))) {
  match_idx <- which(dpcm_raw$LAU_Code[i] == old_codes_padded)
  if (length(match_idx) > 0) {
    dpcm_raw$LAU_Code[i] <- new_codes_padded[match_idx]
  }
}

# Aggregate by new LAU codes to resolve duplicates
dpcm_aggregated <- dpcm_raw %>%
  group_by(LAU_Code) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup()

# Load geographic boundaries and merge with economic data
load("./data/anagrafica_comuni/LAU2021_Eurostat_Italy.RData") 
geo_econ_data <- left_join(x = LAU_IT_sf, y = dpcm_aggregated, by = "LAU_Code")

##############################################
########## Filtering & Base Dataset ##########
##############################################

# Filter for the target region and remove missing observations
regional_data <- geo_econ_data %>%
  filter(NUTS2_Name == TARGET_REGION) %>%
  filter(!is.na(Local_Units_Industry_Active))

# Explicitly handle Milan based on the configuration flag
if (EXCLUDE_MILAN) {
  analysis_dataset <- regional_data %>% filter(LAU_Code != MILAN_LAU_CODE)
} else {
  analysis_dataset <- regional_data
}

##############################################
########## Double Normalization ##############
##############################################

# Step 1: Normalize economic variables by Area (Density)
# Step 2: Robust standardization using Median and IQR
normalized_data <- analysis_dataset %>%
  dplyr::select(Area_km2, LAU_Code, Lon_centr_LAU, Lat_centr_LAU,
                Alt_mean, Alt_min, Alt_max,
                contains("Suspended"), contains("Active")) %>%
  mutate(across(c(where(is.numeric), -contains("Alt"), -contains("Lon"), -contains("Lat")), 
                ~ .x / Area_km2)) %>%
  dplyr::select(-Area_km2) %>%
  mutate(across(c(where(is.numeric), -contains("Lon"), -contains("Lat")), 
                ~ (.x - median(.x, na.rm = TRUE)) / (IQR(.x, na.rm = TRUE))))

##############################################
########## Distance Matrices #################
##############################################

# Separate coordinate data (geographic) from economic data
st_geometry(normalized_data) <- NULL
matrix_geo  <- as.matrix(normalized_data %>% dplyr::select(Lon_centr_LAU, Lat_centr_LAU))
matrix_econ <- as.matrix(normalized_data %>% dplyr::select(-c(1:6)))

# Calculate Euclidean distances
dist_econ <- get_dist(matrix_econ)
dist_geo  <- get_dist(matrix_geo)

# Extract row-wise average distances for mapping
analysis_dataset$AvgDist_Econ <- rowMeans(as.matrix(dist_econ))
analysis_dataset$AvgDist_Geo  <- rowMeans(as.matrix(dist_geo))

##############################################
########## Hierarchical Clustering ###########
##############################################

n_obs <- nrow(normalized_data)
weights <- rep(1 / n_obs, n_obs)

# Model 1: Non-Spatial Clustering (alpha = 0)
tree_non_spatial <- hclustgeo(D0 = dist_econ, wt = weights)

# Model 2: Spatially-Informed Clustering (alpha = ALPHA_SPATIAL)
tree_spatial <- hclustgeo(D0 = dist_econ, D1 = dist_geo, alpha = ALPHA_SPATIAL, wt = weights)

# Assign chosen K clusters to the main spatial dataframe
analysis_dataset$Cluster_NonSpatial <- as.factor(cutree(tree_non_spatial, K_CLUSTERS))
analysis_dataset$Cluster_Spatial    <- as.factor(cutree(tree_spatial, K_CLUSTERS))

##############################################
########## Export & Visualization ############
##############################################

# Colors
cluster_colors <- c("#996633", "#FF9900", "#CC0000", "#0000FF", "#3399FF", "#FFCC33", "#333333", "#66CC33")

# Plot Non-Spatial
plot_non_spatial <- analysis_dataset %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Cluster_NonSpatial)) + 
  scale_fill_manual(name = "Cluster", values = cluster_colors) + 
  theme_bw() + 
  labs(title = "Ward-like hierarchical clustering (Non-Spatial)",
       subtitle = paste0("K = ", K_CLUSTERS, " clusters (alpha = 0)"))

# Plot Spatial
plot_spatial <- analysis_dataset %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Cluster_Spatial)) + 
  scale_fill_manual(name = "Cluster", values = cluster_colors) + 
  theme_bw() + 
  labs(title = "Geographically-informed hierarchical clustering",
       subtitle = paste0("K = ", K_CLUSTERS, " clusters (alpha = ", ALPHA_SPATIAL, ")"))




filename_nospat <- paste0("Map_NonSpatial_", "EXCLUDE_MILAN_", as.character(EXCLUDE_MILAN), ".png")
filename_spat <- paste0("Map_Spatial_", "EXCLUDE_MILAN_", as.character(EXCLUDE_MILAN), ".png")

ggpubr::ggexport(plot_non_spatial, width = 1800, height = 1200, res = 150, filename = filename_nospat)
ggpubr::ggexport(plot_spatial, width = 1800, height = 1200, res = 150, filename = filename_spat)

# save.image(file = "Clustering_Workspace_Clean.RData")