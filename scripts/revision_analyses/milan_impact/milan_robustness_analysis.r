library(readxl)
library(dplyr)
library(stringr)
library(sf)
library(geosphere)
library(factoextra)
library(ClustGeo)
library(mclust)
library(purrr)
library(ggplot2)
library(tidyr)

setwd("/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering/")

# --- 1. SETTINGS & PARAMETERS ---
AOI <- "Lombardy"
# Ensure variabili_clustering is defined, e.g.:
# variabili_clustering <- c("Var1", "Var2", ...) 


variabili_clustering <- c(
  "Local_Units_Industry_Suspended", "Local_Units_Services_Suspended",
  "Value_Added_Industry_Suspended", "Value_Added_Services_Suspended",
  "Employed_Industry_Suspended", "Employed_Services_Suspended",
  "Local_Units_Industry_Active", "Local_Units_Services_Active",
  "Value_Added_Industry_Active", "Value_Added_Services_Active",
  "Employed_Industry_Active", "Employed_Services_Active")

# --- 2. LOAD DATA ---
DPCM <- read_excel("data/DPCM/DPCM_cleaned.xlsx")

ManualMatchLAUs <- read_excel("data/anagrafica_comuni/Manually_matched_LAUs_Oct2023.xlsx", sheet = "Modified")
ManualMatchLAUs <- ManualMatchLAUs %>%
  mutate(across(c(`Denominazione Comune NEW`), toupper),
         across(c(`Denominazione Comune NEW`), stringr::str_to_title)) %>%
  dplyr::select(LAU_Code_old = `Codice Istat Comune`,
                LAU_Code_new = `Codice Istat Comune NEW`) %>%
  distinct()

load("data/anagrafica_comuni/LAU2021_Eurostat_Italy.RData")

# --- 3. DATA CLEANING & MERGING ---
DPCM_red <- DPCM %>%
  dplyr::rename(LAU_Code = ISTAT_Code_LAU) %>%
  dplyr::select(-Name_LAU)
DPCM_red <- DPCM_red %>%
  group_by(LAU_Code) %>%
  mutate(across(is.numeric, function(x) sum(x,na.rm=T))) %>%
  ungroup()
DPCM_red <- left_join(x = LAU_IT_sf, y = DPCM_red, by = c("LAU_Code"))
DPCM_red <- DPCM_red %>% ungroup() %>% filter(!duplicated(LAU_Code))
DPCM_red <- DPCM_red %>%
  filter(NUTS2_Name == "Lombardia")

DPCM_prima_Milano <- DPCM_red %>% filter(LAU_Code != "015146")


# --- 5. FUNCTION TO PREPARE DISTANCE MATRICES ---
prepare_matrices <- function(data, includi_milano, vars_clustering) {
  if (!includi_milano) {
    data <- data %>% filter(LAU_Code != "015146") # Milan code
  }
  
  CData <- data %>%
    filter(!is.na(Local_Units_Industry_Active)) %>%
    dplyr::select(Area_km2, LAU_Code, Lon_centr_LAU, Lat_centr_LAU,
                  Alt_mean, Alt_min, Alt_max,
                  all_of(vars_clustering)) %>%
    mutate(across(c(is.numeric, -contains("Alt"), -contains("Lon"), -contains("Lat")), 
                  function(x) x / Area_km2)) %>%
    dplyr::select(-Area_km2) %>%
    mutate(across(c(is.numeric, -contains("Lon"), -contains("Lat")), 
                  function(x) (x - median(x, na.rm = TRUE)) / (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))))
  
  st_geometry(CData) <- NULL
  
  Xmat_geo <- as.matrix(CData[, c("Lon_centr_LAU", "Lat_centr_LAU")])
  Xmat_nogeo <- as.matrix(CData[, vars_clustering])
  
  D0 <- get_dist(Xmat_nogeo)
  distanza_WGS84 <- function(p1, p2) geosphere::distGeo(p1, p2, a = 6378137, f = 1 / 298.257223563)
  D1 <- as.dist(geosphere::distm(Xmat_geo, fun = distanza_WGS84))
  
  # Standardize and explicitly cast back to dist
  D0_std <- as.dist(as.matrix(D0) / max(D0))
  D1_std <- as.dist(as.matrix(D1) / max(D1))
  
  list(D0 = D0_std, D1 = D1_std, LAU = CData$LAU_Code)
}

# --- 6. GENERATE MATRICES & COMPUTE HCLUSTGEO ---
mat_with <- prepare_matrices(DPCM_red, TRUE, variabili_clustering)
mat_without <- prepare_matrices(DPCM_prima_Milano, FALSE, variabili_clustering)

griglia_alpha <- seq(0, 1, 0.1)
griglia_K <- 2:10

trees_with <- purrr::map(griglia_alpha, ~ hclustgeo(mat_with$D0, mat_with$D1, alpha = .x))
trees_without <- purrr::map(griglia_alpha, ~ hclustgeo(mat_without$D0, mat_without$D1, alpha = .x))

# --- 7. COMPUTE ARI INDEX ---
ari_results <- expand.grid(Alpha = griglia_alpha, K = griglia_K) %>%
  mutate(ARI = NA_real_)

for (i in seq_along(griglia_alpha)) {
  for (k in griglia_K) {
    part_with <- cutree(trees_with[[i]], k)
    part_without <- cutree(trees_without[[i]], k)
    
    names(part_with) <- mat_with$LAU
    names(part_without) <- mat_without$LAU
    
    # Align partitions (subsetting 'with Milan' vector to match 'without Milan' vector)
    part_with_aligned <- part_with[names(part_without)]
    
    idx <- which(ari_results$Alpha == griglia_alpha[i] & ari_results$K == k)
    ari_results$ARI[idx] <- adjustedRandIndex(part_with_aligned, part_without)
  }
}

# --- 8. PLOT HEATMAP ---
dati_heatmap_ari <- ari_results %>%
  mutate(
    K_factor = factor(K, levels = rev(sort(unique(K)))),
    Alpha_factor = factor(Alpha)
  )

min_ari <- min(0, min(dati_heatmap_ari$ARI, na.rm = TRUE))

p_ari_heatmap <- ggplot(dati_heatmap_ari, aes(x = Alpha_factor, y = K_factor, fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(
    label = sprintf("%.3f", ARI),
    colour = ifelse(ARI > 0.65 | ARI < -0.2, "white", "#243746")
  ), size = 3) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#B97572", mid = "#F3F2EF", high = "#365F7C",
    midpoint = 0.5, limits = c(min_ari, 1), name = "ARI"
  ) +
  coord_equal() +
  labs(
    title = "ARI: With vs. Without Milan (Lombardy)",
    subtitle = "Stability of partitions across K and Alpha",
    x = "Alpha",
    y = "Number of Clusters (K)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", colour = "#183746"),
    plot.subtitle = element_text(colour = "#5F6572")
  )

p_ari_heatmap

ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/milan_clustering", "ARI_milan.png"),
       p_ari_heatmap, width = 14, height = 7.5, dpi = 200, bg = "white")
ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/milan_clustering", "ARI_milan.pdf"),
       p_ari_heatmap, width = 14, height = 7.5)