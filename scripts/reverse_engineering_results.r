library(sf)
library(tidyverse)
library(readxl)
library(pracma)
library(cluster)
library(stats)
library(fpc)
library(factoextra) 
library(NbClust)

matrices_diagnostics <- function(X0, X1) {
  
  minmax_norm <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  X0_minmax <- minmax_norm(X0)
  X1_minmax <- minmax_norm(X1)
  
  cat("--- X0 Diagnostics ---\n")
  print(summary(X0))
  print(attributes(X0)$Size)
  print(attributes(X0)$method)
  print(attributes(X0)$class)
  
  cat("\n--- X1 Diagnostics ---\n")
  print(summary(X1))
  print(attributes(X1)$Size)
  print(attributes(X1)$method)
  print(attributes(X1)$class)
  
  cat("\n--- Min-Max Normalized Diagnostics ---\n")
  cat("\n--- X0 ---\n")
  print(summary(X0_minmax))
  cat("\n--- X1 ---\n")
  print(summary(X1_minmax))
  cat("\n--- abs(minmax(X0) - minmax(X1)) ---\n")
  print(summary(abs(X0_minmax - X1_minmax)))
}

base_path <- "/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering"
AOI <- "Lombardy"

setwd(base_path)

load('/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/res1911.RData')
load("data/anagrafica_comuni/LAU2021_Eurostat_Italy.RData")
DPCM <- read_excel("data/DPCM/DPCM_cleaned.xlsx")
ManualMatchLAUs <- read_excel("data/anagrafica_comuni/Manually_matched_LAUs_Oct2023.xlsx", 
                              sheet = "Modified")

ManualMatchLAUs <- ManualMatchLAUs %>%
  mutate(across(c(`Denominazione Comune NEW`), toupper),
         across(c(`Denominazione Comune NEW`), stringr::str_to_title)) %>%
  dplyr::select(LAU_Code_old = `Codice Istat Comune`,
                LAU_Code_new = `Codice Istat Comune NEW`) %>%
  distinct()



DPCM_red <- DPCM %>%
  dplyr::rename(LAU_Code = ISTAT_Code_LAU) %>%
  dplyr::select(-Name_LAU)


codes <- stringr::str_pad(string = ManualMatchLAUs$LAU_Code_old,width = 6,pad = "0")
codes_new <- stringr::str_pad(string = ManualMatchLAUs$LAU_Code_new,width = 6,pad = "0")
for (r in 1:dim(DPCM_red)[1]) {
  pos <- which(DPCM_red$LAU_Code[r] == codes)
  if (!is_empty(pos)) {
    cat(paste0("Row ",r," has a different code\n"))
    DPCM_red[r,]$LAU_Code <- codes_new[pos]
  }
}

DPCM_red <- DPCM_red %>%
  group_by(LAU_Code) %>%
  mutate(across(is.numeric, function(x) sum(x,na.rm=T))) %>%
  ungroup()

DPCM_red <- left_join(x = LAU_IT_sf, y = DPCM_red, by = c("LAU_Code"))


if (AOI == "Lombardy") {
  DPCM_red <- DPCM_red %>%
    filter(NUTS2_Name == "Lombardia")
}


CData_reconstructed <- DPCM_red %>%
  # Filter non-missing data
  filter(
    !is.na(Local_Units_Industry_Active),
    LAU_Code != "015146" # Remove Milan
  ) %>%
  # Select numeric variables
  dplyr::select(Area_km2,LAU_Code,Lon_centr_LAU,Lat_centr_LAU,
                Alt_mean,Alt_min,Alt_max,
                contains("Suspended"),contains("Active")) %>%
  select(-c(Area_km2)) %>%
  ### Standardization 
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) scale(x)[,1]))

st_geometry(CData_reconstructed) <- NULL

Xmat_nogeo_reconstructed <- as.matrix(CData_reconstructed[,-c(1:6)])
D0_reconstructed <- get_dist(Xmat_nogeo_reconstructed, method = "euclidean")


Xmat_geo_reconstructed <- as.matrix(CData_reconstructed[,2:3])
D1_reconstructed <- get_dist(Xmat_geo_reconstructed)


#### Results Diagnostics ####

## Check CData reconstruction
identical(
  CData %>% 
    select(where(is.numeric)) %>% 
    mutate(
      across(where(is.numeric), round, 9)
    ), 
  CData_reconstructed %>% 
    select(where(is.numeric)) %>% 
    mutate(
      across(where(is.numeric), round, 9))
)

## Check Xmat_nogeo and Xmat_geo reconstruction
identical(round(Xmat_nogeo, 9), round(Xmat_nogeo_reconstructed, 9))
identical(Xmat_geo, Xmat_geo_reconstructed)

## Check D0 and D1 reconstruction
identical(D0, D0_reconstructed)
matrices_diagnostics(D0, D0_reconstructed)
identical(D1, D1_reconstructed)
matrices_diagnostics(D1, D1_reconstructed)

