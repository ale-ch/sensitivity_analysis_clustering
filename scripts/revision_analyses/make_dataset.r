library(sf)
library(tidyverse)
library(readxl)
library(pracma)
library(cluster)
library(stats)
library(fpc)
library(factoextra) 
library(NbClust)
library(ClustGeo)
library(viridis)
library(ggpubr)



#### Prepare data ####
base_path <- '/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering'
setwd(base_path)

output_dir <- "/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/output"

'%notin%' <- Negate('%in%')

# AOI <- "Italy"
AOI <- "Lombardy"

##### DPCM data
DPCM <- read_excel("./data/DPCM/DPCM_cleaned.xlsx")

##### Municipalities registry
ManualMatchLAUs <- read_excel("./data/anagrafica_comuni/Manually_matched_LAUs_Oct2023.xlsx", sheet = "Modified")
ManualMatchLAUs <- ManualMatchLAUs %>%
  mutate(across(c(`Denominazione Comune NEW`), toupper),
         across(c(`Denominazione Comune NEW`), stringr::str_to_title)) %>%
  dplyr::select(LAU_Code_old = `Codice Istat Comune`,
                LAU_Code_new = `Codice Istat Comune NEW`) %>%
  distinct()

# load("../../Dati/DataManagement/Data transformed/Unnormalized/Data_CrossSection_C2.RData")
load("./data/anagrafica_comuni/LAU2021_Eurostat_Italy.RData")

DPCM_red <- DPCM %>%
  dplyr::rename(LAU_Code = ISTAT_Code_LAU) %>%
  dplyr::select(-Name_LAU)
# DPCM_red <- data %>%
#   dplyr::rename(LAU_Code = AUX_IstatMunicipalityCode)

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


CData <- DPCM_red %>%
  # Filter Lombardy
  filter(NUTS2_Name == "Lombardia") %>%
  # Filter non-missing data
  filter(!is.na(Local_Units_Industry_Active)) %>%
  # Select numeric variables
  dplyr::select(Area_km2,LAU_Code,Lon_centr_LAU,Lat_centr_LAU,
                Alt_mean,Alt_min,Alt_max,
                contains("Suspended"),contains("Active")) %>%
  mutate(across(c(is.numeric,
                  -contains("Alt"),
                  -contains("Lon"),-contains("Lat")), function(x) x/Area_km2)) %>%
  select(-c(Area_km2)) %>%
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))

st_geometry(CData) <- NULL
Xmat_geo <- as.matrix(CData[,2:3])
Xmat_nogeo <- as.matrix(CData[,-c(1:6)])

D0 <- get_dist(Xmat_nogeo)
D1 <- get_dist(Xmat_geo)
n <- nrow(D0)



rm(list=setdiff(ls(), c("CData", "D0", "D1", "n")))

