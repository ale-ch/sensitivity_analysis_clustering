##############################################
########## Spatial cluster analysis ##########
##############################################

setwd("/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering")

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

'%notin%' <- Negate('%in%')

# AOI <- "Italy"
AOI <- "Lombardy"

source("scripts/original/AuxFuns_NbClust.R", encoding = 'UTF-8')

##### DPCM data
DPCM <- read_excel("data/DPCM/DPCM_cleaned.xlsx")

##### Municipalities registry
ManualMatchLAUs <- read_excel("data/anagrafica_comuni/Manually_matched_LAUs_Oct2023.xlsx", sheet = "Modified")
ManualMatchLAUs <- ManualMatchLAUs %>%
  mutate(across(c(`Denominazione Comune NEW`), toupper),
         across(c(`Denominazione Comune NEW`), stringr::str_to_title)) %>%
  dplyr::select(LAU_Code_old = `Codice Istat Comune`,
                LAU_Code_new = `Codice Istat Comune NEW`) %>%
  distinct()

# load("../../Dati/DataManagement/Data transformed/Unnormalized/Data_CrossSection_C2.RData")
load("data/anagrafica_comuni/LAU2021_Eurostat_Italy.RData")

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


DPCM_red2 <- DPCM_red %>%
  filter(!is.na(Local_Units_Industry_Active))

# png("SuspendedDPCM.png", width = 1800, height = 1200, res = 150)
DPCM_red %>%
  mutate('Share of suspended local units in industry' = Local_Units_Industry_Suspended/Local_Units_Industry_Total*100,
         'Share of suspended local units in services' = Local_Units_Services_Suspended/Local_Units_Services_Total*100) %>%
  select('Share of suspended local units in industry','Share of suspended local units in services') %>%
  pivot_longer(cols = c('Share of suspended local units in industry',
                        'Share of suspended local units in services'),
               names_to = "Var", values_to = "Values") %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Values)) + 
  scale_fill_viridis_b(option = "plasma") + 
  facet_wrap(~ Var) + 
  theme_bw()

DPCM_red %>%
  #filter(NUTS2_Name %in% c("Lombardia")) %>%
  select(Local_Units_Industry_Suspended,Local_Units_Services_Suspended) %>%
  pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended),
               names_to = "Var", values_to = "Values") %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Values)) + 
  scale_fill_viridis_b(option = "plasma",n.breaks = 100) + 
  facet_wrap(~ Var) + 
  theme_bw()



DPCM_red2 %>%
  filter(NUTS2_Name == "Lombardia", !(LAU_Code %in% c("015146"))) %>% 
  select(Local_Units_Industry_Suspended,Value_Added_Industry_Suspended,Employed_Industry_Suspended,
         Local_Units_Services_Suspended,Value_Added_Services_Suspended,Employed_Services_Suspended) %>%
  mutate(across(is.numeric, function(x) x/max(x,na.rm = T))) %>%
  pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended,
                        Value_Added_Industry_Suspended,Value_Added_Services_Suspended,
                        Employed_Industry_Suspended,Employed_Services_Suspended),
               names_to = "Vars", values_to = "Values") %>%
  mutate(Vars = case_when(
    Vars %in% c("Local_Units_Industry_Suspended") ~ "Industry: local units suspended",
    Vars %in% c("Value_Added_Industry_Suspended") ~ "Industry: value adde of local units suspended",
    Vars %in% c("Employed_Industry_Suspended") ~ "Industry: employed persons suspended",
    Vars %in% c("Local_Units_Services_Suspended") ~ "Services: local units suspended",
    Vars %in% c("Value_Added_Services_Suspended") ~ "Services: value adde of local units suspended",
    Vars %in% c("Employed_Services_Suspended") ~ "Services: employed persons suspended",
    TRUE ~ NA_character_
  ),
  Vars = factor(
    Vars,
    levels = c("Industry: local units suspended","Industry: value adde of local units suspended",
               "Industry: employed persons suspended","Services: local units suspended",
               "Services: value adde of local units suspended","Services: employed persons suspended"),
    labels = c("Industry: local units suspended","Industry: value adde of local units suspended",
               "Industry: employed persons suspended","Services: local units suspended",
               "Services: value adde of local units suspended","Services: employed persons suspended"),
    ordered = T)) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Values)) + 
  scale_fill_viridis_c(option = "plasma") + 
  facet_wrap(~ Vars, nrow = 2, ncol = 3) + 
  labs(title = "Normalized DPCM-related variables by industry and services",
       subtitle = "Milan is excluded due to its uniqueness (out-of-scale for every considered variable)",
       y = "", x = "") + 
  theme_bw() + 
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))
