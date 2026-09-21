#####################################
########## DPCM management ##########
#####################################

##### Date: 02/11/2023
##### Source: https://www.istat.it/it/archivio/241341

##### Set working directory
setwd("H:/Il mio Drive/FEEM_COVID19_Lombardia/Dati/DataManagement/Original data/comune/DPCM")

##### Packages
library(readxl)
library(tidyverse)

##### Load data and stacking
DPCM_ind_act <- read_excel("dati_comunali_2017_DPCM_covid19.xlsx", sheet = "Settori attivi_industria")
DPCM_ind_act <- DPCM_ind_act %>%
  mutate(Settore = "Industry", Stato = "Active")
DPCM_ind_sus <- read_excel("dati_comunali_2017_DPCM_covid19.xlsx", sheet = "settori sospesi industria")
DPCM_ind_sus <- DPCM_ind_sus %>%
  mutate(Settore = "Industry", Stato = "Suspended")
DPCM_ser_act <- read_excel("dati_comunali_2017_DPCM_covid19.xlsx", sheet = "settori attivi servizi")
DPCM_ser_act <- DPCM_ser_act %>%
  mutate(Settore = "Services", Stato = "Active")
DPCM_ser_sus <- read_excel("dati_comunali_2017_DPCM_covid19.xlsx", sheet = "settori sospesi servizi")
DPCM_ser_sus <- DPCM_ser_sus %>%
  mutate(Settore = "Services", Stato = "Suspended")
DPCM <- bind_rows(DPCM_ind_act,DPCM_ind_sus,DPCM_ser_act,DPCM_ser_sus)

##### Management
DPCM <- DPCM %>%
  select(ISTAT_Code_LAU = codice_comune, Name_LAU = Denominazione_comune,
         Sector = Settore, Status = Stato,
         Local_Units = Unita_locali, Value_Added = `Valore_aggiunto\r\n(valori in euro)`,
         Revenues = `Fatturato\r\n(valori in euro)`,
         Workers = `Numero Addetti`, Employed = `Numero Dipendenti`) %>%
  mutate(across(c(Local_Units,Value_Added,Revenues,Workers,Employed),
         function(x) case_when(x == "*" ~ NA,
                               TRUE ~ x))) %>%
  mutate(across(c(Local_Units,Value_Added,Revenues,Workers,Employed),
                function(x) as.numeric(x)))


DPCM_Tot <- DPCM %>%
  group_by(ISTAT_Code_LAU,Name_LAU,Sector) %>%
  summarise(across(c(Local_Units,Value_Added,Revenues,Workers,Employed),
            function(x) sum(x, na.rm=T))) %>%
  ungroup() %>%
  mutate(Status = "Total", .after = "Sector")

DPCM <- bind_rows(DPCM,DPCM_Tot)

DPCM <- DPCM %>%
  pivot_wider(names_from = c("Sector","Status"),
              values_from = c("Local_Units","Value_Added","Revenues","Workers","Employed"))

# DPCM %>%
#   ggplot(mapping = aes(x = Local_Units_Industry_Active, y = Local_Units_Industry_Suspended)) + 
#   geom_point(size = 3)
# 
# # 6904 comuni completi
# DPCM_complete <- DPCM[complete.cases(DPCM),]
# apply(X = DPCM, MARGIN = 2, FUN = function(x) sum(is.na(x)))

DPCM <- DPCM %>%
  mutate(across(is.numeric, function(x) replace_na(x,0))) %>%
  arrange(ISTAT_Code_LAU,Name_LAU)

### Export dataset
##### Creo file Excel
wb <- createWorkbook("../DPCM_cleaned")
##### Save Excel sheet
addWorksheet(wb,"DPCM COVID-19")
writeData(wb, sheet = "DPCM COVID-19", DPCM, colNames = T)
##### Salvatggio file Excel
saveWorkbook(wb,"../DPCM_cleaned.xlsx",overwrite = T)


