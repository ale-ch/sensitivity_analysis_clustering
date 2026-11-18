# source("scripts/replication_package/env_preparation.R")

##### DPCM data
# DPCM <- read_excel("data/DPCM/DPCM_cleaned.xlsx")
DPCM <- read_excel(file.path(data_path, "DPCM/DPCM_cleaned.xlsx"))

##### Municipalities registry
# ManualMatchLAUs <- read_excel("data/anagrafica_comuni/Manually_matched_LAUs_Oct2023.xlsx", sheet = "Modified")
ManualMatchLAUs <- read_excel(file.path(data_path, "anagrafica_comuni/Manually_matched_LAUs_Oct2023.xlsx"), sheet = "Modified")
ManualMatchLAUs <- ManualMatchLAUs %>%
  mutate(across(c(`Denominazione Comune NEW`), toupper),
         across(c(`Denominazione Comune NEW`), stringr::str_to_title)) %>%
  dplyr::select(LAU_Code_old = `Codice Istat Comune`,
                LAU_Code_new = `Codice Istat Comune NEW`) %>%
  distinct()

# load("../../Dati/DataManagement/Data transformed/Unnormalized/Data_CrossSection_C2.RData")
load(file.path(data_path, "anagrafica_comuni/LAU2021_Eurostat_Italy.RData"))

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


if (study_area == "Lombardia") {
  DPCM_red <- DPCM_red %>%
    filter(NUTS2_Name == "Lombardia")
}




conta_codici <- function(dati, campo, fase) {
  codici <- as.character(dati[[campo]])
  validi <- !is.na(codici) & nzchar(codici)
  frequenze <- table(codici[validi])
  data.frame(fase = fase, righe = length(codici), codici_mancanti = sum(!validi),
             unique_municipalities = length(frequenze), duplicated_municipalities = sum(frequenze > 1L),
             righe_in_eccesso = sum(pmax(as.integer(frequenze) - 1L, 0L)))
}
stopifnot(!anyDuplicated(LAU_IT_sf$LAU_Code),
          !anyNA(DPCM_red$LAU_Code), all(nzchar(as.character(DPCM_red$LAU_Code))))
colonne_verifica_copie <- c(clustering_variables, "Area_km2", "Lon_centr_LAU", "Lat_centr_LAU")
duplicate_value_check <- sf::st_drop_geometry(DPCM_red) %>%
  group_by(LAU_Code) %>%
  summarise(Comune = first(LAU_Name), righe = n(),
            across(all_of(colonne_verifica_copie), function(x) n_distinct(x),
                   .names = "n_valori_{.col}"), .groups = "drop") %>%
  filter(righe > 1L)
if (nrow(duplicate_value_check)) {
  valori_diversi <- as.matrix(dplyr::select(duplicate_value_check, starts_with("n_valori_")))
  if (any(valori_diversi > 1L)) {
    stop("Copie dello stesso LAU_Code con valori diversi: verificare i dati prima del filtro.")
  }
}
municipalities_duplicated_before_filter <- duplicate_value_check %>%
  dplyr::select(LAU_Code, Comune, righe) %>%
  mutate(righe_rimosse = righe - 1L)
rows_removed_as_duplicates <- sf::st_drop_geometry(DPCM_red) %>%
  mutate(id_riga_prima_filtro = row_number()) %>%
  filter(duplicated(LAU_Code)) %>%
  dplyr::select(id_riga_prima_filtro, LAU_Code, LAU_Name, all_of(clustering_variables))
Audit_prima_filtro <- conta_codici(DPCM_red, "LAU_Code", "Area prima del filtro duplicati")
DPCM_red <- DPCM_red %>% ungroup() %>% filter(!duplicated(LAU_Code))
Audit_dopo_filtro <- conta_codici(DPCM_red, "LAU_Code", "Area dopo il filtro duplicati")
duplicate_filter_audit <- bind_rows(Audit_prima_filtro, Audit_dopo_filtro)
stopifnot(Audit_dopo_filtro$duplicated_municipalities == 0L,
          Audit_dopo_filtro$righe == Audit_prima_filtro$unique_municipalities,
          Audit_prima_filtro$righe - Audit_dopo_filtro$righe == nrow(rows_removed_as_duplicates),
          sum(municipalities_duplicated_before_filter$righe_rimosse) == nrow(rows_removed_as_duplicates))
write.csv(duplicate_filter_audit, file.path(output_dir, "Audit_filtro_duplicati.csv"), row.names = FALSE)
write.csv(municipalities_duplicated_before_filter, file.path(output_dir, "Comuni_duplicati_prima_filtro.csv"), row.names = FALSE)
write.csv(rows_removed_as_duplicates, file.path(output_dir, "Righe_rimosse_duplicati.csv"), row.names = FALSE)
write.csv(duplicate_value_check, file.path(output_dir, "Controllo_duplicati_valori.csv"), row.names = FALSE)
knitr::kable(duplicate_filter_audit, row.names = FALSE)
knitr::kable(municipalities_duplicated_before_filter, row.names = FALSE)




DPCM_prima_Milano <- DPCM_red
if (study_area == "Lombardia" && !include_milan) {
  DPCM_red <- DPCM_red %>% filter(LAU_Code != "015146")
}



CData <- DPCM_red %>%
  # Filter non-missing data
  filter(!is.na(Local_Units_Industry_Active)) %>%
  # Select numeric variables
  dplyr::select(Area_km2,LAU_Code,Lon_centr_LAU,Lat_centr_LAU,
                Alt_mean,Alt_min,Alt_max,
                all_of(clustering_variables)) %>%
  mutate(across(c(is.numeric,
                  -contains("Alt"),
                  -contains("Lon"),-contains("Lat")), function(x) x/Area_km2)) %>%
  select(-c(Area_km2)) %>%
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))



st_geometry(CData) <- NULL
Xmat_geo <- as.matrix(CData[, 2:3])
Xmat_nogeo <- as.matrix(CData[, -c(1:6)])
DPCM_red2 <- DPCM_red %>% filter(!is.na(Local_Units_Industry_Active))

DMat <- get_dist(Xmat_nogeo)
DPCM_red2$AvgDist_D0 <- rowMeans(as.matrix(DMat))
D0 <- DMat
n <- nrow(Xmat_nogeo)
stopifnot(ncol(Xmat_geo) == 2L, all(is.finite(Xmat_geo)),
          all(abs(Xmat_geo[, 1L]) <= 180), all(abs(Xmat_geo[, 2L]) <= 90))
distanza_WGS84 <- function(p1, p2) {
  geosphere::distGeo(p1, p2, a = 6378137, f = 1 / 298.257223563)
}
DMat_geo <- as.dist(geosphere::distm(Xmat_geo, fun = distanza_WGS84))
DPCM_red2$AvgDist_D1 <- rowMeans(as.matrix(DMat_geo))
D1 <- DMat_geo
D0_std <- D0 / max(D0)
D1_std <- D1 / max(D1)
wt <- rep(1 / n, n)

stopifnot(n > max(15L, K_candidati), ncol(Xmat_nogeo) == 12L,
          identical(colnames(Xmat_nogeo), clustering_variables),
          !anyDuplicated(DPCM_red2$LAU_Code),
          nrow(DPCM_red2) == n, attr(D0, "Size") == n, attr(D1, "Size") == n,
          identical(as.character(CData$LAU_Code), as.character(DPCM_red2$LAU_Code)),
          all(is.finite(D0)), all(is.finite(D1)), max(D0) > 0, max(D1) > 0)

# A row ID added AFTER preprocessing prevents ambiguous joins when duplicates are present.
# It is not included in CData, Xmat_nogeo, Xmat_geo, or the distance matrices.
DPCM_red2$id_osservazione <- seq_len(n)
audit_dimensioni <- data.frame(
  Osservazioni_nella_stima = n,
  Codici_comunali_distinti = dplyr::n_distinct(DPCM_red2$LAU_Code),
  Righe_duplicate_per_codice = n - dplyr::n_distinct(DPCM_red2$LAU_Code))
knitr::kable(audit_dimensioni)