# Assign by row ID: joining by municipality code alone would multiply duplicate rows.


# Exact assignment across K clusters (up to 15) using dynamic programming.
# This function is used only for map colors and does not change any clustering.
abbina_colori <- function(riferimento, altra, k = 8L) {
  tab <- table(factor(riferimento, levels = seq_len(k)),
               factor(altra, levels = seq_len(k)))
  ultimo <- 2^k - 1L
  score <- rep(-Inf, ultimo + 1L)
  score[1L] <- 0
  precedente <- scelto <- integer(ultimo + 1L)
  for (mask in 0:(ultimo - 1L)) {
    if (!is.finite(score[mask + 1L])) next
    usati <- as.integer(intToBits(mask))[seq_len(k)]
    riga <- sum(usati) + 1L
    for (col in which(usati == 0L)) {
      nuovo <- bitwOr(as.integer(mask), bitwShiftL(1L, col - 1L))
      valore <- score[mask + 1L] + tab[riga, col]
      if (valore > score[nuovo + 1L]) {
        score[nuovo + 1L] <- valore
        precedente[nuovo + 1L] <- mask
        scelto[nuovo + 1L] <- col
      }
    }
  }
  mappa <- integer(k)
  mask <- ultimo
  for (riga in k:1L) {
    mappa[scelto[mask + 1L]] <- riga
    mask <- precedente[mask + 1L]
  }
  list(etichette = unname(mappa[altra]), mappa = mappa)
}

partizioni_colori <- partizioni[c("id_osservazione", "LAU_Code")]
abbinamenti <- list()
for (i in seq_len(nrow(casi))) {
  nm <- casi$caso[i]
  if (nm == "optimal_alpha") {
    ab <- list(etichette = partizioni[[nm]], mappa = seq_len(K))
  } else ab <- abbina_colori(partizioni$optimal_alpha, partizioni[[nm]], K)
  partizioni_colori[[nm]] <- ab$etichette
  abbinamenti[[nm]] <- data.frame(caso = nm,
                                  etichetta_originale = seq_len(K), etichetta_colore = ab$mappa)
  # Verifica che la rietichettatura sia una permutazione uno a uno.
  stopifnot(setequal(ab$mappa, seq_len(K)))
}
abbinamenti <- dplyr::bind_rows(abbinamenti)

# Assign by row ID: joining by municipality code alone would multiply duplicate rows.
Data_clusters <- DPCM_red2
for (nm in casi$caso) {
  Data_clusters[[paste0(nm, "_originale")]] <- partizioni[[nm]]
  Data_clusters[[nm]] <- partizioni_colori[[nm]]
}
if (study_area == "Lombardia" && !include_milan) {
  Milano <- DPCM_prima_Milano %>% filter(LAU_Code == "015146")
  stopifnot(nrow(Milano) >= 1L)
  Milano$id_osservazione <- n + seq_len(nrow(Milano))
  for (nm in casi$caso) {
    Milano[[paste0(nm, "_originale")]] <- K + 1L
    Milano[[nm]] <- K + 1L
  }
  Milano$D0_P_Kfinale <- Milano$D0D1_P_Kfinale <-
    Milano$D0D1_P_Kfinale_alpha010 <- K + 1L
  Milano$D0_P8 <- Milano$D0D1_P8 <- Milano$D0D1_P8_alpha010 <- 9L
  Data_clusters <- dplyr::bind_rows(Data_clusters, Milano)
}


milano_fuori_stima <- study_area == "Lombardia" && !include_milan
n_colori <- K + as.integer(milano_fuori_stima)
available_palettes <- list(
  sobria = c("#3475AB", "#E18B3C", "#3E9378", "#A56DA5",
             "#B9474E", "#60A4AD", "#A28E3D", "#5F6572",
             "#7B899E", "#B2A181", "#7F8B62", "#977D87", "#759593", "#A37958", "#5F7992"),
  pastello = c("#8BA9C4", "#D4AA84", "#9BB7A6", "#BAA7C4",
               "#C99196", "#94B9BC", "#C3B984", "#A6ABB3",
               "#A3B2C6", "#C8B99B", "#AFB995", "#C3A8B2", "#AAC5C2", "#C6AB91", "#94AFC4"))
multi_k_palette <- setNames(available_palettes[[params$map_palette]], as.character(1:15))
cluster_palette <- unname(multi_k_palette[as.character(seq_len(K))])
if (milano_fuori_stima) cluster_palette <- c(cluster_palette, "#14212B")
cluster_palette <- setNames(cluster_palette,
                            as.character(seq_len(n_colori)))
legenda_cluster <- paste("C", seq_len(n_colori), sep = "")
if (milano_fuori_stima) legenda_cluster[K + 1L] <- paste0("C", K + 1L, ": Milano fuori stima")
nomi_pannelli <- setNames(casi$pannello, casi$caso)


casi <- data.frame(
  caso = c("alpha0", "optimal_alpha", "alpha010", "alpha1"),
  alpha = c(0, optimal_alpha, params$manual_alpha , 1),
  stringsAsFactors = FALSE)
casi$pannello <- c("alpha = 0",
                   paste0("alpha ", if (params$alpha_selection_mode == "automatic") "optimal" else "scelto",
                          " a K0=", K0, ": ", format(optimal_alpha, digits = 5)),
                   "alpha = 0.10", "alpha = 1")
casi$K <- K
casi$K0_scelta_alpha <- K0
casi$criterio_alpha <- params$alpha_criterion

alberi <- vector("list", nrow(casi))
names(alberi) <- casi$caso
partizioni <- data.frame(id_osservazione = seq_len(n),
                         LAU_Code = as.character(DPCM_red2$LAU_Code))
for (i in seq_len(nrow(casi))) {
  # Riutilizza l'albero gia stimato al passaggio 3, quando l'alpha coincide.
  alberi[[i]] <- if (abs(casi$alpha[i] - optimal_alpha) < 1e-15) {
    tree_alpha_scelto
  } else hclustgeo(D0_std, D1_std, alpha = casi$alpha[i], wt = rep(1/n, n))
  partizioni[[casi$caso[i]]] <- cutree(alberi[[i]], K)
  stopifnot(length(unique(partizioni[[casi$caso[i]]])) == K)
}
treeD0D1 <- alberi$optimal_alpha
treeD0D1_alpha010 <- alberi$alpha010
treeD0D1_alpha0 <- alberi$alpha0
treeD0D1_alpha1 <- alberi$alpha1
# Campi principali: K finale; campi storici _P8: realmente otto gruppi.
DPCM_red2$D0_P_Kfinale <- partizioni$alpha0
DPCM_red2$D0D1_P_Kfinale <- partizioni$optimal_alpha
DPCM_red2$D0D1_P_Kfinale_alpha010 <- partizioni$alpha010
DPCM_red2$D0_P8 <- cutree(tree_D0, 8L)
DPCM_red2$D0D1_P8 <- cutree(treeD0D1, 8L)
DPCM_red2$D0D1_P8_alpha010 <- cutree(treeD0D1_alpha010, 8L)
stopifnot(length(unique(DPCM_red2$D0_P8)) == 8L,
          length(unique(DPCM_red2$D0D1_P8)) == 8L)



Mappe <- Data_clusters %>%
  dplyr::select(id_osservazione, LAU_Code, all_of(casi$caso)) %>%
  tidyr::pivot_longer(cols = all_of(casi$caso),
                      names_to = "caso", values_to = "Cluster") %>%
  mutate(Pannello = factor(unname(nomi_pannelli[caso]), levels = casi$pannello),
         Cluster = factor(Cluster, levels = seq_len(n_colori)))

# Main 4-Panel Map (Alpha 0, Optimal, 0.10, 1)
p_mappe <- ggplot(Mappe) +
  geom_sf(aes(fill = Cluster), linewidth = 0.045, colour = "#F5F3EF") +
  facet_wrap(~ Pannello, ncol = 2) +
  scale_fill_manual(values = cluster_palette, labels = legenda_cluster, drop = FALSE, name = "Cluster") +
  labs(title = paste("Ward spaziale -", params$area),
       # subtitle = paste(sottotitolo, "-", descrizione_criterio_alpha), x = "", y = "",
       caption = "Colori abbinati alla partizione dell'alpha scelto; etichette native conservate.") +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom", legend.key.height = grid::unit(0.4, "cm"),
        strip.text = element_text(face = "bold", colour = "#183746", size = 12),
        plot.title = element_text(face = "bold", colour = "#183746", size = 15),
        plot.subtitle = element_text(colour = "#596772"),
        plot.caption = element_text(colour = "#687580", hjust = 0),
        plot.margin = margin(12, 12, 12, 12)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(output_dir, "Mappe_quattro_alpha.png"), p_mappe, width = 13, height = 11, dpi = 200, bg = "white")

# 28-Panel Multi-K Map (K=2 to K=15)

########

K_confronto <- 2:15
Casi_multik <- data.frame(K = rep(K_confronto, each = 2L),
                          caso = rep(c("alpha0", "optimal_alpha"), length(K_confronto)),
                          alpha = rep(c(0, optimal_alpha), length(K_confronto)),
                          stringsAsFactors = FALSE)
Casi_multik$campo <- paste0(Casi_multik$caso, "_K", Casi_multik$K)
Casi_multik$pannello <- paste0("K = ", Casi_multik$K, " | ",
                               ifelse(Casi_multik$caso == "alpha0", "alpha = 0",
                                      paste0("alpha ottimo = ", format(optimal_alpha, digits = 6))))
Partizioni_multik_confronto <- partizioni[c("id_osservazione", "LAU_Code")]
Colori_multik_confronto <- Partizioni_multik_confronto
Abbinamento_colori_multik <- Dimensioni_multik_lista <- ARI_multik_lista <- list()
riferimento_15 <- cutree(treeD0D1, k = 15L)
for (kk in K_confronto) {
  g0 <- cutree(treeD0D1_alpha0, k = kk)
  gs <- cutree(treeD0D1, k = kk)
  rami <- vapply(seq_len(kk), function(gr) {
    as.integer(which.max(table(factor(riferimento_15[gs == gr], levels = 1:15))))
  }, integer(1))
  stopifnot(!anyDuplicated(rami))
  ab0 <- abbina_colori(gs, g0, kk)
  for (cc in c("alpha0", "optimal_alpha")) {
    campo <- paste0(cc, "_K", kk)
    gruppi <- if (cc == "alpha0") g0 else gs
    mappa_rami <- if (cc == "alpha0") rami[ab0$mappa] else rami
    Partizioni_multik_confronto[[campo]] <- gruppi
    Colori_multik_confronto[[campo]] <- unname(mappa_rami[gruppi])
    Abbinamento_colori_multik[[length(Abbinamento_colori_multik) + 1L]] <- data.frame(
      K = kk, caso = cc, cluster_nativo = seq_len(kk), ramo_colore = mappa_rami)
    Dimensioni_multik_lista[[length(Dimensioni_multik_lista) + 1L]] <- data.frame(
      caso = cc, alpha = if (cc == "alpha0") 0 else optimal_alpha, K = kk,
      cluster = seq_len(kk),
      n_comuni = as.integer(table(factor(gruppi, levels = seq_len(kk)))),
      ramo_colore = mappa_rami)
  }
  ARI_multik_lista[[length(ARI_multik_lista) + 1L]] <- data.frame(
    K = kk, alpha_0 = 0, optimal_alpha = optimal_alpha, n_comuni = n,
    ARI = mclust::adjustedRandIndex(g0, gs))
}
Abbinamento_colori_multik <- bind_rows(Abbinamento_colori_multik)
Abbinamento_colori_multik$colore_hex <- unname(multi_k_palette[
  as.character(Abbinamento_colori_multik$ramo_colore)])
Dimensioni_multik <- bind_rows(Dimensioni_multik_lista)
ARI_multik <- bind_rows(ARI_multik_lista)
stopifnot(nrow(Casi_multik) == 28L, all(Dimensioni_multik$n_comuni > 0))



dati_mappa <- Data_clusters %>% dplyr::select(id_osservazione, LAU_Code)
stopifnot(all(dati_mappa$id_osservazione[seq_len(n)] == seq_len(n)))
n_extra <- nrow(dati_mappa) - n
stopifnot(n_extra == (if (milano_fuori_stima) 1L else 0L))
for (campo in Casi_multik$campo) {
  dati_mappa[[campo]] <- c(Colori_multik_confronto[[campo]], rep(16L, n_extra))
}
k_pannelli_multik <- setNames(Casi_multik$K, Casi_multik$campo)
alpha_pannelli_multik <- setNames(Casi_multik$caso, Casi_multik$campo)
n_color <- 15L + as.integer(milano_fuori_stima)
colori <- unname(multi_k_palette[as.character(1:15)])
if (milano_fuori_stima) colori <- c(colori, "#14212B")
names(colori) <- as.character(seq_len(n_color))
legenda <- paste0("R", 1:15)
if (milano_fuori_stima) legenda <- c(legenda, "Milano fuori stima")
dati_mappa <- dati_mappa %>%
  pivot_longer(cols = all_of(Casi_multik$campo), names_to = "campo", values_to = "ramo") %>%
  mutate(K_riga = factor(unname(k_pannelli_multik[campo]), levels = K_confronto,
                         labels = paste0("K = ", K_confronto)),
         Alpha_colonna = factor(unname(alpha_pannelli_multik[campo]),
                                levels = c("alpha0", "optimal_alpha"), labels = c("alpha = 0", "alpha = ottimo")),
         ramo = factor(ramo, levels = seq_len(n_color)))
p <- ggplot(dati_mappa) +
  geom_sf(aes(fill = ramo), linewidth = 0.025, colour = "#F5F3EF") +
  facet_grid(K_riga ~ Alpha_colonna, switch = "y") +
  scale_fill_manual(values = colori, labels = legenda, drop = FALSE,
                    name = "Rami di riferimento a K=15") +
  labs(title = paste("Clustering da K=2 a K=15 -", params$area),
       subtitle = paste("Coppie alpha=0 / alpha ottimo; alpha ottimo =", format(optimal_alpha, digits = 6)),
       caption = "Colori riferiti all'albero alpha ottimo; ID nativi dei cluster esportati separatamente.",
       x = "", y = "") +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10, colour = "#183746"),
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_rect(fill = "#F0F2F3", colour = NA),
        plot.title = element_text(face = "bold", size = 16, colour = "#183746"),
        plot.caption = element_text(hjust = 0, colour = "#687580"),
        plot.margin = margin(12, 12, 12, 12)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
ggsave(file.path(output_dir, "Mappe_alpha0_alphaottimo_K2_K15.png"),
       p, width = 12, height = 48, dpi = 170, bg = "white")
ggsave(file.path(output_dir, "Mappe_alpha0_alphaottimo_K2_K15.pdf"),
       p, width = 12, height = 48)



# Normalized Variables Map
variabili_sei <- c("Local_Units_Industry_Suspended", "Value_Added_Industry_Suspended",
                   "Employed_Industry_Suspended", "Local_Units_Services_Suspended",
                   "Value_Added_Services_Suspended", "Employed_Services_Suspended")
etichette_sei <- c("Industry: local units suspended", "Industry: value added of local units suspended",
                   "Industry: employed persons suspended", "Services: local units suspended",
                   "Services: value added of local units suspended", "Services: employed persons suspended")
colori_continui <- c("#EDF2F3", "#8FB4C1", "#2D607A")
tema_originale <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", colour = "#183746"),
        plot.subtitle = element_text(colour = "#596772"),
        legend.position = "bottom")
dati_norm <- DPCM_red2 %>% dplyr::select(all_of(variabili_sei)) %>%
  mutate(across(is.numeric, function(x) x / max(x, na.rm = TRUE))) %>%
  pivot_longer(cols = all_of(variabili_sei), names_to = "Vars", values_to = "Values") %>%
  mutate(Vars = factor(Vars, levels = variabili_sei, labels = etichette_sei, ordered = TRUE))
p_norm <- ggplot(dati_norm) + geom_sf(aes(fill = Values), colour = NA) +
  scale_fill_gradientn(colours = colori_continui) + facet_wrap(~ Vars, nrow = 2, ncol = 3) +
  labs(title = "Normalized DPCM-related variables by industry and services",
       subtitle = paste(params$area, "- normalizzazione per il massimo solo nel grafico"), x = "", y = "") + tema_originale

ggsave(file.path(output_dir, "NormVariablesD0.png"), p_norm, width = 1800/150, height = 1200/150, dpi = 150, bg = "white")