k_diagnostics <- dplyr::bind_rows(diagnostica_K0, final_k_diagnostic)
tagli_selezionati <- k_diagnostics[k_diagnostics$selezionato, c("fase", "K")]
k_selection_plot <- ggplot(k_diagnostics, aes(K, W_norm)) +
  geom_segment(aes(x = 1, xend = max(K_candidati), y = 1, yend = 0),
               data = tagli_selezionati, inherit.aes = FALSE,
               linetype = 3, colour = "#9AA3A9") +
  geom_line(colour = "#365F7C", linewidth = 0.7) +
  geom_point(aes(shape = candidato), colour = "#365F7C", size = 1.8) +
  geom_vline(data = tagli_selezionati, aes(xintercept = K), colour = "#B46C45", linetype = 2) +
  facet_wrap(~ fase, ncol = 2) +
  scale_x_continuous(breaks = seq_len(max(K_candidati))) +
  scale_shape_manual(values = c(1, 16), guide = "none") +
  labs(title = paste("Selezione del numero di gruppi -", params$area),
       subtitle = "Linea tratteggiata: K usato; corda punteggiata: riferimento del gomito",
       x = "K", y = "Inerzia interna riscalata fra gli estremi") +
  theme_minimal(base_size = 11)
k_selection_plot
ggsave(file.path(output_dir, "Selezione_K_gomiti.png"), k_selection_plot,
       width = 12, height = 5, dpi = 200, bg = "white")

dati_alpha_plot <- curva_alpha_K0 %>%
  dplyr::select(alpha, Q0norm, Q1norm, gap_Qnorm) %>%
  pivot_longer(cols = -alpha, names_to = "statistica", values_to = "valore") %>%
  mutate(pannello = if_else(statistica == "gap_Qnorm",
                            "Differenza assoluta: criterio da minimizzare", "Quote di inerzia spiegata normalizzate"))
alpha_selection_plot <- ggplot(dati_alpha_plot, aes(alpha, valore, colour = statistica)) +
  geom_line(linewidth = 0.7) + facet_wrap(~ pannello, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = optimal_alpha, colour = "#5F6572", linetype = 2) +
  scale_colour_manual(values = c(Q0norm = "#365F7C", Q1norm = "#B46C45", gap_Qnorm = "#527D67")) +
  labs(title = paste("Alpha selezionato a K0 =", K0, "-", params$area),
       subtitle = paste(alpha_criterion_description,
                        "- alpha automatico =", format(automatic_alpha, digits = 6),
                        "; alpha usato =", format(optimal_alpha, digits = 6)),
       x = "alpha", y = NULL, colour = NULL) + theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
alpha_selection_plot
ggsave(file.path(output_dir, "Selezione_alpha_Qnorm.png"), alpha_selection_plot,
       width = 12, height = 5, dpi = 200, bg = "white")



plot_selection_dendrograms <- function() {
  op <- par(mfrow = c(1, 2), mar = c(3, 4, 4, 1))
  on.exit(par(op))
  plot(tree_D0, labels = FALSE, hang = -1, xlab = "", ylab = "Aggregation height",
       main = paste("D0: K0 =", K0), sub = params$area)
  rect.hclust(tree_D0, k = K0, border = "#647B8F")
  plot(tree_alpha_scelto, labels = FALSE, hang = -1, xlab = "", ylab = "Aggregation height",
       main = paste("Misto: K finale =", K), sub = paste("alpha =", optimal_alpha))
  rect.hclust(tree_alpha_scelto, k = K, border = "#647B8F")
}
plot_selection_dendrograms()
png(file.path(output_dir, "Selezione_K_dendrogrammi.png"), width = 2400, height = 1200, res = 200)
tryCatch(plot_selection_dendrograms(), finally = dev.off())




casi <- data.frame(
  caso = c("alpha0", "optimal_alpha", "alpha010", "alpha1"),
  alpha = c(0, optimal_alpha, alpha_fisso, 1),
  stringsAsFactors = FALSE)
casi$pannello <- c("alpha = 0",
                   paste0("alpha ", if (params$alpha_selection_mode == "automatic") "ottimo" else "scelto",
                          " a K0=", K0, ": ", format(optimal_alpha, digits = 5)),
                   "alpha = 0.10", "alpha = 1")
casi$K <- K
casi$k0_alpha_selection <- K0
casi$alpha_criterion <- params$alpha_criterion

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
# Main fields use final K; historical _P8 fields always refer to eight clusters.
DPCM_red2$D0_P_Kfinale <- partizioni$alpha0
DPCM_red2$D0D1_P_Kfinale <- partizioni$optimal_alpha
DPCM_red2$D0D1_P_Kfinale_alpha010 <- partizioni$alpha010
DPCM_red2$D0_P8 <- cutree(tree_D0, 8L)
DPCM_red2$D0D1_P8 <- cutree(treeD0D1, 8L)
DPCM_red2$D0D1_P8_alpha010 <- cutree(treeD0D1_alpha010, 8L)
stopifnot(length(unique(DPCM_red2$D0_P8)) == 8L,
          length(unique(DPCM_red2$D0D1_P8)) == 8L)



T0 <- ClustGeo::inertdiss(D0, wt = wt)
T1 <- ClustGeo::inertdiss(D1, wt = wt)
indici <- vector("list", nrow(casi))
for (i in seq_len(nrow(casi))) {
  g <- partizioni[[casi$caso[i]]]
  s <- fpc::cluster.stats(DMat, g)
  q0 <- 1 - ClustGeo::withindiss(D0, g, wt = wt) / T0
  q1 <- 1 - ClustGeo::withindiss(D1, g, wt = wt) / T1
  indici[[i]] <- data.frame(caso = casi$caso[i], alpha = casi$alpha[i], K = K,
                            n = n, Q0 = q0, Q1 = q1,
                            ASW_economico = s$avg.silwidth, PG = s$pearsongamma,
                            WSS = s$within.cluster.ss, CH = s$ch, Dunn = s$dunn)
}
indici <- dplyr::bind_rows(indici)
final_q0_normalizer <- indici$Q0[indici$caso == "alpha0"]
final_q1_normalizer <- indici$Q1[indici$caso == "alpha1"]
stopifnot(final_q0_normalizer > 0, final_q1_normalizer > 0)
indici$Q0norm <- indici$Q0 / final_q0_normalizer
indici$Q1norm <- indici$Q1 / final_q1_normalizer
indici$gap_Qnorm <- abs(indici$Q0norm - indici$Q1norm)
# When final K equals K0, check consistency with the native choicealpha shares.
if (K == K0) {
  for (ii in seq_len(nrow(casi))) {
    jj <- which(abs(curva_alpha_K0$alpha - casi$alpha[ii]) < 1e-12)
    if (length(jj)) stopifnot(
      abs(indici$Q0norm[ii] - curva_alpha_K0$Q0norm[jj[1L]]) < 1e-8,
      abs(indici$Q1norm[ii] - curva_alpha_K0$Q1norm[jj[1L]]) < 1e-8)
  }
}
selection_summary$alpha_gap_at_final_k <- indici$gap_Qnorm[indici$caso == "optimal_alpha"]
write.csv(selection_summary, file.path(output_dir, "Selezione_K_alpha.csv"), row.names = FALSE)
knitr::kable(indici, digits = 4, row.names = FALSE)



n_casi <- nrow(casi)
ari_matrix <- matrix(NA_real_, n_casi, n_casi,
                     dimnames = list(casi$caso, casi$caso))
ARI_coppie <- list()
for (i in seq_len(n_casi)) {
  for (j in i:n_casi) {
    ari_value <- mclust::adjustedRandIndex(
      partizioni[[casi$caso[i]]], partizioni[[casi$caso[j]]])
    ari_matrix[i, j] <- ari_matrix[j, i] <- ari_value
    if (i < j) {
      ARI_coppie[[length(ARI_coppie) + 1L]] <- data.frame(
        caso_1 = casi$caso[i], alpha_1 = casi$alpha[i],
        caso_2 = casi$caso[j], alpha_2 = casi$alpha[j],
        K = K, k0_alpha_selection = K0, n_osservazioni = n, ARI = ari_value)
    }
  }
}
ARI_coppie <- dplyr::bind_rows(ARI_coppie)
stopifnot(all(is.finite(ari_matrix)), nrow(ARI_coppie) == choose(n_casi, 2))
ari_matrix_report <- ari_matrix
dimnames(ari_matrix_report) <- list(casi$pannello, casi$pannello)
knitr::kable(ari_matrix_report, digits = 4, caption = "ARI fra i quattro alpha")
knitr::kable(ARI_coppie, digits = 4, row.names = FALSE)



ARI_long <- as.data.frame(as.table(ari_matrix), stringsAsFactors = FALSE)
names(ARI_long) <- c("caso_1", "caso_2", "ARI")
etichette_ari <- setNames(stringr::str_wrap(casi$pannello, width = 19), casi$caso)
p_ari <- ggplot(ARI_long, aes(x = caso_2, y = caso_1, fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", ARI),
                colour = ifelse(ARI > 0.65, "white", "#243746")), size = 4) +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#B97572", mid = "#F3F2EF", high = "#365F7C",
                       midpoint = 0, limits = c(-1, 1), name = "ARI") +
  scale_x_discrete(limits = casi$caso, labels = etichette_ari) +
  scale_y_discrete(limits = rev(casi$caso), labels = etichette_ari) +
  coord_equal() +
  labs(title = paste("Concordanza fra partizioni -", params$area),
       subtitle = paste("K =", K, ";", n, "osservazioni di stima"), x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", colour = "#183746"))
p_ari
ggsave(file.path(output_dir, "ARI_mappa_calore.png"), p_ari,
       width = 8, height = 7, dpi = 200, bg = "white")
ggsave(file.path(output_dir, "ARI_mappa_calore.pdf"), p_ari,
       width = 8, height = 7)



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

Mappe <- Data_clusters %>%
  dplyr::select(id_osservazione, LAU_Code, all_of(casi$caso)) %>%
  tidyr::pivot_longer(cols = all_of(casi$caso),
                      names_to = "caso", values_to = "Cluster") %>%
  mutate(Pannello = factor(unname(nomi_pannelli[caso]), levels = casi$pannello),
         Cluster = factor(Cluster, levels = seq_len(n_colori)))

sottotitolo <- if (include_milan) {
  paste("K =", K, "complessivi; Milano incluso nella stima")
} else paste("K =", K, "nella stima; Milano separato come C", K + 1L)
p_mappe <- ggplot(Mappe) +
  geom_sf(aes(fill = Cluster), linewidth = 0.045, colour = "#F5F3EF") +
  facet_wrap(~ Pannello, ncol = 2) +
  scale_fill_manual(values = cluster_palette, labels = legenda_cluster,
                    drop = FALSE, name = "Cluster") +
  labs(title = paste("Ward spaziale -", params$area),
       subtitle = paste(sottotitolo, "-", alpha_criterion_description), x = "", y = "",
       caption = "Colori abbinati alla partizione dell'alpha scelto; etichette native conservate.") +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom", legend.key.height = grid::unit(0.4, "cm"),
        strip.text = element_text(face = "bold", colour = "#183746", size = 12),
        plot.title = element_text(face = "bold", colour = "#183746", size = 15),
        plot.subtitle = element_text(colour = "#596772"),
        plot.caption = element_text(colour = "#687580", hjust = 0),
        plot.margin = margin(12, 12, 12, 12)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
p_mappe

ggsave(file.path(output_dir, "Mappe_quattro_alpha.png"), p_mappe,
       width = 13, height = 11, dpi = 200, bg = "white")
ggsave(file.path(output_dir, "Mappe_quattro_alpha.pdf"), p_mappe,
       width = 13, height = 11)



dimensioni_cluster <- dplyr::bind_rows(lapply(casi$caso, function(nm) {
  data.frame(caso = nm, cluster = seq_len(K),
             osservazioni = as.integer(table(factor(partizioni[[nm]], levels = seq_len(K)))))
}))
knitr::kable(dimensioni_cluster, row.names = FALSE)
if (include_milan) {
  knitr::kable(partizioni[partizioni$LAU_Code == "015146", ], row.names = FALSE)
}




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
      n_municipalities = as.integer(table(factor(gruppi, levels = seq_len(kk)))),
      ramo_colore = mappa_rami)
  }
  ARI_multik_lista[[length(ARI_multik_lista) + 1L]] <- data.frame(
    K = kk, alpha_0 = 0, optimal_alpha = optimal_alpha, n_municipalities = n,
    ARI = mclust::adjustedRandIndex(g0, gs))
}
Abbinamento_colori_multik <- bind_rows(Abbinamento_colori_multik)
Abbinamento_colori_multik$colore_hex <- unname(multi_k_palette[
  as.character(Abbinamento_colori_multik$ramo_colore)])
Dimensioni_multik <- bind_rows(Dimensioni_multik_lista)
ARI_multik <- bind_rows(ARI_multik_lista)
stopifnot(nrow(Casi_multik) == 28L, all(Dimensioni_multik$n_municipalities > 0))
write.csv(Casi_multik, file.path(output_dir, "Casi_multik_K2_K15.csv"), row.names = FALSE)
write.csv(Partizioni_multik_confronto,
          file.path(output_dir, "Partizioni_alpha0_alphaottimo_K2_K15.csv"), row.names = FALSE)
write.csv(Abbinamento_colori_multik,
          file.path(output_dir, "Abbinamento_colori_K2_K15.csv"), row.names = FALSE)
write.csv(Dimensioni_multik, file.path(output_dir, "Dimensioni_gruppi_K2_K15.csv"), row.names = FALSE)
write.csv(ARI_multik, file.path(output_dir, "ARI_alpha0_alphaottimo_K2_K15.csv"), row.names = FALSE)



# Keep the spatial table for the 28 panels local to avoid bloating the workspace.
local({
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
})
knitr::include_graphics(file.path(output_dir, "Mappe_alpha0_alphaottimo_K2_K15.png"))




calcola_ARI_stabilita_K <- function(caso_alpha) {
  campi <- paste0(caso_alpha, "_K", K_confronto)
  stopifnot(all(campi %in% names(Partizioni_multik_confronto)),
            nrow(Partizioni_multik_confronto) == n,
            identical(Partizioni_multik_confronto$id_osservazione, seq_len(n)))
  m <- length(K_confronto)
  etichette_K <- paste0("K = ", K_confronto)
  matrice <- matrix(NA_real_, m, m, dimnames = list(etichette_K, etichette_K))
  for (ii in seq_len(m)) {
    for (jj in seq_len(ii)) {
      valore <- mclust::adjustedRandIndex(
        Partizioni_multik_confronto[[campi[ii]]],
        Partizioni_multik_confronto[[campi[jj]]])
      matrice[ii, jj] <- matrice[jj, ii] <- valore
    }
  }
  stopifnot(all(is.finite(matrice)),
            max(abs(matrice - t(matrice))) < 1e-12,
            all(abs(diag(matrice) - 1) < 1e-12),
            all(matrice >= -1 - 1e-12 & matrice <= 1 + 1e-12))
  matrice
}
ARI_stabilita_K <- setNames(lapply(c("alpha0", "optimal_alpha"),
                                   calcola_ARI_stabilita_K), c("alpha0", "optimal_alpha"))
ARI_stabilita_K_alpha0 <- ARI_stabilita_K$alpha0
ARI_stabilita_K_alphaottimo <- ARI_stabilita_K$optimal_alpha

ARI_stabilita_K_long <- bind_rows(lapply(names(ARI_stabilita_K), function(cc) {
  tab <- expand.grid(K_riga = K_confronto, K_colonna = K_confronto,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # expand.grid varia prima K_riga, come l'ordine per colonne di una matrice R.
  tab$ARI <- as.vector(ARI_stabilita_K[[cc]])
  tab$caso <- cc
  tab$alpha <- if (cc == "alpha0") 0 else optimal_alpha
  tab$n_municipalities <- n
  tab
}))
ARI_stabilita_K_coppie <- ARI_stabilita_K_long %>%
  filter(K_riga < K_colonna) %>%
  mutate(tagli_consecutivi = K_colonna == K_riga + 1L) %>%
  arrange(caso, K_riga, K_colonna)
stopifnot(nrow(ARI_stabilita_K_coppie) == 2L * choose(length(K_confronto), 2L))
write.csv(data.frame(K = K_confronto, ARI_stabilita_K_alpha0, check.names = FALSE),
          file.path(output_dir, "ARI_stabilita_K_alpha0.csv"), row.names = FALSE)
write.csv(data.frame(K = K_confronto, ARI_stabilita_K_alphaottimo, check.names = FALSE),
          file.path(output_dir, "ARI_stabilita_K_alphaottimo.csv"), row.names = FALSE)
write.csv(ARI_stabilita_K_coppie,
          file.path(output_dir, "ARI_stabilita_K_coppie.csv"), row.names = FALSE)




tabella_ARI_stabilita <- function(matrice, didascalia) {
  dati <- data.frame(K = K_confronto, matrice, check.names = FALSE)
  widget <- DT::datatable(dati, rownames = FALSE, class = "stripe hover compact",
                          extensions = "Buttons", caption = didascalia,
                          options = list(dom = "Bt", buttons = list(
                            list(extend = "copy", exportOptions = list(orthogonal = "export")),
                            list(extend = "csv", exportOptions = list(orthogonal = "export")), "print"),
                            paging = FALSE, ordering = FALSE, searching = FALSE, scrollX = TRUE))
  DT::formatRound(widget, columns = colnames(matrice), digits = 3)
}
tabella_ARI_stabilita_alpha0 <- tabella_ARI_stabilita(ARI_stabilita_K_alpha0,
                                                      "ARI tra K=2...15, alpha = 0")
htmlwidgets::saveWidget(tabella_ARI_stabilita_alpha0,
                        file.path(output_dir, "ARI_stabilita_K_alpha0.html"), selfcontained = TRUE,
                        title = "Stabilita rispetto a K: alpha = 0")
tabella_ARI_stabilita_alpha0



tabella_ARI_stabilita_alphaottimo <- tabella_ARI_stabilita(ARI_stabilita_K_alphaottimo,
                                                           paste0("ARI tra K=2...15, alpha ottimo = ", format(optimal_alpha, digits = 6),
                                                                  " (costante, selezionato a K0 = ", K0, ")"))
htmlwidgets::saveWidget(tabella_ARI_stabilita_alphaottimo,
                        file.path(output_dir, "ARI_stabilita_K_alphaottimo.html"), selfcontained = TRUE,
                        title = "Stabilita rispetto a K: alpha ottimo")
tabella_ARI_stabilita_alphaottimo



etichette_ARI_stabilita <- c("alpha = 0",
                             paste0("alpha ottimo = ", format(optimal_alpha, digits = 6)))
dati_ARI_stabilita_grafico <- ARI_stabilita_K_long %>%
  mutate(Pannello = factor(caso, levels = c("alpha0", "optimal_alpha"),
                           labels = etichette_ARI_stabilita))
p_ARI_stabilita <- ggplot(dati_ARI_stabilita_grafico,
                          aes(x = factor(K_colonna, levels = K_confronto),
                              y = factor(K_riga, levels = rev(K_confronto)), fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.3f", ARI),
                colour = ifelse(ARI > 0.65, "white", "#243746")), size = 2.5) +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#B97572", mid = "#F3F2EF", high = "#365F7C",
                       midpoint = 0, limits = c(-1, 1), name = "ARI") +
  facet_wrap(~ Pannello, ncol = 2) + coord_equal() +
  labs(title = paste("Stabilita delle partizioni al variare di K -", params$area),
       subtitle = paste(n, "comuni di stima; alpha costante in ciascuna matrice"),
       x = "K della partizione in colonna", y = "K della partizione in riga") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", colour = "#183746"),
        plot.title = element_text(face = "bold", colour = "#183746"))
p_ARI_stabilita
ggsave(file.path(output_dir, "ARI_stabilita_K_mappe_calore.png"),
       p_ARI_stabilita, width = 14, height = 7.5, dpi = 200, bg = "white")
ggsave(file.path(output_dir, "ARI_stabilita_K_mappe_calore.pdf"),
       p_ARI_stabilita, width = 14, height = 7.5)



# Dedicated descriptive grid: exactly 11 alpha values, including both endpoints.
griglia_alpha_ARI <- (0:10) / 10
etichette_alpha_ARI <- sprintf("%.2f", griglia_alpha_ARI)
id_alpha_ARI <- paste0("alpha_", sprintf("%03d", as.integer(round(100 * griglia_alpha_ARI))))
Griglia_alpha_ARI <- data.frame(alpha = griglia_alpha_ARI, id_alpha = id_alpha_ARI,
                                albero_riutilizzato = FALSE, origine_albero = NA_character_, stringsAsFactors = FALSE)
Partizioni_alpha_ARI <- setNames(vector("list", length(griglia_alpha_ARI)), id_alpha_ARI)
for (ia_ARI in seq_along(griglia_alpha_ARI)) {
  aa_ARI <- griglia_alpha_ARI[ia_ARI]
  # Un solo albero per alpha, poi tutti i tagli K=2...15.
  # Si riusano gli alberi del confronto principale solo ad alpha identico.
  indici_riuso_ARI_alpha <- which(casi$alpha == aa_ARI)
  if (length(indici_riuso_ARI_alpha)) {
    indice_riuso_ARI_alpha <- indici_riuso_ARI_alpha[1L]
    albero_ARI_alpha <- alberi[[indice_riuso_ARI_alpha]]
    Griglia_alpha_ARI$albero_riutilizzato[ia_ARI] <- TRUE
    Griglia_alpha_ARI$origine_albero[ia_ARI] <- casi$caso[indice_riuso_ARI_alpha]
  } else {
    albero_ARI_alpha <- ClustGeo::hclustgeo(D0_std, D1_std, alpha = aa_ARI, wt = wt)
    Griglia_alpha_ARI$origine_albero[ia_ARI] <- "Fit ClustGeo sulla griglia diagnostica"
  }
  tagli_ARI_alpha <- cutree(albero_ARI_alpha, k = K_confronto)
  stopifnot(nrow(tagli_ARI_alpha) == n, ncol(tagli_ARI_alpha) == length(K_confronto),
            identical(colnames(tagli_ARI_alpha), as.character(K_confronto)),
            all(vapply(seq_along(K_confronto), function(ii)
              length(unique(tagli_ARI_alpha[, ii])), integer(1)) == K_confronto))
  Partizioni_alpha_ARI[[ia_ARI]] <- tagli_ARI_alpha
}
write.csv(Griglia_alpha_ARI, file.path(output_dir, "Griglia_alpha_ARI.csv"), row.names = FALSE)
saveRDS(list(comuni = partizioni[c("id_osservazione", "LAU_Code")],
             K = K_confronto, griglia_alpha = Griglia_alpha_ARI, partizioni = Partizioni_alpha_ARI),
        file.path(output_dir, "Partizioni_stabilita_alpha_K2_K15.rds"))



ARI_stabilita_alpha <- setNames(vector("list", length(K_confronto)), as.character(K_confronto))
for (kk_ARI_alpha in K_confronto) {
  n_alpha_ARI <- length(griglia_alpha_ARI)
  alpha_ari_matrix_names <- paste0("alpha = ", etichette_alpha_ARI)
  alpha_ari_matrix <- matrix(NA_real_, n_alpha_ARI, n_alpha_ARI,
                             dimnames = list(alpha_ari_matrix_names, alpha_ari_matrix_names))
  for (ii_ARI_alpha in seq_len(n_alpha_ARI)) {
    for (jj_ARI_alpha in seq_len(ii_ARI_alpha)) {
      alpha_ari_value <- mclust::adjustedRandIndex(
        Partizioni_alpha_ARI[[ii_ARI_alpha]][, as.character(kk_ARI_alpha)],
        Partizioni_alpha_ARI[[jj_ARI_alpha]][, as.character(kk_ARI_alpha)])
      alpha_ari_matrix[ii_ARI_alpha, jj_ARI_alpha] <-
        alpha_ari_matrix[jj_ARI_alpha, ii_ARI_alpha] <- alpha_ari_value
    }
  }
  stopifnot(all(is.finite(alpha_ari_matrix)),
            max(abs(alpha_ari_matrix - t(alpha_ari_matrix))) < 1e-12,
            all(abs(diag(alpha_ari_matrix) - 1) < 1e-12),
            all(alpha_ari_matrix >= -1 - 1e-12 & alpha_ari_matrix <= 1 + 1e-12))
  ARI_stabilita_alpha[[as.character(kk_ARI_alpha)]] <- alpha_ari_matrix
  write.csv(data.frame(alpha_riga = griglia_alpha_ARI, alpha_ari_matrix, check.names = FALSE),
            file.path(output_dir, paste0("ARI_stabilita_alpha_K", kk_ARI_alpha, ".csv")), row.names = FALSE)
}
# Match the main comparison: same alpha values (0, 0.10, 1) and same final K.
indici_riscontro_ARI_alpha <- c(1L, 2L, length(griglia_alpha_ARI))
casi_riscontro_ARI_alpha <- c("alpha0", "alpha010", "alpha1")
stopifnot(isTRUE(all.equal(
  unname(ARI_stabilita_alpha[[as.character(K)]][indici_riscontro_ARI_alpha, indici_riscontro_ARI_alpha]),
  unname(ari_matrix[casi_riscontro_ARI_alpha, casi_riscontro_ARI_alpha]), tolerance = 1e-12)))

ARI_stabilita_alpha_long <- bind_rows(lapply(K_confronto, function(kk) {
  tab <- expand.grid(indice_riga = seq_along(griglia_alpha_ARI),
                     indice_colonna = seq_along(griglia_alpha_ARI), KEEP.OUT.ATTRS = FALSE)
  tab$K <- kk
  tab$alpha_riga <- griglia_alpha_ARI[tab$indice_riga]
  tab$alpha_colonna <- griglia_alpha_ARI[tab$indice_colonna]
  tab$ARI <- as.vector(ARI_stabilita_alpha[[as.character(kk)]])
  tab$n_municipalities <- n
  tab
}))
ARI_stabilita_alpha_coppie <- ARI_stabilita_alpha_long %>%
  filter(indice_riga < indice_colonna) %>%
  mutate(alpha_consecutivi = indice_colonna == indice_riga + 1L) %>%
  arrange(K, alpha_riga, alpha_colonna)
stopifnot(nrow(ARI_stabilita_alpha_coppie) == length(K_confronto) * choose(length(griglia_alpha_ARI), 2L))
write.csv(ARI_stabilita_alpha_coppie,
          file.path(output_dir, "ARI_stabilita_alpha_coppie.csv"), row.names = FALSE)
saveRDS(ARI_stabilita_alpha, file.path(output_dir, "ARI_stabilita_alpha_matrici.rds"))



ARI_stabilita_alpha_tabella <- bind_rows(lapply(K_confronto, function(kk) {
  data.frame(K = kk, alpha_riga = griglia_alpha_ARI,
             ARI_stabilita_alpha[[as.character(kk)]], check.names = FALSE)
}))
dati_tabella_ARI_alpha <- ARI_stabilita_alpha_tabella %>%
  mutate(K = factor(K, levels = K_confronto),
         alpha_riga = factor(sprintf("%.2f", alpha_riga), levels = etichette_alpha_ARI))
tabella_ARI_stabilita_alpha <- DT::datatable(dati_tabella_ARI_alpha, rownames = FALSE,
                                             filter = "top", extensions = "Buttons", class = "stripe hover compact",
                                             caption = "Matrici ARI a K fisso: selezionare K; alpha=0.00, 0.10, ..., 1.00",
                                             options = list(dom = "Blfrtip", buttons = list(
                                               list(extend = "copy", exportOptions = list(orthogonal = "export")),
                                               list(extend = "csv", exportOptions = list(orthogonal = "export")), "print"),
                                               pageLength = 11, lengthMenu = c(11, 22, 55, 154), ordering = FALSE,
                                               scrollX = TRUE, deferRender = TRUE))
tabella_ARI_stabilita_alpha <- DT::formatRound(tabella_ARI_stabilita_alpha,
                                               columns = colnames(ARI_stabilita_alpha[[1L]]), digits = 3)
htmlwidgets::saveWidget(tabella_ARI_stabilita_alpha,
                        file.path(output_dir, "ARI_stabilita_alpha_matrici.html"), selfcontained = TRUE,
                        title = "Stabilita rispetto ad alpha, a K fisso")
tabella_ARI_stabilita_alpha

