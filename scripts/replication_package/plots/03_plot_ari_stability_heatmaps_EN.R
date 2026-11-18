# 4x4 ARI Heatmap (4 Alpha Configurations)

casi <- data.frame(
  caso = c("alpha0", "optimal_alpha", "alpha010", "alpha1"),
  alpha = c(0, optimal_alpha, params$manual_alpha , 1),
  stringsAsFactors = FALSE)
casi$pannello <- c("alpha = 0",
                   paste0("alpha ", if (params$alpha_selection_mode == "automatic") "optimal" else "selected",
                          " at K0=", K0, ": ", format(optimal_alpha, digits = 5)),
                   "alpha = 0.10", "alpha = 1")
casi$K <- K
casi$K0_scelta_alpha <- K0
casi$criterio_alpha <- params$alpha_criterion

griglia_alpha_ARI <- (0:10) / 10
etichette_alpha_ARI <- sprintf("%.2f", griglia_alpha_ARI)
id_alpha_ARI <- paste0("alpha_", sprintf("%03d", as.integer(round(100 * griglia_alpha_ARI))))
Griglia_alpha_ARI <- data.frame(alpha = griglia_alpha_ARI, id_alpha = id_alpha_ARI,
                                albero_riutilizzato = FALSE, origine_albero = NA_character_, stringsAsFactors = FALSE)
Partizioni_alpha_ARI <- setNames(vector("list", length(griglia_alpha_ARI)), id_alpha_ARI)
for (ia_ARI in seq_along(griglia_alpha_ARI)) {
  aa_ARI <- griglia_alpha_ARI[ia_ARI]
  # Single tree per alpha, then all K=2...15 cuts.
  # Main comparison trees are reused only at identical alpha.
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

n_casi <- nrow(casi)
ARI_matrice <- matrix(NA_real_, n_casi, n_casi,
                      dimnames = list(casi$caso, casi$caso))
ARI_coppie <- list()
for (i in seq_len(n_casi)) {
  for (j in i:n_casi) {
    valore_ari <- mclust::adjustedRandIndex(
      partizioni[[casi$caso[i]]], partizioni[[casi$caso[j]]])
    ARI_matrice[i, j] <- ARI_matrice[j, i] <- valore_ari
    if (i < j) {
      ARI_coppie[[length(ARI_coppie) + 1L]] <- data.frame(
        caso_1 = casi$caso[i], alpha_1 = casi$alpha[i],
        caso_2 = casi$caso[j], alpha_2 = casi$alpha[j],
        K = K, K0_scelta_alpha = K0, n_osservazioni = n, ARI = valore_ari)
    }
  }
}
ARI_coppie <- dplyr::bind_rows(ARI_coppie)
stopifnot(all(is.finite(ARI_matrice)), nrow(ARI_coppie) == choose(n_casi, 2))
ARI_matrice_report <- ARI_matrice
dimnames(ARI_matrice_report) <- list(casi$pannello, casi$pannello)
ARI_long <- as.data.frame(as.table(ARI_matrice), stringsAsFactors = FALSE)
names(ARI_long) <- c("caso_1", "caso_2", "ARI")
etichette_ari <- setNames(stringr::str_wrap(casi$pannello, width = 19), casi$caso)

p_ari <- ggplot(ARI_long, aes(x = caso_2, y = caso_1, fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", ARI), colour = ifelse(ARI > 0.65, "white", "#243746")), size = 4) +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#B97572", mid = "#F3F2EF", high = "#365F7C", midpoint = 0, limits = c(-1, 1), name = "ARI") +
  scale_x_discrete(limits = casi$caso, labels = etichette_ari) +
  scale_y_discrete(limits = rev(casi$caso), labels = etichette_ari) +
  coord_equal() +
  labs(title = paste("Partition agreement -", params$area),
       subtitle = paste("K =", K, ";", n, "estimation observations"), x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", colour = "#183746"))

ggsave(file.path(output_dir, "ARI_mappa_calore.png"), p_ari, width = 8, height = 7, dpi = 200, bg = "white")

# K-Stability ARI Heatmaps (14x14)
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
  # expand.grid varies K_riga first, matching R matrix column-major order.
  tab$ARI <- as.vector(ARI_stabilita_K[[cc]])
  tab$caso <- cc
  tab$alpha <- if (cc == "alpha0") 0 else optimal_alpha
  tab$n_comuni <- n
  tab
}))
ARI_stabilita_K_coppie <- ARI_stabilita_K_long %>%
  filter(K_riga < K_colonna) %>%
  mutate(tagli_consecutivi = K_colonna == K_riga + 1L) %>%
  arrange(caso, K_riga, K_colonna)
stopifnot(nrow(ARI_stabilita_K_coppie) == 2L * choose(length(K_confronto), 2L))

etichette_ARI_stabilita <- c("alpha = 0",
                             paste0("optimal alpha = ", format(optimal_alpha, digits = 6)))

dati_ARI_stabilita_grafico <- ARI_stabilita_K_long %>%
  mutate(Pannello = factor(caso, levels = c("alpha0", "optimal_alpha"),
                           labels = etichette_ARI_stabilita))

p_ARI_stabilita <- ggplot(dati_ARI_stabilita_grafico,
                          aes(x = factor(K_colonna, levels = K_confronto),
                              y = factor(K_riga, levels = rev(K_confronto)), fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.3f", ARI), colour = ifelse(ARI > 0.65, "white", "#243746")), size = 2.5) +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#B97572", mid = "#F3F2EF", high = "#365F7C", midpoint = 0, limits = c(-1, 1), name = "ARI") +
  facet_wrap(~ Pannello, ncol = 2) + coord_equal() +
  labs(title = paste("Partition stability across varying K -", params$area),
       subtitle = paste(n, "estimation municipalities; constant alpha in each matrix"),
       x = "K of column partition", y = "K of row partition") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", colour = "#183746"),
        plot.title = element_text(face = "bold", colour = "#183746"))

ggsave(file.path(output_dir, "ARI_stabilita_K_mappe_calore.png"), p_ARI_stabilita, width = 14, height = 7.5, dpi = 200, bg = "white")










# Alpha-Stability ARI Heatmaps (Fixed K)
ARI_stabilita_alpha <- setNames(vector("list", length(K_confronto)), as.character(K_confronto))
for (kk_ARI_alpha in K_confronto) {
  n_alpha_ARI <- length(griglia_alpha_ARI)
  nomi_matrice_ARI_alpha <- paste0("alpha = ", etichette_alpha_ARI)
  matrice_ARI_alpha <- matrix(NA_real_, n_alpha_ARI, n_alpha_ARI,
                              dimnames = list(nomi_matrice_ARI_alpha, nomi_matrice_ARI_alpha))
  for (ii_ARI_alpha in seq_len(n_alpha_ARI)) {
    for (jj_ARI_alpha in seq_len(ii_ARI_alpha)) {
      valore_ARI_alpha <- mclust::adjustedRandIndex(
        Partizioni_alpha_ARI[[ii_ARI_alpha]][, as.character(kk_ARI_alpha)],
        Partizioni_alpha_ARI[[jj_ARI_alpha]][, as.character(kk_ARI_alpha)])
      matrice_ARI_alpha[ii_ARI_alpha, jj_ARI_alpha] <-
        matrice_ARI_alpha[jj_ARI_alpha, ii_ARI_alpha] <- valore_ARI_alpha
    }
  }
  stopifnot(all(is.finite(matrice_ARI_alpha)),
            max(abs(matrice_ARI_alpha - t(matrice_ARI_alpha))) < 1e-12,
            all(abs(diag(matrice_ARI_alpha) - 1) < 1e-12),
            all(matrice_ARI_alpha >= -1 - 1e-12 & matrice_ARI_alpha <= 1 + 1e-12))
  ARI_stabilita_alpha[[as.character(kk_ARI_alpha)]] <- matrice_ARI_alpha
  write.csv(data.frame(alpha_riga = griglia_alpha_ARI, matrice_ARI_alpha, check.names = FALSE),
            file.path(output_dir, paste0("ARI_stabilita_alpha_K", kk_ARI_alpha, ".csv")), row.names = FALSE)
}
griglia_alpha_ARI <- (0:10) / 10
etichette_alpha_ARI <- sprintf("%.2f", griglia_alpha_ARI)
ARI_stabilita_alpha_long <- bind_rows(lapply(K_confronto, function(kk) {
  tab <- expand.grid(indice_riga = seq_along(griglia_alpha_ARI),
                     indice_colonna = seq_along(griglia_alpha_ARI), KEEP.OUT.ATTRS = FALSE)
  tab$K <- kk
  tab$alpha_riga <- griglia_alpha_ARI[tab$indice_riga]
  tab$alpha_colonna <- griglia_alpha_ARI[tab$indice_colonna]
  tab$ARI <- as.vector(ARI_stabilita_alpha[[as.character(kk)]])
  tab$n_comuni <- n
  tab
}))

dati_ARI_alpha_grafico <- ARI_stabilita_alpha_long %>%
  mutate(K_pannello = factor(K, levels = K_confronto,
                             labels = paste0("K = ", K_confronto)))

p_ARI_stabilita_alpha <- ggplot(dati_ARI_alpha_grafico,
                                aes(x = factor(indice_colonna, levels = seq_along(griglia_alpha_ARI), labels = etichette_alpha_ARI),
                                    y = factor(indice_riga, levels = rev(seq_along(griglia_alpha_ARI)), labels = rev(etichette_alpha_ARI)), fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.3f", ARI), colour = ifelse(ARI > 0.65, "white", "#243746")), size = 2.5) +
  scale_colour_identity() +
  scale_fill_gradient2(low = "#B97572", mid = "#F3F2EF", high = "#365F7C", midpoint = 0, limits = c(-1, 1), name = "ARI") +
  facet_wrap(~ K_pannello, ncol = 3) + coord_equal() +
  labs(title = paste("Stability with respect to alpha at fixed K -", params$area),
       subtitle = paste(n, "estimation municipalities; alpha=0.00, 0.10, ..., 1.00"),
       x = "Alpha of column partition", y = "Alpha of row partition") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", colour = "#183746"),
        plot.title = element_text(face = "bold", colour = "#183746"))

ggsave(file.path(output_dir, "ARI_stabilita_alpha_mappe_calore.png"), p_ARI_stabilita_alpha, width = 18, height = 30, dpi = 160, bg = "white")