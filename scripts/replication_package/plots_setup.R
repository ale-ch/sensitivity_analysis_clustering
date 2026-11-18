casi <- data.frame(
  caso = c("alpha0", "optimal_alpha", "alpha010", "alpha1"),
  alpha = c(0, optimal_alpha, alpha_fisso, 1),
  stringsAsFactors = FALSE)
casi$pannello <- c("alpha = 0",
                   paste0("alpha ", if (params$alpha_selection_mode == "automatic") "optimal" else "scelto",
                          " a K0=", K0, ": ", format(optimal_alpha, digits = 5)),
                   "alpha = 0.10", "alpha = 1")
casi$K <- K
casi$K0_alpha_selection_mode <- K0
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
# Campi principali: K finale; campi storici _P8: realmente otto gruppi.
DPCM_red2$D0_P_Kfinale <- partizioni$alpha0
DPCM_red2$D0D1_P_Kfinale <- partizioni$optimal_alpha
DPCM_red2$D0D1_P_Kfinale_alpha010 <- partizioni$alpha010
DPCM_red2$D0_P8 <- cutree(tree_D0, 8L)
DPCM_red2$D0D1_P8 <- cutree(treeD0D1, 8L)
DPCM_red2$D0D1_P8_alpha010 <- cutree(treeD0D1_alpha010, 8L)
stopifnot(length(unique(DPCM_red2$D0_P8)) == 8L,
          length(unique(DPCM_red2$D0D1_P8)) == 8L)




