# In the event of a tie within tolerance, prefer the smaller value.
scegli_minimo <- function(valore, candidato, tolleranza = 1e-12) {
  ok <- is.finite(valore)
  if (!any(ok)) stop("Nessun punteggio finito disponibile per la selezione.")
  migliore <- min(valore[ok])
  min(candidato[ok & valore <= migliore + tolleranza])
}

# Native ClustGeo shares and the paper’s sole selection criterion.
tabella_choicealpha <- function(curva, k) {
  tab <- data.frame(K = k, alpha = curva$range.alpha,
                    Q0 = curva$Q[, 1L], Q1 = curva$Q[, 2L],
                    Q0norm = curva$Qnorm[, 1L], Q1norm = curva$Qnorm[, 2L],
                    gap_Qnorm = abs(curva$Qnorm[, 1L] - curva$Qnorm[, 2L]))
  tab
}




T0_std <- ClustGeo::inertdiss(D0_std, wt = wt)
T1_std <- ClustGeo::inertdiss(D1_std, wt = wt)
stopifnot(T0_std > 0, T1_std > 0,
          T0_std <= 0.5 + 1e-10, T1_std <= 0.5 + 1e-10)

alpha_criterion_description <- "Equilibrio assoluto (paper)"

curva_gomito <- function(albero, alpha, fase) {
  tab <- dplyr::bind_rows(lapply(seq_len(max(K_candidati)), function(kk) {
    g <- cutree(albero, k = kk)
    data.frame(K = kk,
               W0 = ClustGeo::withindiss(D0_std, g, wt = wt),
               W1 = ClustGeo::withindiss(D1_std, g, wt = wt))
  }))
  tab$W_mista <- (1 - alpha) * tab$W0 + alpha * tab$W1
  T_mista <- (1 - alpha) * T0_std + alpha * T1_std
  tab$Q_mista <- 1 - tab$W_mista / T_mista
  escursione <- tab$W_mista[1L] - tail(tab$W_mista, 1L)
  if (!is.finite(escursione) || escursione <= 0) {
    stop("Curva di inerzia priva di escursione: il gomito non e identificabile.")
  }
  tab$x_norm <- (tab$K - 1) / (max(tab$K) - 1)
  tab$W_norm <- (tab$W_mista - tail(tab$W_mista, 1L)) / escursione
  tab$score_gomito <- (1 - tab$x_norm) - tab$W_norm
  tab$candidato <- tab$K %in% K_candidati
  tab$fase <- fase
  tab$alpha <- alpha
  # In un Ward le inerzie interne devono diminuire quando K aumenta.
  stopifnot(all(is.finite(tab$W_mista)),
            all(diff(tab$W_mista) <= 1e-9 * max(1, T_mista)))
  tab
}

proponi_K <- function(tab) {
  candidati <- tab[tab$candidato, ]
  if (max(candidati$score_gomito) <= 1e-10) {
    warning("Gomito non marcato: ispezionare curva e dendrogramma; parita risolte a K minore.")
  }
  as.integer(scegli_minimo(-candidati$score_gomito, candidati$K))
}


## Stage 1: K0 using economic information only

# Also retain the tree on the original D0 scale for historical outputs.
tree_D0 <- hclustgeo(D0 = D0, wt = rep(1/n, n))
diagnostica_K0 <- curva_gomito(tree_D0, alpha = 0, fase = "1. K0 su D0")
K0_automatico <- proponi_K(diagnostica_K0)
K0 <- if (params$k_selection_mode == "automatic") K0_automatico else as.integer(params$manual_k0)
diagnostica_K0$selezionato <- diagnostica_K0$K == K0
knitr::kable(data.frame(K0_usato = K0, K0_gomito_automatico = K0_automatico,
                        Modalita = params$k_selection_mode))





cr_K0 <- choicealpha(D0_std, D1_std, range.alpha = range_alpha,
                     K = K0, wt = wt, graph = FALSE)
curva_alpha_K0 <- tabella_choicealpha(cr_K0, K0)
stopifnot(all(is.finite(curva_alpha_K0$gap_Qnorm)))
automatic_alpha <- scegli_minimo(curva_alpha_K0$gap_Qnorm, curva_alpha_K0$alpha)
optimal_alpha <- if (params$alpha_selection_mode == "automatic") {
  automatic_alpha
} else params$manual_alpha
righe_alpha <- sort(unique(c(1L, nrow(curva_alpha_K0),
                             which(abs(curva_alpha_K0$alpha - automatic_alpha) < 1e-12),
                             which(abs(curva_alpha_K0$alpha - optimal_alpha) < 1e-12),
                             which(abs(curva_alpha_K0$alpha - alpha_fisso) < 1e-12))))
knitr::kable(curva_alpha_K0[righe_alpha, ], digits = 5, row.names = FALSE)



tree_alpha_scelto <- hclustgeo(D0_std, D1_std, alpha = optimal_alpha, wt = wt)
final_k_diagnostic <- curva_gomito(tree_alpha_scelto, alpha = optimal_alpha,
                                   fase = "3. K finale ad alpha scelto")
automatic_final_k <- proponi_K(final_k_diagnostic)
K <- kopt <- if (params$k_selection_mode == "automatic") {
  automatic_final_k
} else as.integer(params$manual_final_k)
final_k_diagnostic$selezionato <- final_k_diagnostic$K == K
stopifnot(K %in% K_candidati, K0 %in% K_candidati)
k_alpha_selection <- data.frame(
  passaggio = c("1. K0 su D0", "2. alpha a K0 fissato", "3. K finale ad alpha fissato"),
  K = c(K0, K0, K), alpha = c(0, optimal_alpha, optimal_alpha),
  modalita = c(params$k_selection_mode, params$alpha_selection_mode, params$k_selection_mode),
  proposta_automatica = c(K0_automatico, automatic_alpha, automatic_final_k),
  criterio = c("Gomito: massima distanza dalla corda; verificare dendrogramma",
               alpha_criterion_description, "Gomito dell'inerzia mista; verificare dendrogramma"))
selection_summary <- data.frame(area = params$area, Milano_nella_stima = include_milan,
                                n_osservazioni = n, unique_municipalities = dplyr::n_distinct(DPCM_red2$LAU_Code),
                                n_variables = length(clustering_variables), unita_statistica = "comune unico",
                                K0 = K0, alpha_scelto = optimal_alpha, final_k = K, alpha_criterion = params$alpha_criterion,
                                K0_automatico = K0_automatico, automatic_alpha = automatic_alpha,
                                automatic_final_k = automatic_final_k,
                                modalita_K = params$k_selection_mode, modalita_alpha = params$alpha_selection_mode,
                                alpha_step = params$alpha_step, gap_minimo_a_K0 = min(curva_alpha_K0$gap_Qnorm),
                                T0_std = T0_std, T1_std = T1_std, alpha_recalculated_after_final_k = FALSE)
write.csv(diagnostica_K0, file.path(output_dir, "Diagnostica_K0.csv"), row.names = FALSE)
write.csv(curva_alpha_K0, file.path(output_dir, "Curve_alpha_selezione_K0.csv"), row.names = FALSE)
write.csv(final_k_diagnostic, file.path(output_dir, "Diagnostica_K_finale.csv"), row.names = FALSE)
write.csv(k_alpha_selection, file.path(output_dir, "Passaggi_selezione.csv"), row.names = FALSE)
write.csv(selection_summary, file.path(output_dir, "Selezione_K_alpha.csv"), row.names = FALSE)
writeLines(c(paste("Area:", params$area), paste("Milano nella stima:", include_milan),
             paste("K0:", K0), paste("Alpha scelto a K0:", format(optimal_alpha, digits = 12)),
             paste("K finale:", K), paste("Modalita K:", params$k_selection_mode),
             paste("Criterio alpha:", params$alpha_criterion), paste("Modalita alpha:", params$alpha_selection_mode),
             "Sequenza: K0 su D0 -> alpha a K0 fissato -> K finale ad alpha fissato.",
             "Il gomito automatico e una regola operativa esplicita; verificare i dendrogrammi.",
             "Nessuna riottimizzazione di alpha dopo il passaggio 3."),
           file.path(output_dir, "Scelta_finale.txt"))



cat(sprintf("**Scelta per questo run: K finale = %d; alpha = %.6g.** K0 = %d.\n\n",
            K, optimal_alpha, K0))
cat(paste("Criterio:", alpha_criterion_description, "\n\n"))
print(knitr::kable(k_alpha_selection, digits = 6, row.names = FALSE))
if (K != K0) cat("\n\nK finale differisce da K0: alpha rimane quello del passaggio 2.\n")



