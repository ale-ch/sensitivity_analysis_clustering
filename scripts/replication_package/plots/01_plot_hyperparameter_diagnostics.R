# Mantiene anche l'albero nella scala D0 originale per gli output storici.
tree_D0 <- hclustgeo(D0 = D0, wt = rep(1/n, n))
diagnostica_K0 <- curva_gomito(tree_D0, alpha = 0, fase = "1. K0 su D0")
K0_automatico <- proponi_K(diagnostica_K0)
K0 <- if (params$k_selection_mode  == "automatic") K0_automatico else as.integer(params$K0_manuale)
diagnostica_K0$selezionato <- diagnostica_K0$K == K0
knitr::kable(data.frame(K0_usato = K0, K0_gomito_automatico = K0_automatico,
                        Modalita = params$k_selection_mode ))

tree_alpha_scelto <- hclustgeo(D0_std, D1_std, alpha = optimal_alpha, wt = wt)
diagnostica_K_finale <- curva_gomito(tree_alpha_scelto, alpha = optimal_alpha,
                                     fase = "3. K finale ad alpha scelto")
K_finale_automatico <- proponi_K(diagnostica_K_finale)
K <- kopt <- if (params$k_selection_mode  == "automatic") {
  K_finale_automatico
} else as.integer(params$K_finale_manuale)
diagnostica_K_finale$selezionato <- diagnostica_K_finale$K == K

diagnostiche_K <- dplyr::bind_rows(diagnostica_K0, diagnostica_K_finale)

tagli_selezionati <- diagnostiche_K[diagnostiche_K$selezionato, c("fase", "K")]

# K Selection Elbow Curves
p_selezione_K <- ggplot(diagnostiche_K, aes(K, W_norm)) +
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

ggsave(file.path(output_dir, "Selezione_K_gomiti.png"), p_selezione_K, width = 12, height = 5, dpi = 200, bg = "white")

# Alpha Selection Curves (Qnorm)
dati_alpha_plot <- curva_alpha_K0 %>%
  dplyr::select(alpha, Q0norm, Q1norm, gap_Qnorm) %>%
  pivot_longer(cols = -alpha, names_to = "statistica", values_to = "valore") %>%
  mutate(pannello = if_else(statistica == "gap_Qnorm",
                            "Differenza assoluta: criterio da minimizzare", "Quote di inerzia spiegata normalizzate"))

p_selezione_alpha <- ggplot(dati_alpha_plot, aes(alpha, valore, colour = statistica)) +
  geom_line(linewidth = 0.7) + facet_wrap(~ pannello, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = optimal_alpha, colour = "#5F6572", linetype = 2) +
  scale_colour_manual(values = c(Q0norm = "#365F7C", Q1norm = "#B46C45", gap_Qnorm = "#527D67")) +
  labs(title = paste("Alpha selezionato a K0 =", K0, "-", params$area),
       subtitle = 
         paste(params$alpha_criterion, "- alpha automatico =", "; alpha usato =", 
               format(optimal_alpha, digits = 6)),
       x = "alpha", y = NULL, colour = NULL) + theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Selezione_alpha_Qnorm.png"), p_selezione_alpha, width = 12, height = 5, dpi = 200, bg = "white")

# Dendrograms
disegna_dendrogrammi_selezione <- function() {
  op <- par(mfrow = c(1, 2), mar = c(3, 4, 4, 1))
  on.exit(par(op))
  plot(tree_D0, labels = FALSE, hang = -1, xlab = "", ylab = "Aggregation height",
       main = paste("D0: K0 =", K0), sub = params$area)
  rect.hclust(tree_D0, k = K0, border = "#647B8F")
  plot(tree_alpha_scelto, labels = FALSE, hang = -1, xlab = "", ylab = "Aggregation height",
       main = paste("Misto: K finale =", K), sub = paste("alpha =", optimal_alpha))
  rect.hclust(tree_alpha_scelto, k = K, border = "#647B8F")
}

png(file.path(output_dir, "Selezione_K_dendrogrammi.png"), width = 2400, height = 1200, res = 200)
tryCatch(disegna_dendrogrammi_selezione(), finally = dev.off())
