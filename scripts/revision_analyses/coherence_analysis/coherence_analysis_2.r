library(tidyverse)

# load('/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/Output_Selezione_K_alpha_RMarkdown/Lombardia_senza_Milano_20261107_104128/env/Res_national_11082024.RData')
load("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/Output_Selezione_K_alpha_RMarkdown/Lombardia_senza_Milano_20261109_164347/env/Res_national_11082024.RData")

etichette_ARI_stabilita <- c("α = 0",
                             paste0("optimum α = ", format(alpha_ottimo, digits = 6)))
dati_ARI_stabilita_grafico <- ARI_stabilita_K_long %>%
  mutate(Pannello = factor(caso, levels = c("alpha0", "alpha_ottimo"),
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
  labs(title = paste("ARI between values of K - Lombardy"),
       subtitle = paste(n, "municipalities; constant alpha on each matrix"),
       x = "K (columns)", y = "K (rows)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", colour = "#183746"),
        plot.title = element_text(face = "bold", colour = "#183746"))
p_ARI_stabilita



ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/spatial_vs_nonspatial_SENZA_MILANO", "ARI_stabilita_K_mappe_calore.png"),
       p_ARI_stabilita, width = 14, height = 7.5, dpi = 200, bg = "white")
ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/spatial_vs_nonspatial_SENZA_MILANO", "ARI_stabilita_K_mappe_calore.pdf"),
       p_ARI_stabilita, width = 14, height = 7.5)