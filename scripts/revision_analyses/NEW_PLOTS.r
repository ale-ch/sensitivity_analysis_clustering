library(tidyverse)
library(sf)

load('/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/Output_Selezione_K_alpha_RMarkdown/Lombardia_senza_Milano_20261109_164347/env/Res_national_11082024.RData')

#### CLUSTERING MAPS (4-panels) ####

crea_mappa_multik2 <- function(prefisso, titolo, sottotitolo_plot) {
  campi <- paste0(prefisso, c(3, 5, 8, 15))
  dati_map <- DPCM_red2 %>% dplyr::select(all_of(campi)) %>%
    pivot_longer(cols = all_of(campi), names_to = "Taglio", values_to = "Clusters") %>%
    mutate(Taglio = factor(Taglio, levels = campi, labels = c("K = 3", "K = 5", "K = 8", "K = 15")),
           Clusters = factor(Clusters, levels = 1:15))
  ggplot(dati_map) + geom_sf(aes(fill = Clusters), linewidth = 0.035, colour = "#F5F3EF") +
    facet_wrap(~ Taglio, ncol = 2) +
    scale_fill_manual(values = palette_multik, name = "Cluster", drop = FALSE) +
    labs(title = titolo, subtitle = sottotitolo_plot, x = "", y = "",
         #caption = "Tagli dello stesso albero; etichette native, senza abbinamento fra diversi K."
         ) +
    tema_originale + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
}


plot_spatial <-crea_mappa_multik2("D0D1_P", "Geographically-informed Ward hierarchical clustering",
                         paste("Mixing parameter of the geographical distance matrix α = 0.67"))
ggpubr::ggexport(plot_spatial, width = 1800, height = 1200, res = 150, filename = "GeoInform.png")


plot_non_spatial <- crea_mappa_multik2("D0_P", "Ward-like hierarchical clustering",
                   paste("Mixing parameter of the geographical distance matrix α = 0"))
ggpubr::ggexport(plot_non_spatial, width = 1800, height = 1200, res = 150, filename = "NonGeoInform.png")

#### DENDROGRAM ####

png("dendrogram_plot.png", width = 800, height = 600, res = 100)

par(mar = c(5, 4, 4, 1)) 
plot(treeD0D1, hang = -1, labels = FALSE, xlab = "", ylab = "Aggregation height",
     main = "Ward dendrogram with both DPCM and geographical features", 
     sub = paste("Mixing parameter of the geographical distance matrix α =", alpha_ottimo))

gruppi <- cutree(treeD0D1, 5L)
ordine <- unique(gruppi[treeD0D1$order])
rect.hclust(treeD0D1, k = 5L, border = unname(palette_multik[as.character(ordine)]))
legend("topright", legend = paste0("C", seq_len(5L)),
       fill = unname(palette_multik[as.character(seq_len(5L))]), border = NA, bty = "n", cex = 0.7)

dev.off()



#### SELECTION OF K ####
library(dplyr)

diagnostiche_K <- diagnostiche_K %>%
  mutate(fase = case_when(
    fase == "1. K0 su D0" ~ "Optimal K on socioeconomic distance matrix",
    fase == "3. K finale ad alpha scelto" ~ "Optimal K at α = 0.67",
    TRUE ~ fase
  ))

tagli_selezionati <- tagli_selezionati %>%
  mutate(fase = case_when(
    fase == "1. K0 su D0" ~ "Optimal K on socioeconomic distance matrix",
    fase == "3. K finale ad alpha scelto" ~ "Optimal K at α = 0.67",
    TRUE ~ fase
  ))



elbow_plot <- ggplot(diagnostiche_K, aes(K, W_norm)) +
  geom_segment(aes(x = 1, xend = max(K_candidati), y = 1, yend = 0),
               data = tagli_selezionati, inherit.aes = FALSE,
               linetype = 3, colour = "#9AA3A9") +
  geom_line(colour = "#365F7C", linewidth = 0.7) +
  geom_point(aes(shape = candidato), colour = "#365F7C", size = 1.8) +
  geom_vline(data = tagli_selezionati, aes(xintercept = K), colour = "#B46C45", linetype = 2) +
  facet_wrap(~ fase, ncol = 2) +
  scale_x_continuous(breaks = seq_len(max(K_candidati))) +
  scale_shape_manual(values = c(1, 16), guide = "none") +
  labs(title = paste("Choice of optimal number of clusters K - Lombardy"),
       # subtitle = "Linea tratteggiata: K usato; corda punteggiata: riferimento del gomito",
       x = "K", y = "Rescaled internal intertia") +
  theme_minimal(base_size = 11)

ggpubr::ggexport(elbow_plot, width = 1800, height = 1200, res = 150, filename = "elbow_plot.png")
