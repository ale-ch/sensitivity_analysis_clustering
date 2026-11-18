# Correlation Matrix Plot
corr <- CData[, -c(1:6)] %>% sf::st_drop_geometry() %>% cor() %>%
  reshape2::melt(na.rm = TRUE) %>% rename(Reg1 = Var1, Reg2 = Var2, Corr2021 = value)
corrmat <- corr %>% pivot_longer(cols = 3:last_col()) %>%
  group_by(name) %>%
  mutate(value_std = scale(x = value, center = TRUE, scale = TRUE)) %>% ungroup()


p_corr <- ggplot(corrmat, aes(x = Reg2, y = Reg1, fill = as.numeric(value_std))) +
  geom_tile(colour = "white") + coord_fixed() +
  geom_text(aes(label = round(value, 2)), colour = "#263C49", size = 2.6) +
  scale_fill_gradient2(low = "#497D9A", mid = "#F2F1EE", high = "#B56F6D", name = "Measure") +
  labs(title = "Linear correlation among DPCM-related variables", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

ggsave(file.path(output_dir, "CorrMat.png"), p_corr, width = 2400/150, height = 2000/150, dpi = 150, bg = "white")

# Distance Matrices Heatmaps
f1 <- fviz_dist(DMat, gradient = list(low = "white", mid = "grey", high = "black"), show_labels = FALSE)
ggsave(file.path(output_dir, "D0_mat.png"), f1, width = 2400/150, height = 2000/150, dpi = 150, bg = "white")

f2 <- fviz_dist(DMat_geo, gradient = list(low = "white", mid = "grey", high = "black"), show_labels = FALSE)
ggsave(file.path(output_dir, "D1_mat.png"), f2, width = 2400/150, height = 2000/150, dpi = 150, bg = "white")

# Average Distance Maps
p_avg0 <- ggplot(DPCM_red2) + geom_sf(aes(fill = AvgDist_D0), colour = NA) +
  scale_fill_gradientn(colours = colori_continui, name = "Eucl. Dist.") +
  labs(title = paste("Average distance on DPCM variables -", params$area), 
       # subtitle = sottotitolo, 
       x = "", y = "") + tema_originale

ggsave(file.path(output_dir, "AvgDist_D0.png"), p_avg0, width = 1800/150, height = 1200/150, dpi = 150, bg = "white")

p_avg1 <- ggplot(DPCM_red2) + geom_sf(aes(fill = AvgDist_D1), colour = NA) +
  scale_fill_gradientn(colours = colori_continui, name = "Meters") +
  labs(title = paste("Average geographical distance -", params$area),
       subtitle = "WGS84 geodetic distance between clustering centroids", x = "", y = "") + tema_originale
p_avg1

ggsave(file.path(output_dir, "AvgDist_D1.png"), p_avg1, width = 1800/150, height = 1200/150, dpi = 150, bg = "white")