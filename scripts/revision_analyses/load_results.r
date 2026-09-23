library(tidyverse)
library(sf)

load("/Users/alessiochiodin/Downloads/res1911.RData")

Data_clusters %>%
  select(# Cluster8_NoGeo,
    Cluster8_Geo,
    D0D1_P8) %>%
  mutate(Cluster8_Geo = as.factor(Cluster8_Geo)) %>%
  filter(Cluster8_Geo %notin% c("Regional median")) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Cluster8_Geo)) + 
  labs(x = "",
       y = "",
       title = "Geographically-informed Ward hierarchical clustering",
       # subtitle = latex2exp::TeX("K = 8 clusters + Milano (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0.10)")
       subtitle = latex2exp::TeX("K = 8 clusters (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0.10)")
  ) +
  scale_fill_manual("", values = colp) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 16),
        plot.subtitle = element_text(face = "bold",size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) + 
  guides(col=guide_legend(nrow=2,byrow=TRUE),
         linetype=guide_legend(nrow=2,byrow=TRUE))
