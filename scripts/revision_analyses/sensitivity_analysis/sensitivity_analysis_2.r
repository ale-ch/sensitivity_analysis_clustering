
load('/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/Output_Selezione_K_alpha_RMarkdown/Lombardia_senza_Milano_20261107_104128/env/Res_national_11082024.RData')

str(dati_ARI_alpha_grafico)

summarized_ARI_data <- dati_ARI_alpha_grafico %>% 
  filter(ARI != 1) %>% 
  group_by(K) %>% 
  mutate(
    ari_at_optim_alpha = if_else(alpha_riga == 0.1 | alpha_riga == 0.1, ARI, NA)
  ) %>% 
  reframe(
    #min(ari_at_optim_alpha, na.rm = TRUE),
    #max(ari_at_optim_alpha, na.rm = TRUE),
    mean_ari_at_optim_alpha = mean(ari_at_optim_alpha, na.rm = TRUE),
    #median(ari_at_optim_alpha, na.rm = TRUE),
    sd_ari_at_optim_alpha = sd(ari_at_optim_alpha, na.rm = TRUE),
    min = min(ARI),
    max = max(ARI),
    mean = mean(ARI),
    median = median(ARI),
    sd = sd(ARI)
  ) %>%
  rename(
    "Mean ARI (α = 0.1)" = mean_ari_at_optim_alpha,
    "SD ARI (α = 0.1)"   = sd_ari_at_optim_alpha,
    "Minimum"                = min,
    "Maximum"                = max,
    "Mean"                   = mean,
    "Median"                 = median,
    "Std Dev"                = sd
  )

str(summarized_ARI_data)





library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Reshape the summary data into long format for the heatmap
dati_heatmap <- summarized_ARI_data %>%
  pivot_longer(
    cols = -K, 
    names_to = "Statistic", 
    values_to = "Value"
  ) %>%
  mutate(
    # Format K as a factor and reverse it so lowest K is at the top of the y-axis
    K_factor = factor(K, levels = rev(sort(unique(K)))),
    # Clean up statistic names for the plot labels
    Statistic = factor(Statistic, levels = colnames(summarized_ARI_data)[-1])
  )

# 2. Generate the heatmap adapting your styling
p_heatmap <- ggplot(dati_heatmap, aes(x = Statistic, y = K_factor, fill = Value)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(
    label = sprintf("%.3f", Value),
    # Adapt text color based on the value for readability against the background
    colour = ifelse(Value > 0.65 | Value < -0.65, "white", "#243746") 
  ), size = 3) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#B97572", mid = "#F3F2EF", high = "#365F7C",
    midpoint = 0, limits = c(-0.1, 1), name = "Value" # Limits adjusted to fit your min/max stats
  ) +
  coord_equal() +
  labs(
    title = "ARI Summary Statistics by K",
    subtitle = "Comparisons among alpha values between 0 and 1, step size 0.1",
    x = "Summary Statistic", 
    y = "K"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", colour = "#183746")
  )

p_heatmap

# 3. Save the plots
ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/sensitivity_analysis_SENZA_MILANO", "ARI_summary_heatmap.png"),
       p_heatmap, width = 10, height = 12, dpi = 160, bg = "white")
ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/sensitivity_analysis_SENZA_MILANO", "ARI_summary_heatmap.pdf"),
       p_heatmap, width = 10, height = 12)