load('/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/Output_Selezione_K_alpha_RMarkdown/Lombardia_senza_Milano_20261107_104128/env/Res_national_11082024.RData')

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# 1. Combine the curve data from the list of results
# Assuming curve_per_K is populated by the loop in the provided code
dati_curve_all_K <- map_df(seq_along(curve_per_K), function(i) {
  curva <- curve_per_K[[i]]
  if (is.null(curva)) return(NULL)
  
  data.frame(
    
    K = i,
    alpha = curva$range.alpha,
    Q0 = curva$Q[, 1],
    Q1 = curva$Q[, 2],
    Q0norm = curva$Qnorm[, 1],
    Q1norm = curva$Qnorm[, 2]
  )
})

# 2. Reshape to long format for plotting
dati_curve_long <- dati_curve_all_K %>%
  pivot_longer(
    cols = c(Q0, Q1, Q0norm, Q1norm),
    names_to = "Stat", 
    values_to = "Value"
  ) %>%
  mutate(
    Value = Value * 100, # Convert to percentage
    Grp = if_else(Stat %in% c("Q0", "Q1"), "Non-normalized pseudo inertia", "Normalized pseudo inertia"),
    Stat = factor(Stat, levels = c("Q0", "Q1", "Q0norm", "Q1norm")),
    K_factor = as.factor(K)
  ) %>%
  # Filter out NA values which might occur for Qnorm when K=1
  filter(is.finite(Value)) 

# 3. Create the plot with overlaid curves for each K
cols <- c(Q0 = "#3475AB", Q1 = "#B9474E", Q0norm = "#3E9378", Q1norm = "#A28E3D")


pastello <- c("#8BA9C4", "#D4AA84", "#9BB7A6", "#BAA7C4",
              "#C99196", "#94B9BC", "#C3B984", "#A6ABB3",
              "#A3B2C6", "#C8B99B", "#AFB995", "#C3A8B2", 
              "#AAC5C2", "#C6AB91", "#94AFC4")


# Calculate the intersection alpha for each K, then find the mean for each group
intersezioni <- dati_curve_all_K %>%
  filter(K != 1) %>% 
  group_by(K) %>%
  summarise(
    # Interpolate the alpha where Q0 - Q1 = 0
    alpha_int_non_norm = approx(x = Q0 - Q1, y = alpha, xout = 0)$y,
    # Interpolate the alpha where Q0norm - Q1norm = 0
    alpha_int_norm = if (any(is.finite(Q0norm) & is.finite(Q1norm))) {
      approx(x = Q0norm - Q1norm, y = alpha, xout = 0)$y
    } else {
      NA_real_
    }
  ) %>%
  ungroup()

# Calculate the mean intersection alpha for the two groups
mean_intersezioni <- data.frame(
  Grp = c("Non-normalized pseudo inertia", "Normalized pseudo inertia"),
  mean_alpha = c(
    mean(intersezioni$alpha_int_non_norm, na.rm = TRUE),
    mean(intersezioni$alpha_int_norm, na.rm = TRUE)
  )
)

# Generate the plot with the vertical line
p_choice_all_K <- ggplot(dati_curve_long, 
                         aes(x = alpha, y = Value, 
                             colour = K_factor, 
                             linetype = Stat, 
                             group = interaction(Stat, K_factor))) +
  geom_line(linewidth = 0.7) + 
  facet_wrap(~ Grp, scales = "free_y") +
  geom_vline(data = mean_intersezioni, aes(xintercept = mean_alpha), 
             colour = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = mean_intersezioni, 
            aes(x = mean_alpha, y = Inf, label = sprintf("Mean optimal α = %.3f", mean_alpha)),
            colour = "red", vjust = 2, hjust = -0.1, 
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_colour_manual(values = pastello) +
  # Manually assign line types: Q1/Q1norm as solid, Q0/Q0norm as dotted
  scale_linetype_manual(values = c("Q0" = "dotted", "Q1" = "solid", 
                                   "Q0norm" = "dotted", "Q1norm" = "solid")) +
  labs(
    title = paste("Explained pseudo-inertia curves, by value of K - Lombardy"),
    x = "Alpha", 
    y = "% of explained pseudo inertia", 
    colour = "Number of Clusters (K)",
    linetype = "Statistic"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", colour = "#183746"),
    plot.title = element_text(face = "bold", colour = "#183746")
  )

p_choice_all_K






##############

# Calculate summary lines (mean of each statistic across all K values at each alpha)
dati_summary <- dati_curve_long %>%
  group_by(alpha, Stat, Grp) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# Generate the plot
p_choice_all_K <- ggplot() +
  # Individual K curves (slightly transparent so the summary lines stand out)
  geom_line(data = dati_curve_long, 
            aes(x = alpha, y = Value, 
                colour = K_factor, 
                linetype = Stat, 
                group = interaction(Stat, K_factor)),
            linewidth = 0.5, alpha = 0.6) + 
  # Summary curves (thicker, black)
  geom_line(data = dati_summary, 
            aes(x = alpha, y = Value, linetype = Stat, group = Stat),
            colour = "blue", linewidth = 1.2) +
  # Stack facets vertically to fix horizontal squishing
  facet_wrap(~ Grp, scales = "free_y", ncol = 1) +
  geom_vline(data = mean_intersezioni, aes(xintercept = mean_alpha), 
             colour = "red", linetype = "dashed", linewidth = 0.8) +
  geom_text(data = mean_intersezioni, 
            aes(x = mean_alpha, y = Inf, label = sprintf("Mean optimal α = %.3f", mean_alpha)),
            colour = "red", vjust = 2, hjust = -0.1, 
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_colour_manual(values = pastello) +
  scale_linetype_manual(values = c("Q0" = "dotted", "Q1" = "solid", 
                                   "Q0norm" = "dotted", "Q1norm" = "solid")) +
  labs(
    title = paste("Explained pseudo-inertia curves, by value of K - Lombardy"),
    x = "Alpha", 
    y = "% of explained pseudo inertia", 
    colour = "Number of Clusters (K)",
    linetype = "Statistic"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", colour = "#183746"),
    plot.title = element_text(face = "bold", colour = "#183746")
  )

p_choice_all_K


# 3. Save the plots
ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/sensitivity_analysis_SENZA_MILANO", "pseudoinertia_curves.png"),
       p_choice_all_K, width = 10, height = 12, dpi = 160, bg = "white")
ggsave(file.path("/Volumes/T7 Shield/FRES/work/sensitivity_analysis_clustering/reply/sensitivity_analysis_SENZA_MILANO", "pseudoinertia_curves.pdf"),
       p_choice_all_K, width = 10, height = 12)
