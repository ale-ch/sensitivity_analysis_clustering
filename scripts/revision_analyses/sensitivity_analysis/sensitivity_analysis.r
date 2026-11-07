##########################################################
########## Spatial cluster coherence analysis ############
##########################################################


library(sf)
library(tidyverse)
library(readxl)
library(pracma)
library(cluster)
library(stats)
library(fpc)
library(factoextra) 
library(NbClust)
library(mclust)

base_path <- '/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering'
setwd(base_path)

source("scripts/utils/make_dataset.r")



# Define baseline parameters based on text
k_opt <- 8
alpha_opt <- 0.10

# Define nearby parameter spaces
alphas_nearby <- c(0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16)
ks_nearby <- c(k_opt - 2, k_opt - 1, k_opt, k_opt + 1, k_opt + 2)

# Generate baseline partition
tree_baseline <- hclustgeo(D0, D1, alpha = alpha_opt, wt = rep(1/n, n))
part_baseline <- cutree(tree_baseline, k_opt)

# Initialize results dataframe
ari_results <- expand.grid(Alpha = alphas_nearby, K = ks_nearby)
ari_results$ARI <- NA

# Run iterations to compare partitions
for (i in 1:nrow(ari_results)) {
  a <- ari_results$Alpha[i]
  k <- ari_results$K[i]
  
  # Generate tree and partition for the modified parameters
  temp_tree <- hclustgeo(D0, D1, alpha = a, wt = rep(1/n, n))
  temp_part <- cutree(temp_tree, k)
  
  # Calculate Adjusted Rand Index compared to baseline (1 = identical)
  ari_results$ARI[i] <- adjustedRandIndex(part_baseline, temp_part)
}

# Display numeric sensitivity results
print(ari_results)

# Visualize sensitivity map
#ggplot(ari_results, aes(x = factor(Alpha), y = factor(K), fill = ARI)) +
#  geom_tile(color = "white") +
#  scale_fill_viridis_c(limits = c(0, 1), option = "mako", name = "Adjusted Rand Index") +
#  geom_text(aes(label = round(ARI, 2)), color = "white") +
#  labs(title = "Clustering Sensitivity to Nearby K and Alpha",
#       x = "Mixing Parameter (Alpha)",
#       y = "Number of Clusters (K)") +
#  theme_bw()




ggplot(ari_results, aes(x = factor(Alpha), y = factor(K), fill = ARI)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(
    label = sprintf("%.2f", ARI), 
    colour = ifelse(ARI > 0.65, "white", "#243746")
  ), size = 3) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#B97572", mid = "#F3F2EF", high = "#365F7C",
    midpoint = 0.5, limits = c(0, 1), name = "Adjusted Rand Index"
  ) +
  coord_equal() +
  labs(
    title = "Clustering Sensitivity to Nearby K and Alpha",
    x = "Mixing Parameter (Alpha)",
    y = "Number of Clusters (K)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(), 
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", colour = "#183746")
  )


# 3. Save the plots
ggsave(file.path(output_dir, "ARI_compare_heatmap.png"),
       p_heatmap, width = 10, height = 12, dpi = 160, bg = "white")
ggsave(file.path(output_dir, "ARI_compare_heatmap.pdf"),
       p_heatmap, width = 10, height = 12)
