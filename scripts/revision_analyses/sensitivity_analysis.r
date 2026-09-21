##########################################################
########## Spatial cluster coherence analysis ############
##########################################################

base_path <- '/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering'
setwd(base_path)

source("scripts/revision_analyses/make_dataset.r")



# Define baseline parameters based on text
k_opt <- 8
alpha_opt <- 0.10

# Define nearby parameter spaces
alphas_nearby <- c(0.05, 0.08, 0.10, 0.12, 0.15)
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
ggplot(ari_results, aes(x = factor(Alpha), y = factor(K), fill = ARI)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(limits = c(0, 1), option = "mako", name = "Adjusted Rand Index") +
  geom_text(aes(label = round(ARI, 2)), color = "white") +
  labs(title = "Clustering Sensitivity to Nearby K and Alpha",
       x = "Mixing Parameter (Alpha)",
       y = "Number of Clusters (K)") +
  theme_bw()
