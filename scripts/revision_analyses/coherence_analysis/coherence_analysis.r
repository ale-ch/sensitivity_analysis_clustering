##########################################################
########## Spatial cluster coherence analysis ############
##########################################################

base_path <- '/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering'
setwd(base_path)

source("scripts/utils/make_dataset.r")


##### Ward-like hierarchical clustering without spatial information
tree_D0 <- hclustgeo(D0 = D0, wt=rep(1/n,n))
P8 <- cutree(tree_D0, 8)


##### Ward-like hierarchical clustering with spatial information
# D1_std <- D1/max(D1)
treeD0D1 <- hclustgeo(D0, D1, alpha=0.1, wt = rep(1/n, n))
D0D1_P8 <- cutree(treeD0D1, 8)


CData$clusters_alpha_0 <- P8
CData$clusters_alpha_01 <- D0D1_P8


# Convert dataframe to a spatial object and project to meters (EPSG: 32632 for Italy)
CData_sf <- st_as_sf(CData, coords = c("Lon_centr_LAU", "Lat_centr_LAU"), crs = 4326) %>%
  st_transform(32632)

distances_geo <- numeric(8)
distances_no_geo <- numeric(8)

for(i in 1:8) {
  # Filter by cluster
  cluster_geo <- CData_sf %>% filter(clusters_alpha_01 == i)
  cluster_nogeo <- CData_sf %>% filter(clusters_alpha_0 == i)
  
  # Calculate distance matrices in meters and find the mean
  distances_geo[i] <- mean(as.vector(st_distance(cluster_geo)), na.rm = TRUE)
  distances_no_geo[i] <- mean(as.vector(st_distance(cluster_nogeo)), na.rm = TRUE)
}


summary(distances_geo)
summary(distances_no_geo)

sd(distances_geo)
sd(distances_no_geo)




# Compile statistics into a single dataframe
summary_table <- data.frame(
  Model = c("Alpha = 0 (No Spatial Constraint)", "Alpha = 0.1 (With Spatial Constraint)"),
  Min = c(min(distances_no_geo[distances_no_geo > 0]), min(distances_geo[distances_geo > 0])),
  Q1 = c(quantile(distances_no_geo, 0.25), quantile(distances_geo, 0.25)),
  Median = c(median(distances_no_geo), median(distances_geo)),
  Mean = c(mean(distances_no_geo), mean(distances_geo)),
  Q3 = c(quantile(distances_no_geo, 0.75), quantile(distances_geo, 0.75)),
  Max = c(max(distances_no_geo), max(distances_geo)),
  SD = c(sd(distances_no_geo), sd(distances_geo))
)

# Convert from meters (EPSG: 32632) to kilometers for readability
summary_table_km <- summary_table %>%
  mutate(across(where(is.numeric), ~ .x / 1000))

# Generate a clean Markdown table ready for copy-pasting
knitr::kable(
  summary_table_km, 
  digits = 2, 
  col.names = c("Model", "Min (km)", "Q1 (km)", "Median (km)", "Mean (km)", "Q3 (km)", "Max (km)", "SD (km)"),
  caption = "Spatial Cluster Coherence: Intra-cluster Distances"
)

