##########################################################
########## Spatial cluster coherence analysis ############
##########################################################

base_path <- '/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering'
setwd(base_path)

source("scripts/revision_analyses/make_dataset.r")


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

