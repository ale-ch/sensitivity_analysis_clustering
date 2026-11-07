##############################################
########## Spatial cluster analysis ##########
##############################################

setwd("K:/Il mio Drive/FEEM_COVID19_Lombardia/Analisi/SpatialCluster_Fragility")

library(sf)
library(tidyverse)
library(readxl)
library(pracma)
library(cluster)
library(stats)
library(fpc)
library(factoextra) 
library(NbClust)

'%notin%' <- Negate('%in%')

# AOI <- "Italy"
AOI <- "Lombardy"

source("AuxFuns_NbClust.R", encoding = 'UTF-8')

### Funzione del kmeans++
kmpp <- function(X, k) {
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  
  for (i in 2:k) {
    dm <- distmat(X, X[C, ])
    pr <- apply(dm, 1, min); pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  kmeans(X, X[C, ], nstart=50)
}

##### DPCM data
DPCM <- read_excel("K:/Il mio Drive/FEEM_COVID19_Lombardia/Dati/DataManagement/Original data/comune/DPCM/DPCM_cleaned.xlsx")

##### Municipalities registry
ManualMatchLAUs <- read_excel("K:/Il mio Drive/FEEM_COVID19_Lombardia/Dati/DataManagement/Original data/anagrafica comuni/Manually_matched_LAUs_Oct2023.xlsx", 
                              sheet = "Modified")
ManualMatchLAUs <- ManualMatchLAUs %>%
  mutate(across(c(`Denominazione Comune NEW`), toupper),
         across(c(`Denominazione Comune NEW`), stringr::str_to_title)) %>%
  dplyr::select(LAU_Code_old = `Codice Istat Comune`,
                LAU_Code_new = `Codice Istat Comune NEW`) %>%
  distinct()

# load("../../Dati/DataManagement/Data transformed/Unnormalized/Data_CrossSection_C2.RData")
load("../../Dati/DataManagement/Original data/anagrafica comuni/LAU2021_Eurostat_Italy.RData")

DPCM_red <- DPCM %>%
  dplyr::rename(LAU_Code = ISTAT_Code_LAU) %>%
  dplyr::select(-Name_LAU)
# DPCM_red <- data %>%
#   dplyr::rename(LAU_Code = AUX_IstatMunicipalityCode)

codes <- stringr::str_pad(string = ManualMatchLAUs$LAU_Code_old,width = 6,pad = "0")
codes_new <- stringr::str_pad(string = ManualMatchLAUs$LAU_Code_new,width = 6,pad = "0")
for (r in 1:dim(DPCM_red)[1]) {
  pos <- which(DPCM_red$LAU_Code[r] == codes)
  if (!is_empty(pos)) {
    cat(paste0("Row ",r," has a different code\n"))
    DPCM_red[r,]$LAU_Code <- codes_new[pos]
  }
}
DPCM_red <- DPCM_red %>%
  group_by(LAU_Code) %>%
  mutate(across(is.numeric, function(x) sum(x,na.rm=T))) %>%
  ungroup()


DPCM_red <- left_join(x = LAU_IT_sf, y = DPCM_red, by = c("LAU_Code"))

png("SuspendedDPCM.png", width = 1800, height = 1200, res = 150)
DPCM_red %>%
  mutate('Share of suspended local units in industry' = Local_Units_Industry_Suspended/Local_Units_Industry_Total*100,
         'Share of suspended local units in services' = Local_Units_Services_Suspended/Local_Units_Services_Total*100) %>%
  select('Share of suspended local units in industry','Share of suspended local units in services') %>%
  pivot_longer(cols = c('Share of suspended local units in industry',
                        'Share of suspended local units in services'),
               names_to = "Var", values_to = "Values") %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Values)) + 
  scale_fill_viridis_b(option = "plasma") + 
  facet_wrap(~ Var) + 
  theme_bw()
dev.off()

DPCM_red %>%
  filter(NUTS2_Name %in% c("Lombardia")) %>%
  select(Local_Units_Industry_Suspended,Local_Units_Services_Suspended) %>%
  pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended),
               names_to = "Var", values_to = "Values") %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Values)) + 
  scale_fill_viridis_b(option = "plasma",n.breaks = 100) + 
  facet_wrap(~ Var) + 
  theme_bw()


Data_plot <- DPCM_red %>%
  # filter(NUTS2_Name %in% c("Lombardia")) %>%
  select(Area_km2,Local_Units_Industry_Suspended,Value_Added_Industry_Suspended,Employed_Industry_Suspended,
         Local_Units_Services_Suspended,Value_Added_Services_Suspended,Employed_Services_Suspended) %>%
  mutate(across(is.numeric, function(x) x/Area_km2)) %>%
  select(-c(Area_km2)) %>%
  ### Normalization w.r.t. the maximum
  # mutate(across(is.numeric, function(x) x/max(x,na.rm = T))) %>%
  ### Min-Max (0,1) normalization
  # mutate(across(is.numeric, function(x) (x - min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))) %>%
  ### Robust normalization
  mutate(across(is.numeric, function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T)))) %>%
  pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended,
                        Value_Added_Industry_Suspended,Value_Added_Services_Suspended,
                        Employed_Industry_Suspended,Employed_Services_Suspended),
               names_to = "Vars", values_to = "Values") %>%
  mutate(
    Vars = case_when(
      Vars %in% c("Local_Units_Industry_Suspended") ~ "Industry: local units suspended",
      Vars %in% c("Value_Added_Industry_Suspended") ~ "Industry: value adde of local units suspended",
      Vars %in% c("Employed_Industry_Suspended") ~ "Industry: employed persons suspended",
      Vars %in% c("Local_Units_Services_Suspended") ~ "Services: local units suspended",
      Vars %in% c("Value_Added_Services_Suspended") ~ "Services: value adde of local units suspended",
      Vars %in% c("Employed_Services_Suspended") ~ "Services: employed persons suspended",
      TRUE ~ NA_character_
    ),
    Vars = factor(
      Vars,
      levels = c("Industry: local units suspended","Industry: value adde of local units suspended",
                 "Industry: employed persons suspended","Services: local units suspended",
                 "Services: value adde of local units suspended","Services: employed persons suspended"),
      labels = c("Industry: local units suspended","Industry: value adde of local units suspended",
                 "Industry: employed persons suspended","Services: local units suspended",
                 "Services: value adde of local units suspended","Services: employed persons suspended"),
      ordered = T)
  )
q <- quantile(Data_plot$Values, probs = seq(0.1,1,by=0.03),na.rm = T)

png("NormVariablesD0_Milano.png", width = 1800, height = 1200, res = 150)
Data_plot %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Values)) + 
  scale_fill_stepsn(colours = colorRampPalette(c("green", "yellow", "red"))(length(q)),
                    breaks = round(unname(q),2),
                    limits = c(min(round(unname(q),2)),max(round(unname(q),2))),
                    values = scales::rescale(round(unname(q),2),
                                             from = c(min(round(unname(q),2)),max(round(unname(q),2))))
  ) +
  facet_wrap(~ Vars, nrow = 2, ncol = 3) + 
  labs(# title = "Normalized DPCM-related variables by industry and services",
       # subtitle = "Here, Milan is used as the maximum value for normalizing",
       title = "Double-normalized DPCM-related variables by industry and services",
       subtitle = "Robust normalization (using quantiles) on variables per km2",
       y = "", x = "") + 
  theme_bw() + 
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        legend.key.height = unit(25, 'mm'), legend.key.width = unit(5, 'mm'))
dev.off()








#############################################
########## Define area of interest ##########
#############################################
if (AOI == "Lombardy") {
  DPCM_red <- DPCM_red %>%
    filter(NUTS2_Name == "Lombardia")
}

apply(DPCM_red, MARGIN = 2, FUN = function(x) sum(is.na(x)))

DPCM_red %>%
  mutate(var = scale(Local_Units_Industry_Suspended/Local_Units_Industry_Total*100,
                     scale = T,center=T)) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = var)) + 
  scale_fill_viridis_b(option = "plasma")



###########################################################
########## Spatially-informed Kmeans++ Algorithm ##########
###########################################################
CData <- DPCM_red %>%
  # Filter non-missing data
  filter(!is.na(Local_Units_Industry_Active)) %>%
  # Select numeric variables
  dplyr::select(Area_km2,LAU_Code,Lon_centr_LAU,Lat_centr_LAU,
                Alt_mean,Alt_min,Alt_max,
                contains("Suspended"),contains("Active")) %>%
  mutate(across(c(is.numeric,
                  -contains("Alt"),
                  -contains("Lon"),-contains("Lat")), function(x) x/Area_km2)) %>%
  select(-c(Area_km2)) %>%
  ### Normalization w.r.t. the maximum
  # mutate(across(is.numeric, function(x) x/max(x,na.rm = T))) %>%
  ### Min-Max (0,1) normalization
  # mutate(across(is.numeric, function(x) (x - min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))) %>%
  ### Standardization
  # mutate(across(c(is.numeric,
  #                 -contains("Lon"),-contains("Lat")), function(x) scale(x)[,1]))
  ### Robust normalization
  mutate(across(c(is.numeric,
                  -contains("Lon"),-contains("Lat")), function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T))))

st_geometry(CData) <- NULL
Xmat_geo <- as.matrix(CData[,2:3])
Xmat_nogeo <- as.matrix(CData[,-c(1:6)])



##### Distance matrix among units: socio-economic variables
DMat <- get_dist(Xmat_nogeo)
# f1 <- fviz_dist(DMat, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
# f1
DMat_geo <- get_dist(Xmat_geo)
# f2 <- fviz_dist(DMat_geo, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
# f2

AvgDist_D0 <- rowMeans(as.matrix(DMat))

DMat_geo_k2 <- mgcv::slanczos(A = as.matrix(DMat_geo),k = 2)

Xmat_full <- cbind(DMat_geo_k2$vectors,Xmat_nogeo)
DMat_full <- get_dist(Xmat_full)


########## Not spatially-informed K-means++ 
clustering <- list()
clust <- seq(from=2,to=2)
indexes <- matrix(NA,length(clust),3)
set.seed(12345)

for (cl in 1:length(clust)){
  cat(paste0("Cluster ",cl," of ",length(clust),"\n"))
  clustering[[cl]] <- kmpp(X = Xmat_nogeo, k = clust[cl])
  stats <- cluster.stats(DMat, clustering[[cl]]$cl)
  indexes[cl,] <- c(stats$avg.silwidth, stats$pearsongamma,stats$within.cluster.ss)
}
# Performance indices
colnames(indexes) <- c("ASW", "PG","WSS")
rownames(indexes) <- paste0("cl = ",clust)
indexes
# Plot
win.graph()
par(mfrow = c(3,1))
plot(clust,indexes[,1], type = "b")
plot(clust,indexes[,2], type = "b")
plot(clust,indexes[,3], type = "b")
par(mfrow = c(3,1))








res_kmpp_nogeo <- NbClust_mod(data = Xmat_nogeo,
                              distance = "euclidean", min.nc = 2, max.nc = 30,
                              method = "kmpp", 
                              index ="all")
attr(res_kmpp_nogeo$Best.nc, "class") <- c("matrix")
res_kmpp_geo <- NbClust_mod(data = Xmat_full,
                            distance = "euclidean", min.nc = 2, max.nc = 10,
                            method = "kmpp", 
                            index ="all")
attr(res_kmpp_geo$Best.nc, "class") <- c("matrix")


## "Regola della maggioranza"
p1 <- fviz_nbclust(Xmat, FUNcluster = kmpp, method = "gap_stat", nboot = 2, k.max = 30) +
  labs(subtitle = "Gap statistic",
       title = "K-means++ with geographic coordinates",
       x = "# clusters")
p2 <- fviz_nbclust(Xmat, FUNcluster = kmpp, method = "silhouette", nboot = 2, k.max = 30) + 
  labs(subtitle = "Silhouette",
       title = "K-means++ with geographic coordinates",
       x = "# clusters")
p3 <- fviz_nbclust(Xmat, FUNcluster = kmpp, method = "wss", nboot = 2, k.max = 30) + 
  labs(subtitle = "Total Within Sum of Squares",
       title = "K-means++ with geographic coordinates",
       x = "# clusters")

p3 <- fviz_nbclust(x = res_kmpp_geo, method = "silhouette") +
  labs(subtitle = "Total Within Sum of Squares",
       title = "K-means++ with geographic coordinates",
       x = "# clusters")
p3 <- fviz_nbclust(x = res_kmpp_nogeo) +
  labs(subtitle = "Total Within Sum of Squares",
       title = "K-means++ with geographic coordinates",
       x = "# clusters")

class(nb_kmeans_JJI$Best.nc)
nb_kmeans_JJI$All.index


DPCM_red2 <- DPCM_red %>%
  filter(LAU_Code %notin% c("015146"))

DPCM_red2$kmpp_nogeo <- clustering[[3]]$cluster
DPCM_red2$kmpp_geo <- clustering_sp[[4]]$cluster

DPCM_red2 %>%
  mutate(var = as.factor(kmpp_nogeo)) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = var)) + 
  scale_fill_viridis_d(option = "plasma")
DPCM_red2 %>%
  mutate(var = as.factor(kmpp_geo)) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = var)) + 
  scale_fill_viridis_d(option = "plasma")

DPCM_red %>%
  group_by(kmpp_geo) %>%
  summarise(across(is.numeric, function(x) mean(x,na.rm=T))) %>%
  ungroup() %>%
  View(title = "Geo")
DPCM_red %>%
  group_by(kmpp_nogeo) %>%
  summarise(across(is.numeric, function(x) mean(x,na.rm=T))) %>%
  ungroup() %>%
  View(title = "No geo")

  

##### Follow the guide: https://geographicdata.science/book/notebooks/10_clustering_and_regionalization.html  
# 1. Clustering without geography.
#     Attention: the scale matters! Clustering is based on the concept of distances, and distances
#     (such as variance) are scale dependent/proportional (i.e., variables with large values will
#     have large variabilities) --> Standardize
# 2. Global spatial autocorrelation measures (e.g. Moran's I) to detect which variables
#     contain spatial structure: if positive spatial autocorrelation is estimated,
#     a clear spatial structure is embedded in the data.
#     The fact that all of the clustering variables are positively autocorrelated does not 
#     say much about how attributes co-vary over space. To explore cross-attribute relationships,
#     we need to consider the spatial correlation between variables.

  




##### Spatially adjacency augmented k-means clustering method (SAAKM)
# https://link.springer.com/article/10.1007/s11042-019-08009-x
# https://igraph.org/r/html/1.3.4/laplacian_matrix.html
# https://rdrr.io/github/abylayzhumekenov/lanczos/man/lanczos.html
# https://en.wikipedia.org/wiki/Adjacency_matrix
# https://en.wikipedia.org/wiki/Laplacian_matrix
# https://en.wikipedia.org/wiki/Graph_theory

dh <- sf::st_distance(DPCM_red2,which = "Great Circle")
dh_mat <- as.matrix(dh)
attr(dh_mat,"class") <- NULL
range_KW <- max(dh_mat)
adj_mat <- exp(-dh_mat/(2*range_KW))

library(factoextra)
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
pca_dhmat <- prcomp(adj_mat,center = TRUE, scale. = TRUE)
fviz_eig(pca_dhmat,ncp = 10)
# dh_lanc <- mgcv::slanczos(A = adj_mat,k = 1)
# plot(pca_dhmat$x[,1],dh_lanc$vectors[,1])

# https://igraph.org/r/html/1.3.1/graph_from_adjacency_matrix.html
# https://search.r-project.org/CRAN/refmans/diffudist/html/getSpectralDecomp.html
library(igraph)
g <- graph_from_adjacency_matrix(adjmatrix = adj_mat, weighted=TRUE, mode="undirected")
library(diffudist)
lapl <- get_spectral_decomp(g, type = "Normalized Laplacian", verbose = TRUE)
lapl <- get_laplacian(g, type = "Normalized Laplacian")
# lapl <- laplacian_matrix(graph = g, normalized = T)
pca_lapl <- prcomp(lapl,center = F, scale. = F)
fviz_eig(pca_lapl,ncp = 10)


# Xmat_full <- cbind(pca_lapl$x[,1:10],as.matrix(CData[,-c(1:6)]))
spcomp <- apply(pca_dhmat$x[,1:3],MARGIN = 2,FUN = function(x) scale(x,center = T,scale = T))
Xmat_full <- cbind(spcomp,as.matrix(CData[,-c(1:6)]))
DMat_full <- get_dist(Xmat_full)

########## Spatially-informed K-means++ 
clustering_sp <- list()
clust_sp <- seq(from=2,to=30)
indexes_sp <- matrix(NA,length(clust),3)
set.seed(12345)

for (cl in 1:length(clust_sp)){
  cat(paste0("Cluster ",cl," of ",length(clust_sp),"\n"))
  clustering_sp[[cl]] <- kmpp(X = Xmat_full, k = clust_sp[cl])
  stats_sp <- cluster.stats(DMat_full, clustering_sp[[cl]]$cl)
  indexes_sp[cl,] <- c(stats_sp$avg.silwidth, stats_sp$pearsongamma,stats_sp$within.cluster.ss)
}
# Performance indices
colnames(indexes_sp) <- c("ASW", "PG","WSS")
rownames(indexes_sp) <- paste0("cl = ",clust)
indexes
# Plot
win.graph()
par(mfrow = c(3,1))
plot(clust_sp,indexes_sp[,1], type = "b")
plot(clust_sp,indexes_sp[,2], type = "b")
plot(clust_sp,indexes_sp[,3], type = "b")
par(mfrow = c(3,1))






str(adj_mat)
DMat_m <- as.matrix(DMat)
attr(DMat_m,"class") <- NULL
range_DMat_m <- max(DMat_m)
adj_DMat <- exp(-DMat_m/(2*range_DMat_m))

alpha <- 0.1

WMat <- (1-alpha)*adj_DMat + alpha*adj_mat


res_kmpp_nogeo <- NbClust_mod(diss = WMat,
                              distance = NULL,
                              min.nc = 2, max.nc = 4,
                              method = "kmpp", 
                              index ="all")




#################################################################
library(ClustGeo)

dh <- sf::st_distance(DPCM_red,which = "Great Circle")
# dh_mat <- as.matrix(dh)
# attr(dh_mat,"class") <- NULL
# range_KW <- max(dh_mat)
# adj_mat <- exp(-dh_mat/(2*range_KW))

D0 <- DMat
# D0 <- D0/max(D0)
D1 <- DMat_geo
# D1 <- D1/max(D1)
n <- nrow(DMat)

tree_D0 <- hclustgeo(D0 = D0, wt=rep(1/n,n))
plot(tree_D0,hang=-1, label=FALSE, xlab="", sub="", main="")
rect.hclust(tree_D0, k=10, border=c(4, 5, 3, 2, 1))
legend("topright", legend=paste("cluster", 1:8), fill=1:8, bty="n", border="white")

P5 <- cutree(tree_D0, 5) # cut the dendrogram to get the partition in 5-8 clusters
P6 <- cutree(tree_D0, 6) # cut the dendrogram to get the partition in 5-8 clusters
P7 <- cutree(tree_D0, 7) # cut the dendrogram to get the partition in 5-8 clusters
P8 <- cutree(tree_D0, 8) # cut the dendrogram to get the partition in 5-8 clusters
sp::plot(estuary$map, border="grey", col=P5) # plot an object of class sp
legend("topleft", legend=paste("cluster", 1:5), fill=1:5, bty="n", border="white")

inertdiss(D, wt=rep(1/n, n)) # the pseudo-inertia of the data
sum(tree$height)



DPCM_red2 <- DPCM_red %>%
  filter(LAU_Code %notin% c("015146"))
# DPCM_red2 <- DPCM_red %>%
#   filter(!is.na(Local_Units_Industry_Active))

DPCM_red2$D0_P5 <- cutree(tree_D0, 5)
DPCM_red2$D0_P6 <- cutree(tree_D0, 6)
DPCM_red2$D0_P7 <- cutree(tree_D0, 7)
DPCM_red2$D0_P8 <- cutree(tree_D0, 8)

cr1 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=1, graph=TRUE) # alpha = 0.10
cr2 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=2, graph=TRUE) # alpha = 0.10
cr3 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=3, graph=TRUE) # alpha = 0.10
cr4 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=4, graph=TRUE) # alpha = 0.10
cr5 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=5, graph=TRUE) # alpha = 0.10
cr6 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=6, graph=TRUE) # alpha = 0.10
cr7 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=7, graph=TRUE) # alpha = 0.10
cr8 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=8, graph=TRUE) # alpha = 0.10
cr9 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=9, graph=TRUE) # alpha = 0.10
cr10 <- choicealpha(D0, D1, range.alpha=seq(0, 1, 0.1), K=10, graph=TRUE) # alpha = 0.10


##### Distance matrix among units: socio-economic variables
DMat <- get_dist(Xmat_nogeo)
f1 <- fviz_dist(DMat, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
DPCM_red2$AvgDist_D0 <- rowMeans(as.matrix(DMat))

##### Distance matrix among units: geographical variables
DMat_geo <- get_dist(Xmat_geo)
f2 <- fviz_dist(DMat_geo, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
DPCM_red2$AvgDist_D1 <- rowMeans(as.matrix(DMat_geo))

##### Mapping average distances
png("AvgDist_D0.png", width = 1800, height = 1200, res = 150)
DPCM_red2 %>%
  select(LAU_Code,AvgDist_D0) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = AvgDist_D0)) + 
  scale_fill_viridis_c("Eucl. Dist.",option = "plasma") + 
  labs(title = "Average Euclidean distance on DPCM-related variables",
       subtitle = "Milan is excluded due to its uniqueness (out-of-scale for every considered variable)") + 
  theme_bw() + 
  theme(plot.title = element_text(size = 14, face = "bold"))
dev.off()


cr1_Q0 <- cbind(cr1$Q[,1],cr1$Qnorm[,1])
cr2_Q0 <- cbind(cr2$Q[,1],cr2$Qnorm[,1])
cr3_Q0 <- cbind(cr3$Q[,1],cr3$Qnorm[,1])
cr4_Q0 <- cbind(cr4$Q[,1],cr4$Qnorm[,1])
cr5_Q0 <- cbind(cr5$Q[,1],cr5$Qnorm[,1])
cr6_Q0 <- cbind(cr6$Q[,1],cr6$Qnorm[,1])
cr7_Q0 <- cbind(cr7$Q[,1],cr7$Qnorm[,1])
cr8_Q0 <- cbind(cr8$Q[,1],cr8$Qnorm[,1])
cr9_Q0 <- cbind(cr9$Q[,1],cr9$Qnorm[,1])
cr10_Q0 <- cbind(cr10$Q[,1],cr10$Qnorm[,1])
cr_out <- rbind(cr1_Q0[1,],cr2_Q0[1,],cr3_Q0[1,],cr4_Q0[1,],cr5_Q0[1,],
                cr6_Q0[1,],cr7_Q0[1,],cr8_Q0[1,],cr9_Q0[1,],cr10_Q0[1,])[1,]
diff(cr_out)










Cl_D0 <- NbClust(data = Xmat_geo, method = "ward.D", distance = "euclidean")
plot(Cl_D0$All.index[,"Silhouette"], type = "b")
plot(Cl_D0$All.index[,"Duda"], type = "b")
plot(Cl_D0$All.index[,"KL"], type = "b")
plot(Cl_D0$All.index[,"CH"], type = "b")

treeD0D1 <- hclustgeo(D0, D1, alpha=0.1, wt = rep(1/n, n))

cols <- c("Q0" = "#CC0000",
          "Q1" = "#3399FF",
          "Q0norm" = "#FF9933",
          "Q1norm" = "#0000FF")

kopt <- 8
cr8_out <- as.data.frame(cr8$Q)
cr8_out$Q0norm <- cr8$Qnorm[,1]
cr8_out$Q1norm <- cr8$Qnorm[,2]
cr8_out$alpha <- cr8$range.alpha
png("ChoiceAlpha.png", width = 1800, height = 1200, res = 150)
cr8_out %>%
  pivot_longer(cols = c("Q0","Q1","Q0norm","Q1norm"), names_to = "Stat", values_to = "Value") %>%
  mutate(Value = Value*100,
         Grp = case_when(Stat %in% c("Q0","Q1") ~ "Non-normalized pseudo inertia",
                         Stat %in% c("Q0norm","Q1norm") ~ "Normalized pseudo inertia",
                         TRUE ~ NA_character_),
         Stat = factor(Stat, levels = c("Q0","Q1","Q0norm","Q1norm"),
                       labels = c("Q0","Q1","Q0norm","Q1norm"),ordered = T)) %>%
  ggplot(mapping = aes(x = alpha, y = Value, col = Stat)) + 
  geom_point(size = 2) + 
  geom_line(linewidth = 1) + 
  facet_wrap(~ Grp) + 
  labs(title = latex2exp::TeX("Spatial hierarchical clustering: proportion of explained pseudo inertia using K=8 clusters"),
       y = "% of explained pseudo inertia", x = latex2exp::TeX("$\\alpha$")) + 
  scale_y_continuous(breaks=seq(0,100,5)) + 
  scale_x_continuous(breaks=seq(0,1,0.10)) + 
  scale_color_manual("",values = cols) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"))
dev.off()





png("Dendrogram.png", width = 1800, height = 1200, res = 150)
plot(treeD0D1,hang=-1, label=FALSE, xlab="",ylab = "Aggregation height",
     sub = latex2exp::TeX("Mixing parameter of the geographical distance matrix: $\\alpha$ = 0"),
     main="Ward dendrogram with both DPCM and geographical features")
par(lwd=4, mar=c(0,0,0,0))
rect.hclust(treeD0D1, k=kopt, border=c(1:kopt))
legend("topright", legend=paste("cluster", 1:kopt), fill=1:kopt, bty="n", border="white", pt.cex = 14)
dev.off()




DPCM_red2$D0D1_P5 <- cutree(treeD0D1, 5)
DPCM_red2$D0D1_P6 <- cutree(treeD0D1, 6)
DPCM_red2$D0D1_P7 <- cutree(treeD0D1, 7)
DPCM_red2$D0D1_P8 <- cutree(treeD0D1, 8)
DPCM_red2$D0D1_P9 <- cutree(treeD0D1, 11)

DPCM_red2 %>%
  select(D0D1_P5,D0D1_P6,D0D1_P7,D0D1_P8,D0D1_P9) %>%
  pivot_longer(cols = c(D0D1_P5,D0D1_P6,D0D1_P7,D0D1_P8,D0D1_P9),
               names_to = "D0D1", values_to = "Clusters") %>%
  mutate(Clusters = as.factor(Clusters)) %>%
  # filter(D0D1 %in% c("D0D1_P8")) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Clusters)) + 
  scale_fill_viridis_d(option = "plasma") + 
  facet_wrap(~ D0D1)


png("NormVariablesD0.png", width = 1800, height = 1200, res = 150)
DPCM_red2 %>%
  select(Local_Units_Industry_Suspended,Value_Added_Industry_Suspended,Employed_Industry_Suspended,
         Local_Units_Services_Suspended,Value_Added_Services_Suspended,Employed_Services_Suspended) %>%
  mutate(across(is.numeric, function(x) x/max(x,na.rm = T))) %>%
  pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended,
                        Value_Added_Industry_Suspended,Value_Added_Services_Suspended,
                        Employed_Industry_Suspended,Employed_Services_Suspended),
               names_to = "Vars", values_to = "Values") %>%
  mutate(Vars = case_when(
    Vars %in% c("Local_Units_Industry_Suspended") ~ "Industry: local units suspended",
    Vars %in% c("Value_Added_Industry_Suspended") ~ "Industry: value adde of local units suspended",
    Vars %in% c("Employed_Industry_Suspended") ~ "Industry: employed persons suspended",
    Vars %in% c("Local_Units_Services_Suspended") ~ "Services: local units suspended",
    Vars %in% c("Value_Added_Services_Suspended") ~ "Services: value adde of local units suspended",
    Vars %in% c("Employed_Services_Suspended") ~ "Services: employed persons suspended",
    TRUE ~ NA_character_
  ),
  Vars = factor(
    Vars,
    levels = c("Industry: local units suspended","Industry: value adde of local units suspended",
               "Industry: employed persons suspended","Services: local units suspended",
               "Services: value adde of local units suspended","Services: employed persons suspended"),
    labels = c("Industry: local units suspended","Industry: value adde of local units suspended",
               "Industry: employed persons suspended","Services: local units suspended",
               "Services: value adde of local units suspended","Services: employed persons suspended"),
    ordered = T)) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Values)) + 
  scale_fill_viridis_c(option = "plasma") + 
  facet_wrap(~ Vars, nrow = 2, ncol = 3) + 
  labs(title = "Normalized DPCM-related variables by industry and services",
       subtitle = "Milan is excluded due to its uniqueness (out-of-scale for every considered variable)",
       y = "", x = "") + 
  theme_bw() + 
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))
dev.off()



##### Add Milano
Milano <- DPCM_red %>% 
  filter(LAU_Code %in% c("015146")) %>%
  mutate(D0D1_P5 = 9,D0D1_P6 = 9,D0D1_P7 = 9,D0D1_P8 = 9,
         D0_P5 = 9, D0_P6 = 9, D0_P7 = 9, D0_P8 = 9)
Data_clusters <- bind_rows(DPCM_red2,Milano)

##### Add regional mean
Regione <- Data_clusters %>%
  summarise(n = sum(n),
            across(is.numeric, function(x) median(x,na.rm=T))) %>%
  mutate(D0D1_P5 = 10,D0D1_P6 = 10,D0D1_P7 = 10,D0D1_P8 = 10,
         D0_P5 = 10, D0_P6 = 10, D0_P7 = 10, D0_P8 = 10)
Data_clusters <- bind_rows(Data_clusters,Regione)

##### Labelling clusters
Data_clusters <- Data_clusters %>%
  mutate(Cluster8_Geo = case_when(D0D1_P8 == 1 ~ "Rural Pavia and South-West",
                                  D0D1_P8 == 2 ~ "Southern agro Milano, Lodi and Cremona",
                                  D0D1_P8 == 3 ~ "Varese, Lecco e Como",
                                  D0D1_P8 == 4 ~ "Large cities and urban areas",
                                  D0D1_P8 == 5 ~ "Cities of Brescia and San Donato",
                                  D0D1_P8 == 6 ~ "Rural Bergamo and Brescia",
                                  D0D1_P8 == 7 ~ "Alpine Sondrio, Brescia and Bergamo",
                                  D0D1_P8 == 8 ~ "Southern Brescia and Mantova",
                                  D0D1_P8 == 9 ~ "Milano",
                                  D0D1_P8 == 10 ~ "Regional median"),
         Cluster8_NoGeo = case_when(D0_P8 == 1 ~ "Low density rural areas",
                                    D0_P8 == 2 ~ "High density rural areas",
                                    D0_P8 == 3 ~ "Milano suburban areas",
                                    D0_P8 == 4 ~ "City of San Donato Milanese",
                                    D0_P8 == 5 ~ "Cities of Segrate and Assago",
                                    D0_P8 == 6 ~ "Sparse urban areas and cities",
                                    D0_P8 == 7 ~ "Cities of Bergamo, Monza and Como",
                                    D0_P8 == 8 ~ "City of Brescia",
                                    D0_P8 == 9 ~ "Milano",
                                    D0_P8 == 10 ~ "Regional median"))


Data_clusters %>%
  dplyr::select(Pop2021,PopDens_km2,DegUrba,Cluster8_Geo,D0D1_P8,contains("Local_Units"),
                contains("Workers"),contains("Value_Added")) %>%
  mutate(Sh_LU_Ind = Local_Units_Industry_Suspended / Local_Units_Industry_Total*100,
         Sh_LU_Serv = Local_Units_Services_Suspended / Local_Units_Services_Total*100,
         Sh_W_Ind = Workers_Industry_Suspended / Workers_Industry_Total*100,
         Sh_W_Serv = Workers_Services_Suspended / Workers_Services_Total*100,
         Sh_VA_Ind = Value_Added_Industry_Suspended / Value_Added_Industry_Total*100,
         Sh_VA_Serv = Value_Added_Services_Suspended / Value_Added_Services_Total*100) %>%
  group_by(Cluster8_Geo,D0D1_P8) %>%
  summarise(n = n(),
            across(is.numeric, function(x) mean(x,na.rm=T))) %>%
  ungroup() %>%
  select(Cluster8_Geo,D0D1_P8,n,Pop2021,PopDens_km2,DegUrba,contains("Sh_")) %>%
  View(title = "Geo")
Data_clusters %>%
  dplyr::select(Pop2021,PopDens_km2,DegUrba,Cluster8_NoGeo,D0_P8,contains("Local_Units"),
                contains("Workers"),contains("Value_Added")) %>%
  mutate(Sh_LU_Ind = Local_Units_Industry_Suspended / Local_Units_Industry_Total*100,
         Sh_LU_Serv = Local_Units_Services_Suspended / Local_Units_Services_Total*100,
         Sh_W_Ind = Workers_Industry_Suspended / Workers_Industry_Total*100,
         Sh_W_Serv = Workers_Services_Suspended / Workers_Services_Total*100,
         Sh_VA_Ind = Value_Added_Industry_Suspended / Value_Added_Industry_Total*100,
         Sh_VA_Serv = Value_Added_Services_Suspended / Value_Added_Services_Total*100) %>%
  group_by(Cluster8_NoGeo) %>%
  summarise(n = n(),
            across(is.numeric, function(x) mean(x,na.rm=T))) %>%
  ungroup() %>%
  select(Cluster8_NoGeo,D0_P8,n,Pop2021,PopDens_km2,DegUrba,contains("Sh_")) %>%
  View(title = "No geo")


##### Colors
colp <- c(
  "Alpine Sondrio, Brescia and Bergamo" = "#996633",
  "Rural Pavia and South-West" = "#FF9900",
  "Varese, Lecco e Como" = "#CC0000",
  "Cities of Brescia and San Donato" = "#0000FF",
  "Large cities and urban areas" = "#3399FF",
  "Rural Bergamo and Brescia"  = "#FFCC33",
  "Milano" = "#333333",
  "Southern agro Milano, Lodi and Cremona" = "#66CC33",
  "Southern Brescia and Mantova" = "#339900",
  
  "Low density rural areas" = "#996633",
  "City of San Donato Milanese" = "#FF9900",
  "Cities of Bergamo, Monza and Como" = "#CC0000",
  "City of Brescia" = "#0000FF",
  "Sparse urban areas and cities" = "#3399FF",
  "Cities of Segrate and Assago"  = "#FFCC33",
  "Milano" = "#333333",
  "Milano suburban areas" = "#66CC33",
  "High density rural areas" = "#339900"
)

p1 <- Data_clusters %>%
  select(Cluster8_Geo,D0D1_P8) %>%
  mutate(Cluster8_Geo = as.factor(Cluster8_Geo)) %>%
  filter(Cluster8_Geo %notin% c("Regional median")) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Cluster8_Geo)) + 
  labs(x = "",
       y = "",
       title = "Geographically-informed Ward hierarchical clustering",
       subtitle = latex2exp::TeX("K = 8 clusters + Milano (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0.10)")) +
  scale_fill_manual("", values = colp) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 16),
        plot.subtitle = element_text(face = "bold",size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) + 
  guides(col=guide_legend(nrow=2,byrow=TRUE),
         linetype=guide_legend(nrow=2,byrow=TRUE))
ggpubr::ggexport(p1,width = 1800, height = 1200, res = 150, filename = "GeoInform_k8.png")


p2 <- Data_clusters %>%
  select(Cluster8_NoGeo,D0_P8) %>%
  mutate(Cluster8_NoGeo = as.factor(Cluster8_NoGeo)) %>%
  filter(Cluster8_NoGeo %notin% c("Regional median")) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Cluster8_NoGeo)) + 
  labs(x = "",
       y = "",
       title = "Ward-like hierarchical clustering",
       subtitle = latex2exp::TeX("K = 8 clusters + Milano (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0)")) +
  scale_fill_manual("", values = colp) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold",size = 16),
        plot.subtitle = element_text(face = "bold",size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10)) + 
  guides(col=guide_legend(nrow=2,byrow=TRUE),
         linetype=guide_legend(nrow=2,byrow=TRUE))
ggpubr::ggexport(p2,width = 1800, height = 1200, res = 150, filename = "NonGeoInform_k8.png")



# save.image(file = "res1911.RData")
