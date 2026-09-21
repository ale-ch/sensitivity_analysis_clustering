##############################################
########## Spatial cluster analysis ##########
##############################################

setwd("H:/Il mio Drive/FEEM_COVID19_Lombardia/Analisi/SpatialCluster_Fragility")

library(sf)
library(tidyverse)
library(readxl)
library(pracma)
library(cluster)
library(stats)
library(fpc)
library(factoextra) 
library(NbClust)
library(ClustGeo)
library(viridis)
library(ggpubr)

'%notin%' <- Negate('%in%')

AOI <- "Italy"
# AOI <- "Lombardy"

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
DPCM <- read_excel("DPCM/DPCM_cleaned.xlsx")

##### Municipalities registry
ManualMatchLAUs <- read_excel("anagrafica comuni/Manually_matched_LAUs_Oct2023.xlsx", sheet = "Modified")
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

# png("SuspendedDPCM.png", width = 1800, height = 1200, res = 150)
# DPCM_red %>%
#   mutate('Share of suspended local units in industry' = Local_Units_Industry_Suspended/Local_Units_Industry_Total*100,
#          'Share of suspended local units in services' = Local_Units_Services_Suspended/Local_Units_Services_Total*100) %>%
#   select('Share of suspended local units in industry','Share of suspended local units in services') %>%
#   pivot_longer(cols = c('Share of suspended local units in industry',
#                         'Share of suspended local units in services'),
#                names_to = "Var", values_to = "Values") %>%
#   ggplot() + 
#   geom_sf(mapping = aes(fill = Values)) + 
#   scale_fill_viridis_b(option = "plasma") + 
#   facet_wrap(~ Var) + 
#   theme_bw()
# dev.off()
# 
# DPCM_red %>%
#   filter(NUTS2_Name %in% c("Lombardia")) %>%
#   select(Local_Units_Industry_Suspended,Local_Units_Services_Suspended) %>%
#   pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended),
#                names_to = "Var", values_to = "Values") %>%
#   ggplot() + 
#   geom_sf(mapping = aes(fill = Values)) + 
#   scale_fill_viridis_b(option = "plasma",n.breaks = 100) + 
#   facet_wrap(~ Var) + 
#   theme_bw()
# 
# 
# Data_plot <- DPCM_red %>%
#   # filter(NUTS2_Name %in% c("Lombardia")) %>%
#   select(Area_km2,Local_Units_Industry_Suspended,Value_Added_Industry_Suspended,Employed_Industry_Suspended,
#          Local_Units_Services_Suspended,Value_Added_Services_Suspended,Employed_Services_Suspended) %>%
#   mutate(across(is.numeric, function(x) x/Area_km2)) %>%
#   select(-c(Area_km2)) %>%
#   ### Normalization w.r.t. the maximum
#   # mutate(across(is.numeric, function(x) x/max(x,na.rm = T))) %>%
#   ### Min-Max (0,1) normalization
#   # mutate(across(is.numeric, function(x) (x - min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))) %>%
#   ### Robust normalization
#   mutate(across(is.numeric, function(x) (x - median(x,na.rm = T))/(quantile(x,0.75,na.rm = T)-quantile(x,0.25,na.rm = T)))) %>%
#   pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended,
#                         Value_Added_Industry_Suspended,Value_Added_Services_Suspended,
#                         Employed_Industry_Suspended,Employed_Services_Suspended),
#                names_to = "Vars", values_to = "Values") %>%
#   mutate(
#     Vars = case_when(
#       Vars %in% c("Local_Units_Industry_Suspended") ~ "Industry: local units suspended",
#       Vars %in% c("Value_Added_Industry_Suspended") ~ "Industry: value adde of local units suspended",
#       Vars %in% c("Employed_Industry_Suspended") ~ "Industry: employed persons suspended",
#       Vars %in% c("Local_Units_Services_Suspended") ~ "Services: local units suspended",
#       Vars %in% c("Value_Added_Services_Suspended") ~ "Services: value adde of local units suspended",
#       Vars %in% c("Employed_Services_Suspended") ~ "Services: employed persons suspended",
#       TRUE ~ NA_character_
#     ),
#     Vars = factor(
#       Vars,
#       levels = c("Industry: local units suspended","Industry: value adde of local units suspended",
#                  "Industry: employed persons suspended","Services: local units suspended",
#                  "Services: value adde of local units suspended","Services: employed persons suspended"),
#       labels = c("Industry: local units suspended","Industry: value adde of local units suspended",
#                  "Industry: employed persons suspended","Services: local units suspended",
#                  "Services: value adde of local units suspended","Services: employed persons suspended"),
#       ordered = T)
#   )









#############################################
########## Define area of interest ##########
#############################################
if (AOI == "Lombardy") {
  DPCM_red <- DPCM_red %>%
    filter(NUTS2_Name == "Lombardia")
}



###########################################################
########## Spatially-informed Kmeans++ Algorithm ##########
###########################################################
CData <- DPCM_red %>%
  # Filter Lombardy
  filter(NUTS2_Name == "Lombardia") %>%
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

##### Plot variables
Data_plot <- CData %>%
  pivot_longer(cols = c(Local_Units_Industry_Suspended,Local_Units_Services_Suspended,
                        Value_Added_Industry_Suspended,Value_Added_Services_Suspended,
                        Employed_Industry_Suspended,Employed_Services_Suspended),
               names_to = "Vars", values_to = "Values") %>%
  mutate(
    Vars = case_when(
      Vars %in% c("Local_Units_Industry_Suspended") ~ "Industry: local units suspended",
      Vars %in% c("Value_Added_Industry_Suspended") ~ "Industry: value added of local units suspended",
      Vars %in% c("Employed_Industry_Suspended") ~ "Industry: employed persons suspended",
      Vars %in% c("Local_Units_Services_Suspended") ~ "Services: local units suspended",
      Vars %in% c("Value_Added_Services_Suspended") ~ "Services: value added of local units suspended",
      Vars %in% c("Employed_Services_Suspended") ~ "Services: employed persons suspended",
      TRUE ~ NA_character_
    ),
    Vars = factor(
      Vars,
      levels = c("Industry: local units suspended","Industry: value added of local units suspended",
                 "Industry: employed persons suspended","Services: local units suspended",
                 "Services: value added of local units suspended","Services: employed persons suspended"),
      labels = c("Industry: local units suspended","Industry: value added of local units suspended",
                 "Industry: employed persons suspended","Services: local units suspended",
                 "Services: value added of local units suspended","Services: employed persons suspended"),
      ordered = T)
  )
q <- quantile(Data_plot$Values, probs = seq(0.1,1,by=0.03),na.rm = T)
png("DoubleNormVariablesD0_Lombardy.png", width = 1800, height = 1200, res = 150)
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


### Compute Pearson's linear correlation
corr <- CData[,-c(1:6)] %>%
  sf::st_drop_geometry() %>%
  cor() %>%
  reshape2::melt(na.rm = TRUE) %>%
  rename(Reg1 = Var1, Reg2 = Var2, Corr2021 = value)
corrmat <- corr %>%
  pivot_longer(cols = 3:last_col()) %>%
  # Standardize for fill/color scaling in tile
  group_by(name) %>%
  mutate(value_std = scale(x = value, center = T, scale = T)) %>%
  ungroup()

### Plot
p_corr <- corrmat %>%
  ggplot(mapping = aes(x = Reg2, y = Reg1, fill = value_std)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option="magma", alpha = 0.6, name = "Measure", begin = 1, end = 0.60) + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed() + 
  geom_text(aes(x = Reg2, y = Reg1, label = round(value,2)), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) + 
  labs(title = "Linear correlation among DPCM-related variables")
ggexport(p_corr,width = 2400, height = 2000, res = 150, filename = "CorrMat.png")









# ##### Distance
# st_geometry(CData) <- NULL
# Xmat_geo <- as.matrix(CData[,2:3])
# Xmat_nogeo <- as.matrix(CData[,-c(1:6)])
# 
# ##### Distance matrix among units: socio-economic variables
# DMat <- get_dist(Xmat_nogeo,method = "euclidean",stand = FALSE)
# f1 <- fviz_dist(DMat, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
# ggexport(f1,width = 2400, height = 2000, res = 150, filename = "D0_mat.png")
# DMat_geo <- get_dist(Xmat_geo)
# f2 <- fviz_dist(DMat_geo, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
# ggexport(f2,width = 2400, height = 2000, res = 150, filename = "D1_mat.png")
# 
# # AvgDist_D0 <- rowMeans(as.matrix(DMat))
# # DMat_geo_k2 <- mgcv::slanczos(A = as.matrix(DMat_geo),k = 2)
# # Xmat_full <- cbind(DMat_geo_k2$vectors,Xmat_nogeo)
# # DMat_full <- get_dist(Xmat_full)
# 
# 
# ########## Not spatially-informed K-means++ 
# clustering <- list()
# clust <- seq(from=2,to=2)
# indexes <- matrix(NA,length(clust),3)
# set.seed(12345)
# 
# for (cl in 1:length(clust)){
#   cat(paste0("Cluster ",cl," of ",length(clust),"\n"))
#   clustering[[cl]] <- kmpp(X = Xmat_nogeo, k = clust[cl])
#   stats <- cluster.stats(DMat, clustering[[cl]]$cl)
#   indexes[cl,] <- c(stats$avg.silwidth, stats$pearsongamma,stats$within.cluster.ss)
# }
# # Performance indices
# colnames(indexes) <- c("ASW", "PG","WSS")
# rownames(indexes) <- paste0("cl = ",clust)
# indexes
# # Plot
# win.graph()
# par(mfrow = c(3,1))
# plot(clust,indexes[,1], type = "b")
# plot(clust,indexes[,2], type = "b")
# plot(clust,indexes[,3], type = "b")
# par(mfrow = c(3,1))
# 
# 
# 
# 
# ########## Spatially-informed K-means++ 
# clustering_sp <- list()
# clust_sp <- seq(from=2,to=30)
# indexes_sp <- matrix(NA,length(clust),3)
# set.seed(12345)
# 
# for (cl in 1:length(clust_sp)){
#   cat(paste0("Cluster ",cl," of ",length(clust_sp),"\n"))
#   clustering_sp[[cl]] <- kmpp(X = Xmat_full, k = clust_sp[cl])
#   stats_sp <- cluster.stats(DMat_full, clustering_sp[[cl]]$cl)
#   indexes_sp[cl,] <- c(stats_sp$avg.silwidth, stats_sp$pearsongamma,stats_sp$within.cluster.ss)
# }
# # Performance indices
# colnames(indexes_sp) <- c("ASW", "PG","WSS")
# rownames(indexes_sp) <- paste0("cl = ",clust)
# indexes
# # Plot
# win.graph()
# par(mfrow = c(3,1))
# plot(clust_sp,indexes_sp[,1], type = "b")
# plot(clust_sp,indexes_sp[,2], type = "b")
# plot(clust_sp,indexes_sp[,3], type = "b")
# par(mfrow = c(3,1))








#################################################################
############## Hierarchical spatial clustering ##################
#################################################################

##### Distance
st_geometry(CData) <- NULL
Xmat_geo <- as.matrix(CData[,2:3])
Xmat_nogeo <- as.matrix(CData[,-c(1:6)])

# DPCM_red2 <- DPCM_red %>%
#   filter(LAU_Code %notin% c("015146"))
DPCM_red2 <- DPCM_red %>%
  filter(!is.na(Local_Units_Industry_Active))


##### Distance matrix among units: socio-economic variables
DMat <- get_dist(Xmat_nogeo)
f1 <- fviz_dist(DMat, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
ggexport(f1,width = 2400, height = 2000, res = 150, filename = "D0_mat.png")
DPCM_red2$AvgDist_D0 <- rowMeans(as.matrix(DMat))
D0 <- DMat
n <- nrow(DMat)

##### Distance matrix among units: geographical variables
DMat_geo <- get_dist(Xmat_geo)
f2 <- fviz_dist(DMat_geo, gradient = list(low="white", mid="grey", high="black"),show_labels = F)
ggexport(f2,width = 2400, height = 2000, res = 150, filename = "D1_mat.png")
DPCM_red2$AvgDist_D1 <- rowMeans(as.matrix(DMat_geo))
D1 <- DMat_geo

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

##### Standardize distances
D0_std <- D0/max(D0)
D1_std <- D1/max(D1)

##### Ward-like hierarchical clustering without spatial information
tree_D0 <- hclustgeo(D0 = D0, wt=rep(1/n,n))
ind_sil <- NbClust(data=Xmat_nogeo,diss=as.dist(D0),distance=NULL,
                   min.nc = 2, max.nc = 15, method="ward.D2",index="silhouette")
ind_ch <- NbClust(data=Xmat_nogeo,diss=as.dist(D0_std),distance=NULL,
                  min.nc = 2, max.nc = 15, method="ward.D2",index="ch")
ind_dunn <- NbClust(data=Xmat_nogeo,diss=as.dist(D0_std),distance=NULL,
                    min.nc = 2, max.nc = 15, method="ward.D2",index="dunn")
ind_mc <- NbClust(data=Xmat_nogeo,diss=as.dist(D0_std),distance=NULL,
                  min.nc = 2, max.nc = 15, method="ward.D2",index="mcclain")
ind_c <- NbClust(data=Xmat_nogeo,diss=as.dist(D0_std),distance=NULL,
                 min.nc = 2, max.nc = 15, method="ward.D2",index="cindex")

plot(tree_D0,hang=-1, label=FALSE, xlab="", sub="", main="")
rect.hclust(tree_D0, k=8, border=c(4, 5, 3, 2, 1))
legend("topright", legend=paste("cluster", 1:8), fill=1:8, bty="n", border="white")
# 
# P5 <- cutree(tree_D0, 5) # cut the dendrogram to get the partition in 5-8 clusters
# P6 <- cutree(tree_D0, 6) # cut the dendrogram to get the partition in 5-8 clusters
# P7 <- cutree(tree_D0, 7) # cut the dendrogram to get the partition in 5-8 clusters
# P8 <- cutree(tree_D0, 8) # cut the dendrogram to get the partition in 5-8 clusters
# sp::plot(estuary$map, border="grey", col=P5) # plot an object of class sp
# legend("topleft", legend=paste("cluster", 1:5), fill=1:5, bty="n", border="white")
# 
# inertdiss(D, wt=rep(1/n, n)) # the pseudo-inertia of the data
# sum(tree$height)
DPCM_red2$D0_P3 <- cutree(tree_D0, 3)
DPCM_red2$D0_P4 <- cutree(tree_D0, 4)
DPCM_red2$D0_P5 <- cutree(tree_D0, 5)
DPCM_red2$D0_P6 <- cutree(tree_D0, 6)
DPCM_red2$D0_P7 <- cutree(tree_D0, 7)
DPCM_red2$D0_P8 <- cutree(tree_D0, 8)
DPCM_red2$D0_P9 <- cutree(tree_D0, 9)
DPCM_red2$D0_P10 <- cutree(tree_D0, 10)
DPCM_red2$D0_P11 <- cutree(tree_D0, 11)
DPCM_red2$D0_P12 <- cutree(tree_D0, 12)
DPCM_red2$D0_P13 <- cutree(tree_D0, 13)
DPCM_red2$D0_P14 <- cutree(tree_D0, 14)
DPCM_red2$D0_P15 <- cutree(tree_D0, 15)


clustering <- list()
clust <- seq(from=2,to=15)
indexes <- matrix(NA,length(clust),5)
set.seed(12345)
for (cl in 1:length(clust)){
  cat(paste0("Cluster ",cl," of ",length(clust),"\n"))
  clustering[[cl]] <- cutree(tree_D0, clust[cl])
  stats <- cluster.stats(DMat, clustering[[cl]])
  indexes[cl,] <- c(stats$avg.silwidth, stats$pearsongamma,stats$within.cluster.ss,
                    stats$ch,stats$dunn)
}
# Performance indices
colnames(indexes) <- c("ASW", "PG","WSS","CH","Dunn")
rownames(indexes) <- paste0("cl = ",clust)
indexes
# Plot
win.graph()
par(mfrow = c(2,3))
plot(clust,indexes[,1], type = "b") ; which.max(abs(diff(indexes[,1])))
plot(clust,indexes[,2], type = "b") ; which.max(abs(diff(indexes[,2])))
plot(clust,indexes[,3], type = "b") ; which.max(abs(diff(indexes[,3])))
plot(clust,indexes[,4], type = "b") ; which.max(abs(diff(indexes[,4])))
plot(clust,indexes[,5], type = "b") ; which.max(abs(diff(indexes[,5])))
par(mfrow = c(1,1))






#### Choice of the mixing parameter alpha
Estim_alpha <- FALSE
if (Estim_alpha == TRUE) {
  cr1 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 0.1, 0.01), K=1, graph=TRUE) # alpha = 0.10
  cr2 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=2, graph=TRUE) # alpha = 0.10
  cr3 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 0.1, 0.01), K=3, graph=TRUE) # alpha = 0.10
  cr4 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=4, graph=TRUE) # alpha = 0.10
  cr5 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=5, graph=TRUE) # alpha = 0.10
  cr6 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=6, graph=TRUE) # alpha = 0.10
  cr7 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=7, graph=TRUE) # alpha = 0.10
  cr8 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=8, graph=TRUE) # alpha = 0.10
  cr9 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=9, graph=TRUE) # alpha = 0.10
  cr10 <- choicealpha(D0_std, D1_std, range.alpha=seq(0, 1, 0.1), K=10, graph=TRUE) # alpha = 0.10
  save(cr1,cr2,cr3,cr4,cr5,cr6,cr7,cr8,cr9,cr10,D0,D1,D0_std,D1_std,
       file = "C:/Users/paulm/OneDrive/Desktop/Res_Nat_AlphaGivenK.RData")
} else {
  load("H:/Il mio Drive/FEEM_COVID19_Lombardia/Analisi/SpatialCluster_Fragility/Res_Nat_AlphaGivenK.RData")
}


plot(cr2$Qnorm[,2] - cr2$Qnorm[,1])

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
                cr6_Q0[1,],cr7_Q0[1,],cr8_Q0[1,],cr9_Q0[1,],cr10_Q0[1,])

cr_out

diff(cr_out)

Cl_D0 <- NbClust(data = Xmat_geo, method = "ward.D", distance = "euclidean")
plot(Cl_D0$All.index[,"Silhouette"], type = "b")
plot(Cl_D0$All.index[,"Duda"], type = "b")
plot(Cl_D0$All.index[,"KL"], type = "b")
plot(Cl_D0$All.index[,"CH"], type = "b")

##### Ward-like hierarchical clustering with spatial information
# D1_std <- D1/max(D1)
treeD0D1 <- hclustgeo(D0, D1, alpha=0.02, wt = rep(1/n, n))

cols <- c("Q0" = "#CC0000",
          "Q1" = "#3399FF",
          "Q0norm" = "#FF9933",
          "Q1norm" = "#0000FF")

kopt <- 3
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




DPCM_red2$D0D1_P3 <- cutree(treeD0D1, 3)
DPCM_red2$D0D1_P4 <- cutree(treeD0D1, 4)
DPCM_red2$D0D1_P5 <- cutree(treeD0D1, 5)
DPCM_red2$D0D1_P6 <- cutree(treeD0D1, 6)
DPCM_red2$D0D1_P7 <- cutree(treeD0D1, 7)
DPCM_red2$D0D1_P8 <- cutree(treeD0D1, 8)
DPCM_red2$D0D1_P9 <- cutree(treeD0D1, 9)
DPCM_red2$D0D1_P10 <- cutree(treeD0D1, 10)
DPCM_red2$D0D1_P11 <- cutree(treeD0D1, 11)
DPCM_red2$D0D1_P12 <- cutree(treeD0D1, 12)
DPCM_red2$D0D1_P13 <- cutree(treeD0D1, 13)
DPCM_red2$D0D1_P14 <- cutree(treeD0D1, 14)
DPCM_red2$D0D1_P15 <- cutree(treeD0D1, 15)

Pa <- DPCM_red2 %>%
  select(D0D1_P3,D0D1_P5,D0D1_P6,D0D1_P7,D0D1_P8,D0D1_P9,D0D1_P10,D0D1_P11,D0D1_P12,
         D0D1_P13,D0D1_P14,D0D1_P15) %>%
  pivot_longer(cols = c(D0D1_P3,D0D1_P5,D0D1_P6,D0D1_P7,D0D1_P8,D0D1_P9,D0D1_P10,D0D1_P11,D0D1_P12,
                        D0D1_P13,D0D1_P14,D0D1_P15),
               names_to = "D0D1", values_to = "Clusters") %>%
  mutate(Clusters = as.factor(Clusters)) %>%
  filter(D0D1 %in% c("D0D1_P3","D0D1_P8","D0D1_P15")) %>%
  mutate(D0D1 = case_when(D0D1 == "D0D1_P3" ~ "K = 3",
                          D0D1 == "D0D1_P8" ~ "K = 8",
                          D0D1 == "D0D1_P15" ~ "K = 15",
                          TRUE ~ D0D1),
         D0D1 = factor(D0D1,levels = c("K = 3","K = 8","K = 15"),
                       labels = c("K = 3","K = 8","K = 15"),
                       ordered = TRUE)) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Clusters),linewidth = 0) + 
  scale_fill_viridis_d(option = "turbo") + 
  facet_wrap(~ D0D1, ncol = 3) + 
  # theme(legend.position = "none") + 
  theme_bw() + 
  labs(x = "",
       y = "",
       title = "Geographically-informed Ward hierarchical clustering",
       subtitle = latex2exp::TeX("K = 3 or 8 or 15 clusters (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0.02)")
  )
ggpubr::ggexport(Pa,width = 1800, height = 1200, res = 150, filename = "ClusterGeoInform.png")

Pb <- DPCM_red2 %>%
  select(D0_P3,D0_P5,D0_P6,D0_P7,D0_P8,D0_P9,D0_P10,D0_P11,D0_P12,
         D0_P13,D0_P14,D0_P15) %>%
  pivot_longer(cols = c(D0_P3,D0_P5,D0_P6,D0_P7,D0_P8,D0_P9,D0_P10,D0_P11,D0_P12,
                        D0_P13,D0_P14,D0_P15),
               names_to = "D0", values_to = "Clusters") %>%
  mutate(Clusters = as.factor(Clusters)) %>%
  filter(D0 %in% c("D0_P3","D0_P8","D0_P15")) %>%
  mutate(D0 = case_when(D0 == "D0_P3" ~ "K = 3",
                        D0 == "D0_P8" ~ "K = 8",
                        D0 == "D0_P15" ~ "K = 15",
                        TRUE ~ D0),
         D0 = factor(D0,levels = c("K = 3","K = 8","K = 15"),
                     labels = c("K = 3","K = 8","K = 15"),
                     ordered = TRUE)) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Clusters),linewidth = 0) + 
  scale_fill_viridis_d(option = "turbo") + 
  facet_wrap(~ D0, ncol = 3) + 
  # theme(legend.position = "none") + 
  theme_bw() + 
  labs(x = "",
       y = "",
       title = "Ward-like hierarchical clustering",
       subtitle = latex2exp::TeX("K = 3 or 8 or 15 clusters (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0.02)")
  )
ggpubr::ggexport(Pb,width = 2400, height = 2000, res = 225, filename = "ClusterNonGeoInform.png")


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
# Milano <- DPCM_red %>% 
#   filter(LAU_Code %in% c("015146")) %>%
#   mutate(D0D1_P5 = 9,D0D1_P6 = 9,D0D1_P7 = 9,D0D1_P8 = 9,
#          D0_P5 = 9, D0_P6 = 9, D0_P7 = 9, D0_P8 = 9)
# Data_clusters <- bind_rows(DPCM_red2,Milano)

##### 
Data_clusters <- DPCM_red2

##### Add regional/national median
Regione <- Data_clusters %>%
  summarise(n = sum(n),
            across(is.numeric, function(x) median(x,na.rm=T))) %>%
  mutate(D0D1_P5 = 10,D0D1_P6 = 10,D0D1_P7 = 10,D0D1_P8 = 10,
         D0_P5 = 10, D0_P6 = 10, D0_P7 = 10, D0_P8 = 10)
Data_clusters <- bind_rows(Data_clusters,Regione)

##### Labelling clusters
# Data_clusters <- Data_clusters %>%
#   mutate(
#     Cluster8_Geo = case_when(D0D1_P8 == 1 ~ "Rural Pavia and South-West",
#                              D0D1_P8 == 2 ~ "Southern agro Milano, Lodi and Cremona",
#                              D0D1_P8 == 3 ~ "Varese, Lecco e Como",
#                              D0D1_P8 == 4 ~ "Large cities and urban areas",
#                              D0D1_P8 == 5 ~ "Cities of Brescia and San Donato",
#                              D0D1_P8 == 6 ~ "Rural Bergamo and Brescia",
#                              D0D1_P8 == 7 ~ "Alpine Sondrio, Brescia and Bergamo",
#                              D0D1_P8 == 8 ~ "Southern Brescia and Mantova",
#                              D0D1_P8 == 9 ~ "Milano",
#                              D0D1_P8 == 10 ~ "Regional median"),
#     Cluster8_NoGeo = case_when(D0_P8 == 1 ~ "Low density rural areas",
#                                D0_P8 == 2 ~ "High density rural areas",
#                                D0_P8 == 3 ~ "Milano suburban areas",
#                                D0_P8 == 4 ~ "City of San Donato Milanese",
#                                D0_P8 == 5 ~ "Cities of Segrate and Assago",
#                                D0_P8 == 6 ~ "Sparse urban areas and cities",
#                                D0_P8 == 7 ~ "Cities of Bergamo, Monza and Como",
#                                D0_P8 == 8 ~ "City of Brescia",
#                                D0_P8 == 9 ~ "Milano",
#                                D0_P8 == 10 ~ "Regional median")
#   )
Data_clusters <- Data_clusters %>%
  mutate(
    Cluster8_Geo = case_when(D0D1_P8 == 1 ~ "South and Sicily",
                             D0D1_P8 == 2 ~ "Piemonte and western Liguria",
                             D0D1_P8 == 3 ~ "Center-South",
                             D0D1_P8 == 4 ~ "Triveneto",
                             D0D1_P8 == 5 ~ "Lombardy and Emilia",
                             D0D1_P8 == 6 ~ "Center",
                             D0D1_P8 == 7 ~ "Large urban centers",
                             D0D1_P8 == 8 ~ "Sardinia",
                             D0D1_P8 == 10 ~ "National median"),
        Cluster8_NoGeo = case_when(D0_P8 == 1 ~ "Low density municipalities",
                                   D0_P8 == 2 ~ "High density urban centers (sparse)",
                                   D0_P8 == 3 ~ "Low density urban centers (sparse)",
                                   D0_P8 == 4 ~ "High density municipalities (sparse)",
                                   D0_P8 == 5 ~ "Milano hinterland and Torino and Napoli",
                                   D0_P8 == 6 ~ "City of San Donato Milanese",
                                   D0_P8 == 7 ~ "Large urban centers (sparse)",
                                   D0_P8 == 8 ~ "Cities of Milano and Assago",
                                   D0_P8 == 10 ~ "National median")
  )


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


# ##### Colors
# colp <- c(
#   "Alpine Sondrio, Brescia and Bergamo" = "#996633",
#   "Rural Pavia and South-West" = "#FF9900",
#   "Varese, Lecco e Como" = "#CC0000",
#   "Cities of Brescia and San Donato" = "#0000FF",
#   "Large cities and urban areas" = "#3399FF",
#   "Rural Bergamo and Brescia"  = "#FFCC33",
#   "Milano" = "#333333",
#   "Southern agro Milano, Lodi and Cremona" = "#66CC33",
#   "Southern Brescia and Mantova" = "#339900",
#   
#   "Low density rural areas" = "#996633",
#   "City of San Donato Milanese" = "#FF9900",
#   "Cities of Bergamo, Monza and Como" = "#CC0000",
#   "City of Brescia" = "#0000FF",
#   "Sparse urban areas and cities" = "#3399FF",
#   "Cities of Segrate and Assago"  = "#FFCC33",
#   "Milano" = "#333333",
#   "Milano suburban areas" = "#66CC33",
#   "High density rural areas" = "#339900"
# )

##### Colors
colp <- c(
  "South and Sicily" = "#996633",
  "Piemonte and western Liguria" = "#FF9900",
  "Center-South" = "#CC0000",
  "Triveneto" = "#0000FF",
  "Lombardy and Emilia" = "#3399FF",
  "Center"  = "#FFCC33",
  "Large urban centers" = "#333333",
  "Sardinia" = "#66CC33",
  
  "Low density municipalities" = "#996633",
  "High density urban centers (sparse)" = "#FF9900",
  "Low density urban centers (sparse)" = "#CC0000",
  "High density municipalities (sparse)" = "#0000FF",
  "Milano hinterland and Torino and Napoli" = "#3399FF",
  "City of San Donato Milanese"  = "#FFCC33",
  "Large urban centers (sparse)" = "#333333",
  "Cities of Milano and Assago" = "#66CC33"
)

p1 <- Data_clusters %>%
  select(Cluster8_Geo,D0D1_P8) %>%
  mutate(Cluster8_Geo = as.factor(Cluster8_Geo)) %>%
  filter(Cluster8_Geo %notin% c("National median")) %>%
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
ggpubr::ggexport(p1,width = 1800, height = 1200, res = 150, filename = "GeoInform_k8.png")


p2 <- Data_clusters %>%
  select(Cluster8_NoGeo,D0_P8) %>%
  mutate(Cluster8_NoGeo = as.factor(Cluster8_NoGeo)) %>%
  filter(Cluster8_NoGeo %notin% c("National median")) %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = Cluster8_NoGeo)) + 
  labs(x = "",
       y = "",
       title = "Ward-like hierarchical clustering",
       # subtitle = latex2exp::TeX("K = 8 clusters + Milano (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0)")
       subtitle = latex2exp::TeX("K = 8 clusters (Mixing parameter of the geographical distance matrix: $\\alpha$ = 0)")
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
ggpubr::ggexport(p2,width = 1800, height = 1200, res = 150, filename = "NonGeoInform_k8.png")



save.image(file = "Res_national_11082024.RData")
