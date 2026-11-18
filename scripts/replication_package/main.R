params <- NULL
params$area <- "Lombardia"
params$study_area <- "Lombardia"
params$k_selection_mode <- "automatic"
params$k_min <- 1L
params$base_path <- '/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering'
params$include_milan <- FALSE
params$k_selection_mode <- "automatic"
params$k_min <- 2
params$k_max <- 15
params$manual_k0 <- 8
params$manual_final_k <- 8
params$alpha_selection_mode <- "automatic"
params$alpha_criterion <- "paper_absolute"
params$alpha_step <- 0.3
params$manual_alpha <- 0.1
params$map_palette <- "sobria"
params$save_distances <- TRUE



base_path <- "/Volumes/T7 Shield/FRES/repos/sensitivity_analysis_clustering"
data_path <- file.path(base_path, "data")

source(file.path(base_path, "scripts/replication_package/01_env_preparation.R"))
source(file.path(base_path, "scripts/replication_package/02_data_preparation.R"))
source(file.path(base_path, "scripts/replication_package/03_params_optim.R"))
source(file.path(base_path, "scripts/replication_package/04_plots_setup.R"))
source(file.path(base_path, "scripts/replication_package/plots/01_plot_hyperparameter_diagnostics.R"))
source(file.path(base_path, "scripts/replication_package/plots/02_plot_spatial_clustering_maps.R"))
source(file.path(base_path, "scripts/replication_package/plots/03_plot_ari_stability_heatmaps.R"))
source(file.path(base_path, "scripts/replication_package/plots/04_plot_exploratory_eda.R"))