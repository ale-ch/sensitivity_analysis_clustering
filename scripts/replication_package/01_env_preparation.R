# install.packages(c("rmarkdown", "knitr", "sf", "tidyverse", "readxl",
#                   "pracma", "cluster", "fpc", "factoextra", "NbClust",
#                   "ClustGeo", "viridis", "ggpubr", "mclust", "reshape2",
#                   "geosphere", "DT", "htmlwidgets"))


knitr::opts_chunk$set(echo = TRUE, warning = TRUE, message = FALSE,
                      cache = FALSE, fig.width = 12, fig.height = 10)
stopifnot(params$area %in% c("Italia", "Lombardia"),
          is.logical(params$include_milan), length(params$include_milan) == 1L,
          !is.na(params$include_milan),
          params$k_selection_mode %in% c("automatic", "manual"),
          is.numeric(params$k_min), length(params$k_min) == 1L,
          is.finite(params$k_min), params$k_min == as.integer(params$k_min),
          is.numeric(params$k_max), length(params$k_max) == 1L,
          is.finite(params$k_max), params$k_max == as.integer(params$k_max),
          params$k_min >= 2, params$k_max <= 15, params$k_min < params$k_max,
          is.numeric(params$manual_k0), length(params$manual_k0) == 1L,
          is.finite(params$manual_k0), params$manual_k0 == as.integer(params$manual_k0),
          is.numeric(params$manual_final_k), length(params$manual_final_k) == 1L,
          is.finite(params$manual_final_k), params$manual_final_k == as.integer(params$manual_final_k),
          params$alpha_selection_mode %in% c("automatic", "manual"),
          identical(params$alpha_criterion, "paper_absolute"),
          is.numeric(params$alpha_step), length(params$alpha_step) == 1L,
          is.finite(params$alpha_step),
          params$alpha_step > 0, params$alpha_step <= 1,
          is.numeric(params$manual_alpha), length(params$manual_alpha) == 1L,
          is.finite(params$manual_alpha),
          params$manual_alpha >= 0, params$manual_alpha <= 1,
          params$map_palette %in% c("sobria", "pastel"))

base_path <- normalizePath(path.expand(params$base_path), mustWork = TRUE)
# Equivalent to the original setwd(base_path) for subsequent chunks.
knitr::opts_knit$set(root.dir = base_path)
study_area <- if (params$area == "Italia") "Italia" else "Lombardia"
include_milan <- if (study_area == "Italia") TRUE else params$include_milan
clustering_variables <- c(
  "Local_Units_Industry_Suspended", "Local_Units_Services_Suspended",
  "Value_Added_Industry_Suspended", "Value_Added_Services_Suspended",
  "Employed_Industry_Suspended", "Employed_Services_Suspended",
  "Local_Units_Industry_Active", "Local_Units_Services_Active",
  "Value_Added_Industry_Active", "Value_Added_Services_Active",
  "Employed_Industry_Active", "Employed_Services_Active")
stopifnot(length(clustering_variables) == 12L,
          !any(grepl("^(Revenues|Workers)_", clustering_variables)))
K_candidati <- seq.int(as.integer(params$k_min), as.integer(params$k_max))
if (params$k_selection_mode == "manual") {
  stopifnot(params$manual_k0 %in% K_candidati,
            params$manual_final_k %in% K_candidati)
}
alpha_fisso <- 0.10
# Explicitly adding 1 also handles steps such as 0.03 or 0.3.
# Rounding removes floating-point noise from decimal steps.
range_alpha <- sort(unique(c(0,
                             round(seq(0, 1, by = params$alpha_step), digits = 12L), 1)))
stopifnot(range_alpha[1L] == 0, tail(range_alpha, 1L) == 1,
          all(range_alpha >= 0 & range_alpha <= 1))

run_name <- paste(params$area,
                  if (include_milan) "con_Milano" else "senza_Milano",
                  format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "_")
output_dir <- file.path(base_path, "Output_Selezione_K_alpha_RMarkdown", run_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines("IN_CORSO: l'attestato finale viene scritto solo a esecuzione completata.",
           file.path(output_dir, "STATO.txt"))

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
if (!requireNamespace("mclust", quietly = TRUE)) {
  stop("Installare il pacchetto mclust per calcolare gli ARI.")
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
  stop("Installare il pacchetto reshape2, usato dalla correlazione originale.")
}
if (!requireNamespace("geosphere", quietly = TRUE)) {
  stop("Installare geosphere: install.packages('geosphere').")
}

if (!requireNamespace("DT", quietly = TRUE) ||
    !requireNamespace("htmlwidgets", quietly = TRUE)) {
  stop("Installare DT e htmlwidgets: install.packages(c('DT', 'htmlwidgets')).")
}

'%notin%' <- Negate('%in%')
# Version compatibility: use the same empty-vector check as in the loop.
if (!exists("is_empty", mode = "function")) is_empty <- rlang::is_empty
if (file.exists(file.path(base_path, "AuxFuns_NbClust.R"))) {
  source(file.path(base_path, "AuxFuns_NbClust.R"), encoding = "UTF-8")
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_iniziale.txt"))
knitr::kable(data.frame(Area = params$area, Selezione_K = params$k_selection_mode,
                        k_min = min(K_candidati), k_max = max(K_candidati),
                        Milano_nella_stima = include_milan,
                        Selezione_alpha = params$alpha_selection_mode,
                        Criterio_alpha = params$alpha_criterion,
                        Passo_alpha = params$alpha_step, Numero_alpha = length(range_alpha),
                        Palette_mappe = params$map_palette,
                        Riferimento_trasformazioni = "dataset selezionato"))





