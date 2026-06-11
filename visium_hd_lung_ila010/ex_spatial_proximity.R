# Installing Giotto
# (Had to update Rcpp and install extra packages on system for terra)
# remotes::install_github('drieslab/GiottoVisuals', ref='R4.1.0')
# remotes::install_github('drieslab/GiottoClass', ref='R4.1.0')
# remotes::install_github('drieslab/GiottoUtils', ref='R4.1.0')
# devtools::install_github('drieslab/Giotto', ref='R4.1.0')
# remotes::install_github('drieslab/GiottoData')

library(Seurat)
library(tidyverse)
library(Giotto)
#### GIOTTO SPATIAL PROXIMITY ANALYSIS ####
## DOESN'T WORK!
## Even when using ROI object, it crashes
data_path <- paste0("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_")
proportions_path <- paste0("visium_hd_lung_ila010/Visium_HD_Lung_ILA010_")
ext <- "_filtered100UMI"
visium_obj <- readRDS(paste0(data_path, "008um", ext, "_pp.rds"))
deconv_props <- read.table(paste0(proportions_path, "008um", ext,
                                  "/proportions_rctd_Visium_HD_Lung_ILA010_008um", ext),
                           header = TRUE)
roi_ext <- "inset"
if (roi_ext == "inset"){
  roi <- c(xmin = 5700, xmax = 8000, ymin = 16800, ymax = 18300)
  width <- 6
  height <- 7
} else {
  roi <- c(xmin = 5200, xmax = 9100, ymin = 12200, ymax = 18800)
  width <- 11
  height <- 8
}

cells_roi <- GetTissueCoordinates(visium_obj) %>%
  filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)

visium_obj_roi <- visium_obj %>% .[, colnames(.) %in% cells_roi]
#square_size <- visium_obj_roi@images[["slice1.008um"]]@scale.factors$spot

deconv_props_roi <- deconv_props %>% `rownames<-`(colnames(visium_obj)) %>% 
  .[cells_roi, ] %>% 
  rownames_to_column("cell_ID")

giotto_obj <- createGiottoObject(raw_exprs = as.matrix(GetAssayData(visium_obj_roi, layer="counts")),
                                 spatial_locs = GetTissueCoordinates(visium_obj_roi))

# giotto_obj@spatial_enrichment$cell$rna$RCTD@enrichDT,
enrObj <- createSpatEnrObj(
  deconv_props_roi,
  name = "DWLS",
  method = "DWLS",
  spat_unit = "cell",
  feat_type = "rna",
  provenance = "cell",
  misc = list()
)

# create spatial enrichment object
giotto_obj <- setGiotto(giotto_obj, enrObj)
giotto_obj <- createSpatialNetwork(giotto_obj, minimum_k = 0)

# saveRDS(giotto_obj, paste0(data_path, "008um", ext, "_pp_giotto.rds"))
# giotto_obj <- readRDS("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_008um_filtered100UMI_pp_giotto.rds")

proximity_table <- cellProximityEnrichmentSpots(giotto_obj,
                                                spatial_network_name = "Delaunay_network",
                                                number_of_simulations = 10)
