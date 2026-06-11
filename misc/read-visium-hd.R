library(Seurat)
library(ggplot2)
library(dplyr)

###############################################################################
## Fixes for Seurat to work with VisiumHD
###############################################################################

Read10X_Image <- function(image.dir, filter.matrix = TRUE, ...) {
  image <- png::readPNG(source = file.path(image.dir, 'tissue_lowres_image.png'))
  scale.factors <- jsonlite::fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions.path <- Sys.glob(paths = file.path(image.dir, 'tissue_positions*'))

  # The format changed from csv -> parquet
  # They used `zstd` for compression, which is not installed in the default install of the R package "arrow"
  # This can normally be solved by doing a "full" install, but I didn't have any luck:
  # Sys.setenv("LIBARROW_MINIMAL" = FALSE); install.packages("arrow")
  # Current solution: use conda environment with r-Seurat and r-arrow installed (the conda install of arrow includes all compression algorithms already)
  tissue.positions <- arrow::read_parquet(tissue.positions.path)

  colnames(tissue.positions) <- c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol')
  tissue.positions <- tissue.positions %>% tibble::column_to_rownames("barcodes")

  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))

  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$spot_diameter_fullres,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}

# This is very simplified version of the default Seurat function, but enough to get the data loaded
Load10X_Spatial <- function(data.dir, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = T, to.upper = F, image = NULL, ...) {
  data <- Read10X_h5(filename = file.path(data.dir, filename))

  object <- CreateSeuratObject(counts = data, assay = assay)
  image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- "Spatial"
  object[["slice1"]] <- image

  return(object)
}

##

# NOTE: use the outputs from the binned_outputs directory!
#       I am not sure what the matrices in the root of the spaceranger output are, but there are missing files in the spatial folder which makes the data not usable.
# NOTE: Only got binsize of 48um working. Smaller bin size all result in SCTransform crashes due to NA values in UMI counts.
#       This can probably be resolved by adding a filtering step to remove empty bins.
# 
# data.dir <- "~/topsecret/spaceranger/binned_outputs/square_048um/"
# 
# seuratObj <- Load10X_Spatial(data.dir)
# 
# # The default Seurat point size for spatial plots does not work well with Visium HD. The pt.size scale factor should be set to make the plots somewhat readable.
# # This scale factor can probably be deduced from something, but didn't look into this yet so now its just manual trial and error :)
# 
# VlnPlot(seuratObj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
# SpatialFeaturePlot(seuratObj, features = "nCount_Spatial", pt.size.factor = 0.07) + theme(legend.position = "right")
# 
# seuratObj <- SCTransform(seuratObj, assay = "Spatial", verbose = FALSE)
# seuratObj <- RunPCA(seuratObj, assay = "SCT", verbose = FALSE)
# seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = 1:30)
# seuratObj <- FindClusters(seuratObj, verbose = FALSE)
# seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:30)
# 
# SpatialDimPlot(seuratObj, pt.size.factor = 0.25)
