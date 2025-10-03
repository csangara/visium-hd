library(Seurat)
library(dplyr)
library(schard)

# Don't forget to run var_names_make_unique() on the h5ad data before this
visium_obj <- schard::h5ad2seurat("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_056um.h5ad")
visium_obj
rownames(visium_obj)
colnames(visium_obj)

visium_obj <- SeuratObject::RenameAssays(visium_obj, assay.name = "RNA", new.assay.name = "Spatial.056um")

# Rename cells
visium_obj <- RenameCells(visium_obj, new.names = visium_obj$cell_ID)

# For resolutions given by SpaceRanger, we can juse use this function

# img <- Read10X_Image("data/Visium_HD_Liver/binned_outputs/square_016um/spatial",
#                      assay = "Spatial.016um",
#                      slice = "slice1.016um")

# But for custom resolutions, we need to modify the code of Read10X_Image itself
image.dir <- "data/Visium_HD_Liver_CAW009/binned_outputs/square_008um/spatial"

# Read in the H&E stain image.
image <- png::readPNG(
  source = file.path(
    image.dir,
    "tissue_lowres_image.png"
  )
)

# Read in the scale factors,
scale.factors <- Read10X_ScaleFactors(
  filename = file.path(image.dir, "scalefactors_json.json")
)

# Multiply scale.factors by 7 (8*7=56)
scale.factors[["spot"]] <- scale.factors[["spot"]]*7

coords <- read.csv("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_056um_centroids.csv", header=TRUE, row.names=1) %>% 
  select(y, x) %>% rename(imagerow=y, imagecol=x)

key <- Key("slice1.056um", quiet = TRUE)

# Create an `sp` compatible `FOV` instance.
fov <- CreateFOV(
  coords[, c("imagerow", "imagecol")],
  type = "centroids",
  radius = scale.factors[["spot"]],
  assay = "Spatial.056um",
  key = key
)

# Build the final `VisiumV2` instance, essentially just adding `image` and
# `scale.factors` to the `fov`.
visium.v2 <- new(
  Class = "VisiumV2",
  boundaries = fov@boundaries,
  molecules = fov@molecules,
  assay = fov@assay,
  key = fov@key,
  image = image,
  scale.factors = scale.factors
)

img <- visium.v2[Cells(visium_obj)]
visium_obj[["slice1.056um"]] <- img

SpatialPlot(visium_obj)

visium_obj # 9515 x 19059

GetTissueCoordinates(visium_obj) %>% head
saveRDS(visium_obj, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_056um.rds")

library(ggplot2)
# Check count distribution
ggplot(data.frame(value = colSums(GetAssayData(visium_obj)))) +
  geom_histogram(aes(x = value), bins = 100) +
  labs(x = "Log10 Count", y = "Frequency") +
  theme_minimal()


