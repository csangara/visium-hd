library(Seurat)

# Directly load the data (this will give gene symbols)
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Liver/", bin.size = 8)
saveRDS(visium_obj, "data/Visium_HD_Liver/Visium_HD_Liver_008um.rds")

visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Liver/", bin.size = 16)
saveRDS(visium_obj, "data/Visium_HD_Liver/Visium_HD_Liver_016um.rds")
