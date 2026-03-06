library(Seurat)

# Directly load the data (this will give gene symbols)
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Liver_CAW009/", bin.size = 8)
saveRDS(visium_obj, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_008um.rds")

# bin.size 16
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Liver_CAW009/", bin.size = 16)
saveRDS(visium_obj, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_016um.rds")

# bin size 32
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Liver_CAW009/", bin.size = 32)
saveRDS(visium_obj, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_032um.rds")