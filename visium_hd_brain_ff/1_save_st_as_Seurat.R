library(Seurat)

# Directly load the data (this will give gene symbols)
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_MouseBrain_FF/", bin.size = 8)
saveRDS(visium_obj, "data/Visium_HD_MouseBrain_FF/Visium_HD_MouseBrain_FF_008um.rds")

