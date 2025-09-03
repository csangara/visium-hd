library(Seurat)

# Directly load the data (this will give gene symbols)
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Lung_ILA010/", bin.size = 8)
saveRDS(visium_obj, "data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_008um.rds")

# bin.size 16
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Liver_CAW009/", bin.size = 16)
saveRDS(visium_obj, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_016um.rds")

# bin size 32
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Liver_CAW009/", bin.size = 32)
saveRDS(visium_obj, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_032um.rds")

obj <- readRDS("data/scref_Lung_UBla/lung_combined_macs_DCs_seurat.rds")
obj$annot_ident <- Idents(obj)
saveRDS(obj, "data/scref_Lung_UBla/lung_combined_macs_DCs_seurat_annot_ident.rds")
obj$annot_ident

# Rownames of Seurat object not in visium
not_common_genes <- setdiff(rownames(obj), rownames(visium_obj))
common_genes <- intersect(rownames(visium_obj), rownames(obj))

length(not_common_genes)
length(common_genes)

GetAssayData(obj, slot = "counts")[common_genes, ] %>% colSums %>% hist(breaks = 100)
GetAssayData(obj, slot = "counts")[not_common_genes, ] %>% colSums %>% hist(breaks = 100)

# Explore data
# 8um
visium_obj_8um <- readRDS("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_008um.rds")
dim(visium_obj_8um)

sum(visium_obj_8um$nCount_Spatial.008um > 100)

visium_obj_8um_filtered <- subset(visium_obj_8um, subset = nCount_Spatial.008um > 100)
saveRDS(visium_obj_8um_filtered, "data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_008um_filtered100UMI.rds")
rm(visium_obj_8um, visium_obj_8um_filtered); gc()

# 16um
visium_obj_16um <- readRDS("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_016um.rds")
dim(visium_obj_16um)

sum(visium_obj_16um$nCount_Spatial.016um > 100)

visium_obj_16um_filtered <- subset(visium_obj_16um, subset = nCount_Spatial.016um > 100)
saveRDS(visium_obj_16um_filtered, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_016um_filtered100UMI.rds")
rm(visium_obj_16um, visium_obj_16um_filtered); gc()

# 32um
visium_obj_32um <- readRDS("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_032um.rds")
dim(visium_obj_32um)

sum(visium_obj_32um$nCount_Spatial.032um > 100)

visium_obj_32um_filtered <- subset(visium_obj_32um, subset = nCount_Spatial.032um > 100)
saveRDS(visium_obj_32um_filtered, "data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_032um_filtered100UMI.rds")
rm(visium_obj_32um, visium_obj_32um_filtered); gc()