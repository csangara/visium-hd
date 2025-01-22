library(Seurat)

# To use ENSEMBL ids - construct object yourself
counts_matrix <- Read10X(data.dir = "data/Visium_HD_MouseBrain_FF/binned_outputs/square_008um/filtered_feature_bc_matrix/",
                         gene.column = 1)

seurat_obj <- CreateSeuratObject(counts = counts_matrix, assay="Spatial.008um")

img <- Read10X_Image("data/Visium_HD_MouseBrain_FF/binned_outputs/square_008um/spatial",
                     assay = "Spatial.008um",
                     slice = "slice1.008um")
img <- img[Cells(seurat_obj)]
seurat_obj[["slice1.008um"]] <- img

saveRDS(seurat_obj, "data/Visium_HD_MouseBrain_FF/Visium_HD_MouseBrain_FF_008um.rds")
