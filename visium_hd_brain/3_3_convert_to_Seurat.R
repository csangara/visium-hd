# To convert the data from AnnData to Seurat, use sceasy
# conda activate python_r
# R
# library(sceasy)
# library(reticulate)
# use_condaenv('python_r')
# sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
#                       outFile='filename.rds')


# Or directly load the data (this will give gene symbols)
visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_MouseBrain/", bin.size = 8)

# To use ENSEMBL ids - construct object yourself
counts_matrix <- Read10X(data.dir = "data/Visium_HD_MouseBrain/binned_outputs/square_008um/filtered_feature_bc_matrix/",
                      gene.column = 1)

seurat_obj <- CreateSeuratObject(counts = counts_matrix, assay="Spatial.008um")

img <- Read10X_Image("data/Visium_HD_MouseBrain/binned_outputs/square_008um/spatial",
                     assay = "Spatial.008um",
                     slice = "slice1.008um")
img <- img[Cells(seurat_obj)]
seurat_obj[["slice1.008um"]] <- img

saveRDS(seurat_obj, "data/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds")
