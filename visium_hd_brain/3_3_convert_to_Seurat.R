# To convert the data from AnnData to Seurat, use sceasy
# conda activate python_r
# R
# library(sceasy)
# library(reticulate)
# use_condaenv('python_r')
# sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
#                       outFile='filename.rds')