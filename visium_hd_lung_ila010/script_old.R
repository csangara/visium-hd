library(Seurat)
library(tidyverse)

#### Load in data and filter ####
# for (bin_size in c(8, 16, 32)) {
#   bin_size_str <- sprintf("%03dum", bin_size)
#   
#   visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Lung_ILA010/", bin.size = bin_size)
#   bins_to_keep <- visium_obj@meta.data[paste0("nCount_Spatial.", bin_size_str)] > 100
#   visium_obj_filtered <- visium_obj[, bins_to_keep]
#   
#   saveRDS(visium_obj_filtered, paste0("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_",
#                                       bin_size_str, "_filtered100UMI.rds"))
#   
# }


#### Deconvolution is done on the HPC ####
# Results: proportion files

#### Load in deconvolution results ####
bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)
celltypes <- c('Car4_BEC', 'Cap_BEC', 'Venous_BEC', 'Art_BEC', 'Lymph_endo',
               'Mesothelium', 'Cilliated', 'Secretory', 'AT2_cells', 'AT1_cells',
               'Fibroblasts', 'Macrophages', 'Dendriticcells', 'NK_CD8', 'ILC2s',
               'Resident_Tcells', 'Tregs', 'CD4_Tcells', 'Tcells',
               'Plasma_cells', 'GC_Bcells', 'Bcells')

deconv_props <- read.table(paste0("visium_hd_lung_ila010/proportions_rctd_lung_ila010_",
                          bin_size_str, "_filtered100UMI.tsv"),
                   sep = "\t", header = TRUE)

# Calculate cooccurrence matrix based on presence/absence,
# The actual proportion is not considered
cooccurrence <- deconv_props %>% 
  rownames_to_column(var = "bin") %>% 
  pivot_longer(cols = -bin, names_to = "celltype", values_to = "proportion") %>% 
  filter(proportion > 0) %>% 
  select(-proportion) %>% 
  table() %>% crossprod() %>% 
  .[rev(celltypes), rev(celltypes)]
diag(cooccurrence) <- 0

pheatmap::pheatmap(cooccurrence, color = RColorBrewer::brewer.pal(9, "Blues"),
                   cluster_rows = FALSE, cluster_cols = FALSE)