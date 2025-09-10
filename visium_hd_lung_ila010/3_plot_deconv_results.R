library(Seurat)
library(tidyverse)

color_palette <- c("Bcells" = "#498db4",
                   "GC_Bcells" = "#80b654",
                   "Plasma_cells" = "#74b88b",
                   "Tcells" = "#792c68",
                   'CD4_Tcells' = '#d87743',
                   'Tregs' = '#41377d',
                   'Resident_Tcells' = '#94679a',
                   'ILC2s' = '#b7cf64',
                   'NK_CD8' = '#f6de5e',
                   'Dendriticcells' = '#c8779f',
                   'Macrophages' = '#6f4a86',
                   'Fibroblasts' = '#7bbc9d',
                   'AT1_cells' = '#bfa373',
                   'AT2_cells' = '#604988',
                   'Secretory' = '#293877',
                   'Cilliated' = '#475e99',
                   'Mesothelium' = '#7498bf',
                   'Lymph_endo' = '#ae4a75',
                   'Art_BEC' = '#7e432b',
                   'Venous_BEC' = '#6eaf5d',
                   'Cap_BEC' = '#4a9850',
                   'Car4_BEC' = '#544b89'
                   )

data_path <- paste0("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_")
proportions_path <- paste0("visium_hd_lung_ila010/Visium_HD_Lung_ILA010_")
ext <- "_filtered100UMI"
plot_path <- paste0("visium_hd_lung_ila010/plots/")

bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)

for (bin_size in c(8, 16, 32)) {
  print(bin_size)
  bin_size_str <- sprintf("%03dum", bin_size)
  
  visium_obj <- readRDS(paste0(data_path, bin_size_str, ext, ".rds"))
  dim(visium_obj) # 19070 x 143657 spots
  
  deconv_props <- read.table(paste0(proportions_path, bin_size_str, ext,
                                    "/proportions_rctd_Visium_HD_Lung_ILA010_",
                                    bin_size_str, ext),
                             header = TRUE)
  dim(deconv_props) #143657 spots
  
  # Check if removed rows + leftover rows == total rows (yes)
  dim(deconv_props)[1] == dim(visium_obj)[2]
  
  # Add rownames to deconv_props
  rownames(deconv_props) <- colnames(visium_obj)
  all(rownames(deconv_props) == colnames(visium_obj))
  
  # Save deconv_props with rownames
  write.table(deconv_props, file = paste0("visium_hd_lung_ila010/proportions_rctd_lung_ila010_",
                                          bin_size_str, ext, ".tsv"),
              sep = "\t", quote = FALSE)
  
  # Calculate cooccurrence matrix
  cooccurrence <- deconv_props %>% 
    rownames_to_column(var = "bin") %>% 
    pivot_longer(cols = -bin, names_to = "celltype", values_to = "proportion") %>% 
    filter(proportion > 0) %>% 
    select(-proportion) %>% 
    table() %>% 
    crossprod() %>% 
    # Order by names(color_palette)
    .[names(color_palette), names(color_palette)]
  diag(cooccurrence) <- 0
  
  pheatmap::pheatmap(cooccurrence, color = RColorBrewer::brewer.pal(9, "Blues"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     filename = paste0(plot_path, "cooccurrence_heatmap_", bin_size_str, ".png"))
  
  
}


# Assign cell type with max proportion to each spot
# visium_obj$celltype <- factor(colnames(deconv_props)[max.col(deconv_props)],
#                               levels = names(color_palette))
# visium_obj$celltype %>% table %>% sort
# 
# p_celltype <- SpatialDimPlot(visium_obj, group.by = "celltype",
#                              image.alpha = 0, stroke=NA) +
#   scale_fill_manual(values = color_palette) +
#   guides(fill = guide_legend(override.aes = list(size = 3))) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(0.2, "cm"))
# p_celltype
# ggsave(paste0(plot_path, "spatialdimplot_celltype_", bin_size_str, ".png"), p_celltype,
#        width = 10, height = 12, bg = "white")


