library(Seurat)
library(tidyverse)
library(patchwork)
library(sf)
library(ggnewscale)
library(ggplotify)

if (use_color_palette == "_k"){
  color_palette <- c(
    "Bcells" = "#1f1fad",
    "GC_Bcells" = "#0000FF",
    "Plasma_cells" = "#F6C644",
    "Tcells" = "#008000",
    "CD4_Tcells" = "#008000",
    "Tregs" = "#008000",
    "Resident_Tcells" = "#008000",
    "ILC2s" = "#016c59",
    "NK_CD8" = "#016c59",
    "Dendriticcells" = "#c77c7c",
    "Macrophages" = "#c77c7c",
    "Fibroblasts" = "#ff000d",
    "AT1_cells" = "#00FF00",
    "AT2_cells" = "#FF00FF",
    "Secretory" = "#fee3cd",
    "Cilliated" = "#e09252",
    "Mesothelium" = "#008585",
    "Lymph_endo" = "#008585",
    "Art_BEC" = "#abc4ff",
    "Venous_BEC" = "#abc4ff",
    "Cap_BEC" = "#abc4ff",
    "Car4_BEC" = "#abc4ff"
  )
  
} else {
  color_palette <- c(
    "Bcells" = "#41b6c4",
    "GC_Bcells" = "#7fcdbb",
    "Plasma_cells" = "#fb9a99",
    "Tcells" = "#1a9850",
    "CD4_Tcells" = "#1a9850",
    "Tregs" = "#1a9850",
    "Resident_Tcells" = "#1a9850",
    "ILC2s" = "#016c59",
    "NK_CD8" = "#016c59",
    "Dendriticcells" = "#b15929",
    "Macrophages" = "#8c510a",
    "Fibroblasts" = "#e31a1c",
    "AT1_cells" = "#cab2d6",
    "AT2_cells" = "#6a3d9a",
    "Secretory" = "#fec44f",
    "Cilliated" = "#fee391",
    "Mesothelium" = "#253494",
    "Lymph_endo" = "#253494",
    "Art_BEC" = "#4575b4",
    "Venous_BEC" = "#4575b4",
    "Cap_BEC" = "#4575b4",
    "Car4_BEC" = "#4575b4"
  )
}
proper_celltype_names <- 
  c("Bcells" = "B cells",
    "GC_Bcells" = "GC B cells",
    "Plasma_cells" = "ASCs",
    "Tcells" = "T cells",
    "CD4_Tcells" = "CD4 T cells",
    "Tregs" = "Tregs",
    "Resident_Tcells" = "Resident T cells",
    "ILC2s" = "ILC2s",
    "NK_CD8" = "NK/CD8 T cells",
    "Dendriticcells" = "Dendritic cells",
    "Macrophages" = "Macrophages",
    "Fibroblasts" = "Fibroblasts",
    "AT1_cells" = "AT1 cells",
    "AT2_cells" = "AT2 cells",
    "Secretory" = "Secretory",
    "Cilliated" = "Ciliated",
    "Mesothelium" = "Mesothelium",
    "Lymph_endo" = "Lymphatics",
    "Art_BEC" = "Arterial BEC",
    "Venous_BEC" = "Venous BEC",
    "Cap_BEC" = "Capillary BEC",
    "Car4_BEC" = "Car4+ BEC"
  )

data_path <- paste0("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_")
proportions_path <- paste0("visium_hd_lung_ila010/Visium_HD_Lung_ILA010_")
ext <- "_filtered100UMI"
plot_path <- paste0("visium_hd_lung_ila010/plots/")

umap_theme <- theme_classic(base_size=8) +
  theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6))