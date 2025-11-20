library(Seurat)
library(tidyverse)
library(patchwork)
library(sf)
library(ggnewscale)
library(ggplotify)

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

#### SIMPLE COOCCURRENCE HEATMAP ####
bin_size <- 8
save_cooccurrence_plot <- FALSE
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
  
  cooccurrence <- cooccurrence %>% `row.names<-`(proper_celltype_names) %>% 
    `colnames<-`(proper_celltype_names)
  
  if (save_cooccurrence_plot){
    pheatmap::pheatmap(cooccurrence, color = RColorBrewer::brewer.pal(9, "Blues"),
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       filename = paste0(plot_path, "cooccurrence_heatmap_", bin_size_str, ".pdf"))
  } else {
  
    p_heatmap <- ggplotify::as.ggplot(
      pheatmap::pheatmap(cooccurrence, color = RColorBrewer::brewer.pal(9, "Blues"),
                         cluster_rows = FALSE, cluster_cols = FALSE)
  )
  }
  
  
}

# First run
# Normalize data & cluster data
# visium_obj <- NormalizeData(visium_obj) %>% FindVariableFeatures() %>% ScaleData() %>% 
#   RunPCA() %>% RunUMAP(dims = 1:30)
# visium_obj <- FindNeighbors(visium_obj, dims = 1:30) %>% FindClusters(resolution = 0.5)

# Save
# saveRDS(visium_obj, paste0(data_path, bin_size_str, ext, "_pp.rds"))

visium_obj <- readRDS(paste0(data_path, "008um", ext, "_pp.rds"))


#### PLOT DECONVOLUTION BARPLOT ####
bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)

visium_obj <- readRDS(paste0(data_path, bin_size_str, ext, "_pp.rds"))
deconv_props <- read.table(paste0(proportions_path, bin_size_str, ext,
                                  "/proportions_rctd_Visium_HD_Lung_ILA010_",
                                  bin_size_str, ext),
                           header = TRUE)

# Just for general use
deconv_props_df <- deconv_props %>% `rownames<-`(colnames(visium_obj)) %>% 
  rownames_to_column("spot") %>%
  pivot_longer(cols = -spot, names_to = "celltype",
               values_to = "proportion")
# Barplot
deconv_props_summ <- deconv_props %>%
  pivot_longer(cols = everything(), names_to = "celltype", values_to = "proportion") %>% 
  group_by(celltype) %>% 
  summarise(agg_proportion = mean(proportion)) %>% 
  mutate(celltype = factor(celltype, levels = names(color_palette)))


p_barplot <- ggplot(deconv_props_summ, aes(x = agg_proportion, y = forcats::fct_rev(celltype), fill = celltype)) +
  geom_bar(stat = "identity", width = 1) +
  # reverse fill order
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(labels = proper_celltype_names) +
  labs(x = "Mean proportion across all bins") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
p_barplot
ggsave(paste0(plot_path, "deconv_barplot_", bin_size_str, ".pdf"),
       plot = p_barplot, width = 4, height = 4, bg = "white")


# Barplot from scRNA-seq data
sc_data <- readRDS("data/scref_Lung_UBla/lung_combined_macs_DCs_seurat_annot_ident.rds")
sc_data@meta.data %>% select(annot_ident) %>% 
  group_by(annot_ident) %>% 
  summarise(n = n()) %>%
  mutate(prop = n / sum(n),
         # Remove space from annot_ident
         annot_ident = gsub(" ", "", annot_ident)) %>% 
  mutate(annot_ident = factor(annot_ident, levels = names(color_palette))) %>% 
  ggplot(aes(x = prop, y = forcats::fct_rev(annot_ident), fill = annot_ident)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = color_palette) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(labels = proper_celltype_names) +
  labs(x = "Proportion in scRNA-seq data") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# Assign cell type with max proportion to each spot
visium_obj$celltype <- factor(colnames(deconv_props)[max.col(deconv_props)],
                              levels = names(color_palette))
visium_obj$celltype %>% table %>% sort

# Instead do clustering on deconv proportions
deconv_kmeans <- mbkmeans::mbkmeans(as.matrix(t(deconv_props)), clusters=10)
deconv_kmeans_20 <- mbkmeans::mbkmeans(as.matrix(t(deconv_props)), clusters=20)
visium_obj$kmeans_10 <- deconv_kmeans$Clusters
visium_obj$kmeans_20 <- deconv_kmeans_20$Clusters

# Plot proportion per kmean cluster
deconv_props %>% 
  mutate(kmeans = deconv_kmeans_20$Clusters) %>%
  pivot_longer(cols = -kmeans, names_to = "celltype", values_to = "proportion") %>% 
  group_by(kmeans, celltype) %>%
  summarise(agg_proportion = mean(proportion)) %>%
  mutate(celltype = factor(celltype, levels = names(color_palette))) %>% 
  # Stacked bar plot
  ggplot(aes(x = factor(kmeans), y = agg_proportion, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_classic()

SpatialDimPlot(visium_obj, group.by="kmeans_10",
               image.alpha = 0, stroke=NA, pt.size.factor = 6)
DimPlot(visium_obj, group.by="kmeans_20",
        label=TRUE, label.size = 3)

### UMAP PLOTS ####
visium_obj$fibroblast_props <- deconv_props$Fibroblasts
visium_obj$plasma_props <- deconv_props$Plasma_cells
visium_obj$macrophage_props <- deconv_props$Macrophages
visium_obj$dendritic_props <- deconv_props$Dendriticcells

umap_meta_df <- Embeddings(visium_obj, "umap") %>% as.data.frame() %>%
  rownames_to_column("spot") %>% 
  inner_join(visium_obj@meta.data %>% rownames_to_column("spot"),
             by = "spot")

umap_median <- umap_meta_df %>% 
  group_by(celltype) %>% 
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

p_umap_celltype <- ggplot(umap_meta_df) +
  geom_point(aes(x = umap_1, y = umap_2, color = celltype), 
            size=0.2, stroke=0, shape=16) +
  scale_color_manual(values = color_palette) +
  guides(color = guide_legend(override.aes = list(size=1.5))) +
  ggtitle(paste0("UMAP of ", bin_size, "\u00b5m bins colored by most abundant celltype")) +
  theme_classic(base_size=8) +
  theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"))

# For reference
p_umap_celltype_text <- p_umap_celltype + geom_label(data = umap_median,
             aes(x = umap_1, y = umap_2, label = celltype, fill = celltype),
             color = "black", size=3) +
  scale_fill_manual(values = color_palette)

ggsave(paste0(plot_path, "umap/umap_celltype_", bin_size_str, ".pdf"),
       p_umap_celltype, 
       width = 8, height = 8, bg = "white")
ggsave(paste0(plot_path, "umap/umap_celltype_with_labels_", bin_size_str, ".pdf"),
       p_umap_celltype_text, 
       width = 8, height = 8, bg = "white")

p_umap_leiden <- ggplot(umap_meta_df) +
  geom_point(aes(x = umap_1, y = umap_2, color = as.character(Spatial.008um_snn_res.0.5)),
             size=0.4, stroke=0, shape=16) +
  geom_label(data = umap_meta_df %>% group_by(Spatial.008um_snn_res.0.5) %>%
                    summarise(umap_1 = median(umap_1), umap_2 = median(umap_2)),
                    aes(x = umap_1, y = umap_2,
                        label = Spatial.008um_snn_res.0.5, fill = Spatial.008um_snn_res.0.5),
                    color = "black", size=3) +
  ggtitle(paste0("UMAP of ", bin_size, "\u00b5m bins colored by Leiden clusters at 0.5 resolution")) +
  theme_classic(base_size=8) +
  theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
ggsave(paste0(plot_path, "umap/umap_leiden_", bin_size_str, ".pdf"),
       p_umap_leiden, 
       width = 8, height = 8, bg = "white")

props_vec <- c(fibroblast_props = "Fibroblast", plasma_props = "Plasma cell",
                   macrophage_props = "Macrophage", dendritic_props = "Dendritic cell")

p_feature_props <- lapply(names(props_vec), function(props_name){
  ggplot(umap_meta_df %>% arrange(!!sym(props_name)) %>% 
           mutate(prop = !!sym(props_name))) +
    geom_point(aes(x = umap_1, y = umap_2, color = prop), size=0.4, stroke=0, shape=16) +
    scale_color_viridis_c(limits=c(0,1)) +
    ggtitle(paste0(props_vec[props_name], " proportions")) +
    theme_classic(base_size=8) +
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 6))
  
})

p_feature_fibro_plasma <- p_feature_props[[1]] + p_feature_props[[2]] +
  plot_layout(ncol=2, guides="collect") &
  theme(legend.position = "bottom")
ggsave(paste0(plot_path, "umap/umap_props_fibro_plasma_", bin_size_str, ".pdf"),
       p_feature_fibro_plasma,
       width = 8, height = 6, bg = "white")

p_feature_macrophage_dendritic <- p_feature_props[[3]] + p_feature_props[[4]] +
  plot_layout(ncol=2, guides="collect") &
  theme(legend.position = "bottom")
ggsave(paste0(plot_path, "umap/umap_props_macrophage_dendritic_", bin_size_str, ".pdf"),
       p_feature_macrophage_dendritic,
       width = 8, height = 6, bg = "white")


# Il6r doesn't exist so we use IL6ra instead -> Ask Ines
# macs_ligands <- c("Il1a", "Il1b")
fib_ligands <- c("Tnfsf13b", "Cxcl12", "Il6")
fib_receptors <- c("Il1r1", "Il1rap")
# asc_receptors <- c("Tnfrsf13b", "Tnfrsf13c", "Tnfrsf17", "Cxcr4", "Il6ra", "Il6st")
# ligrecs_oi <- list(macs_ligands, c(fib_ligands, fib_receptors), asc_receptors) %>% 
#   setNames(c("macs_ligs", "fibro_ligrecs", "asc_recs"))

ligrec_pairs <- list(
  c("Tnfsf13b", "Tnfrsf13b"),
  c("Tnfsf13b", "Tnfrsf13c"),
  c("Tnfsf13b", "Tnfrsf17"),
  c("Cxcl12", "Cxcr4"),
  c("Il6", "Il6ra"),
  c("Il6", "Il6st"),
  c("Il1b", "Il1r1"),
  c("Il1b", "Il1rap"),
  c("Il1a", "Il1r1"),
  c("Il1a", "Il1rap")
)
unique(unlist(ligrec_pairs)) %in% rownames(visium_obj)

for (ligrec in ligrec_pairs){
  # Plot expression of ligands receptors, etc.
  feature_plots <- lapply(ligrec, function(lig){
    umap_meta_df %>% 
      inner_join(GetAssayData(visium_obj)[c(lig),,drop=FALSE] %>% t() %>% as.data.frame() %>%
                   rownames_to_column("spot"),
                 by = "spot") %>%
      pivot_longer(cols = all_of(c(lig)), names_to = "gene", values_to = "expression") %>% 
      arrange(expression) %>%
      ggplot(aes(x=umap_1, y=umap_2, color=expression)) +
      geom_point(size=0.4, stroke=0, shape=16) +
      scale_color_viridis_c() +
      ggtitle(lig) +
      theme_classic(base_size=8) 
  })
  
  feature_plots_wrapped <- wrap_plots(feature_plots, nrow=1) &
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(hjust = 0.5, size=8, face="bold"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  ggsave(paste0(plot_path, "umap/ligand_receptor_expr/umap_expr_", paste0(ligrec, collapse="_"), "_", bin_size_str, ".pdf"),
         feature_plots_wrapped,
         width = 8, height = 4, bg = "white")
}


# Plot average fibroblast proportions on cluster
cluster_order <- visium_obj@meta.data %>% 
  group_by(Spatial.008um_snn_res.0.5) %>%
  summarise(mean_fibro = mean(fibroblast_props)) %>% 
  arrange(desc(mean_fibro)) %>% pull(Spatial.008um_snn_res.0.5) %>% as.character()

ggplot(visium_obj@meta.data %>% 
         mutate(Spatial.008um_snn_res.0.5 = factor(Spatial.008um_snn_res.0.5,
                                                  levels = cluster_order)),
       aes(x=Spatial.008um_snn_res.0.5, y=fibroblast_props)) +
  geom_boxplot() +
  theme_classic()

visium_obj_top5_fib_clusters <- visium_obj %>% 
  subset(Spatial.008um_snn_res.0.5 %in% c(cluster_order[1:5], "4"))
SpatialDimPlot(visium_obj_top5_fib_clusters, group.by="Spatial.008um_snn_res.0.5",
               image.alpha = 0, stroke=NA, pt.size.factor = 6)

# Calculate distances between fibroblasts and plasma cells
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
fibro_points <- GetTissueCoordinates(visium_obj) %>%
  right_join(deconv_props_df %>% filter(celltype == "Fibroblasts" & proportion > 0),
             by = c("cell" = "spot")) %>%
  st_as_sf(coords = c("x", "y"))

plasma_points <- GetTissueCoordinates(visium_obj) %>%
  right_join(deconv_props_df %>% filter(celltype == "Plasma_cells" & proportion > 0),
             by = c("cell" = "spot")) %>%
  st_as_sf(coords = c("x", "y"))

# For each fibroblast, find the nearest plasma cell
fibro_to_plasma <- st_nearest_feature(fibro_points, plasma_points)

# Calculate distance to nearest plasma cell
fibro_points$dist_to_plasma <- st_distance(fibro_points, plasma_points[fibro_to_plasma, ], by_element = TRUE)
fib_ligands_expr <- GetAssayData(visium_obj, layer="counts")[c(fib_ligands, fib_receptors),] %>% t() %>%
  as.data.frame() %>% rownames_to_column("spot")

fibro_points_df <- fibro_points %>%
  left_join(fib_ligands_expr, by = c("cell" = "spot")) %>% 
  pivot_longer(cols = all_of(c(fib_ligands, fib_receptors)), names_to = "gene", values_to = "expression") %>% 
  # Convert distance to micron (8um = 18.08939 pixels)
  mutate(dist_to_plasma_um = as.numeric(dist_to_plasma) * (bin_size / square_size))

ggplot(fibro_points_df %>%
         mutate(gene = factor(gene, levels = c(fib_ligands, fib_receptors))),
       aes(x = dist_to_plasma_um, y = expression)) +
  facet_wrap(~gene, scales = "free_y") +
  geom_jitter(size=0.4, stroke=0, shape=16) +
  labs(y="Gene expression (raw counts)", x = "Distance to nearest plasma cell (\u00b5m)") +
  theme_classic(base_size=8)
ggsave(paste0(plot_path, "ligand_expression_vs_dist_to_plasma_", bin_size_str, ".pdf"),
       width = 8, height = 5, bg = "white")

# To get coordinates
# ISpatialDimPlot(visium_obj)


#### GIOTTO SPATIAL PROXIMITY ###
bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)
roi_ext <- ""

if (roi_ext == "inset_"){
  roi <- c(xmin = 5700, xmax = 8000, ymin = 16800, ymax = 18300)
  width <- 6
  height <- 7
} else {
  roi <- c(xmin = 5200, xmax = 9100, ymin = 12200, ymax = 18800)
  width <- 11
  height <- 8
}

cells_roi <- GetTissueCoordinates(visium_obj) %>%
  filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)

visium_obj_roi <- visium_obj %>% .[, colnames(.) %in% cells_roi]
square_size <- visium_obj_roi@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

deconv_props_roi <- deconv_props %>% `rownames<-`(colnames(visium_obj)) %>% 
  .[cells_roi, ] %>% 
  rownames_to_column("cell_ID")
#
library(Giotto)
# visium_obj <- readRDS(paste0(data_path, "008um", ext, "_pp.rds"))

# deconv_props <- read.table(paste0("visium_hd_lung_ila010/proportions_rctd_lung_ila010_008um", ext, ".tsv"),
#                            header = TRUE, sep = "\t")

giotto_obj <- createGiottoObject(raw_exprs = as.matrix(GetAssayData(visium_obj_roi, layer="counts")),
                                 spatial_locs = GetTissueCoordinates(visium_obj_roi))

enrObj <- createSpatEnrObj(
  deconv_props_roi,
  name = "DWLS",
  method = "DWLS",
  spat_unit = "cell",
  feat_type = "rna",
  provenance = "cell",
  misc = list()
)

giotto_obj

# create spatial enrichment object
giotto_obj <- setGiotto(giotto_obj, enrObj)
giotto_obj <- createSpatialNetwork(giotto_obj, minimum_k = 0)

saveRDS(giotto_obj, paste0(data_path, "008um", ext, "_pp_giotto.rds"))

proximity_table <- cellProximityEnrichmentSpots(giotto_obj,
                             spatial_network_name = "Delaunay_network")

#### SPATIAL SCATTERBARPLOTS ####
roi_ext <- "inset_"
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

if (roi_ext == "inset_"){
  roi <- c(xmin = 5700, xmax = 8000, ymin = 16800, ymax = 18300)
  width <- 6
  height <- 7
} else {
  roi <- c(xmin = 5200, xmax = 9100, ymin = 12200, ymax = 18800)
  width <- 11
  height <- 8
}

cells_roi <- GetTissueCoordinates(visium_obj) %>%
  filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)

visium_obj_roi <- visium_obj %>% .[, colnames(.) %in% cells_roi]

deconv_props_roi <- deconv_props %>% `rownames<-`(colnames(visium_obj)) %>% 
  .[cells_roi, ] %>% 
  rownames_to_column("spot") %>%
  pivot_longer(cols = -spot, names_to = "celltype",
               values_to = "proportion") %>% 
  filter(proportion > 0) %>% 
  # merge with coordinates
  left_join(GetTissueCoordinates(visium_obj_roi), by = c("spot" = "cell")) %>% 
  # the coords are x and y centroid, so we want x1, y1, x2, y2 as corners of squares
  # get whether there are one or two cell types
  mutate(n_celltypes = n(), .by = "spot",
         square_size = square_size,
         bin_size = bin_size)

deconv_props_df_square <- deconv_props_roi %>%
  filter(n_celltypes == 1) %>% 
  mutate(group = spot) %>% 
  # Draw squares for each
  mutate(x1 = y - square_size / 2,
         y1 = x - square_size / 2,
         x2 = y + square_size / 2,
         y2 = x - square_size / 2,
         x3 = y + square_size / 2,
         y3 = x + square_size / 2,
         x4 = y - square_size / 2,
         y4 = x + square_size / 2) %>%
  rename(coord_x = x, coord_y = y) %>% 
  pivot_longer(cols = c(x1, y1, x2, y2, x3, y3, x4, y4),
               names_to = c(".value", "corner"),
               names_pattern = "(x|y)([1-4])")

# Create barplot in squares
deconv_props_df_barplot <- deconv_props_roi %>% 
  filter(n_celltypes == 2) %>% 
  group_by(spot) %>% arrange(spot, desc(proportion)) %>% 
  mutate(rank = row_number(),
         group = paste0(spot, "_", rank)) %>% 
  mutate(x1 = case_when(rank == 1 ~ y - (square_size / 2) + (proportion*square_size),
                        rank == 2 ~ y + (square_size / 2) - (proportion*square_size)),
         y1 = x - square_size / 2,
         x2 = case_when(rank == 1 ~ y - (square_size / 2) + (proportion*square_size),
                        rank == 2 ~ y + (square_size / 2) - (proportion*square_size)),
         y2 = x + square_size / 2,
         x3 = case_when(rank == 1 ~ y - square_size / 2,
                        rank == 2 ~ y + square_size / 2),
         y3 = x + square_size / 2,
         x4 = case_when(rank == 1 ~ y - square_size / 2,
                        rank == 2 ~ y + square_size / 2),
         y4 = x - square_size / 2
  ) %>% 
  rename(coord_x = x, coord_y = y) %>%
  pivot_longer(cols = c(x1, y1, x2, y2, x3, y3, x4, y4),
               names_to = c(".value", "corner"),
               names_pattern = "(x|y)([1-4])")

deconv_props_df_shapes <- bind_rows(deconv_props_df_square, deconv_props_df_barplot) %>% 
  mutate(celltype = factor(celltype, levels = names(color_palette)))

p_scatterbar <- ggplot(deconv_props_df_shapes,
       aes(x = x, y = y)) +
  geom_polygon(aes(fill = celltype, group = group), show.legend = TRUE) +
  # White border
  geom_tile(data = deconv_props_roi %>% distinct(x, y),
            aes(x = y, y = x), height = square_size, width = square_size,
            fill = NA, color = "white", inherit.aes = FALSE) +
  scale_fill_manual(values = color_palette) +
  theme_void(base_size=8) +
  scale_y_reverse() +
  coord_fixed() +
  guides(fill = guide_legend(ncol=1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"))

ggsave(paste0(plot_path, "spatial_proportions/spatial_scatterbarplot_", roi_ext, bin_size_str, ".pdf"),
       p_scatterbar,
       width = width, height = height, bg = "white")

# Barplot of abundance in roi
deconv_props_roi %>% ungroup() %>% 
  tidyr::complete(spot, celltype, fill=list(proportion=0)) %>%
  group_by(celltype) %>% 
  summarise(agg_proportion = mean(proportion)) %>% 
  mutate(celltype = factor(celltype, levels = names(color_palette))) %>%
  ggplot(aes(x = agg_proportion, y = celltype, fill = celltype)) +
  geom_bar(stat = "identity", width = 1) +
  # reverse fill order
  scale_fill_manual(values = color_palette) +
  theme_classic()

#### CO-OCCURRENCE PLOT: Plasma_cells and Fibroblasts ####
fibro_plasma <- deconv_props_roi %>%
  filter(celltype %in% c("Plasma_cells", "Fibroblasts")) %>% 
  # Group by spot and arrange by descending proportion
  group_by(spot) %>% arrange(spot, desc(proportion))

fibro_plasma_scaled <- fibro_plasma %>% filter(n() == 2) %>% 
  # get the higher proportion cell type
  slice_max(proportion, n=1) %>% ungroup() %>% 
  # If it's Fibroblasts, scale from 0 to 1, otherwise from 0 to -1
  mutate(scaled_proportion = ifelse(celltype == "Fibroblasts",
                                    scales::rescale(proportion, to = c(0, 1), from=c(0.5, 1)),
                                    scales::rescale(proportion, to = c(0, -1), from=c(0.5, 1))))

p_spatial_cooccurrence <- ggplot() +
  geom_tile(aes(fill = proportion, x=y, y=x), color = "white", height=square_size, width=square_size,
            data = fibro_plasma %>% filter(celltype == 'Fibroblasts', n() == 1)) +
  scale_fill_gradientn(name = "Fibroblasts", colors = RColorBrewer::brewer.pal(9, "Greens")[3:8]) +
  new_scale_fill() +
  geom_tile(aes(fill = proportion, x=y, y=x), color = "white", height=square_size, width=square_size,
            data = fibro_plasma %>% filter(celltype == 'Plasma_cells', n() == 1)) +
  scale_fill_gradientn(name = "Plasma cells", colors = RColorBrewer::brewer.pal(9, "Purples")[3:8]) +
  new_scale_fill() +
  geom_tile(aes(fill = scaled_proportion, x=y, y=x), color = "white", height=square_size, width=square_size,
            data = fibro_plasma_scaled) +
  scale_fill_gradientn(name = "Co-occurring", colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")[2:10]),
                       limits = c(-1, 1), labels=c("Plasma cells", "Equal", "Fibroblasts"), breaks=c(-1, 0, 1)) +
  theme_void() +
  scale_y_reverse() +
  coord_fixed() +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=6))
ggsave(paste0(plot_path, "spatial_proportions/spatial_cooccurrence_fibro_plasma_", roi_ext, bin_size_str, ".pdf"),
       p_spatial_cooccurrence,
       width = width, height = height, bg = "white")

### COMBINE PLOTS FOR SUPPLEMENTARY ###
library(patchwork)
bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)
visium_obj_ori <- readRDS(paste0(data_path, bin_size_str, ".rds"))
square_size <- visium_obj_ori@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
ncounts_features_df <- visium_obj_ori@meta.data %>% 
  rownames_to_column("cell") %>%
  left_join(GetTissueCoordinates(visium_obj_ori), by="cell") %>% 
  mutate(square_size = square_size)

# Get outliers
quant <- visium_obj_ori$nCount_Spatial.008um %>% 
  # Get 90th percentile
  quantile(0.99)

p_scatterbar <- ggplot(mapping = aes(x = y, y = x)) +
  geom_tile(mapping = aes(fill = label),
            height=square_size, width=square_size,
            data = ncounts_features_df %>% filter(nCount_Spatial.008um <= 100) %>% 
              mutate(label="\u2264100")) +
  scale_fill_manual(values = "lightgrey", name = NULL) +
  new_scale_fill() + 
  geom_tile(mapping = aes(fill = nCount_Spatial.008um),
            height=square_size, width=square_size,
            data = ncounts_features_df %>% filter(nCount_Spatial.008um > 100)) +
  scale_fill_viridis_c(limits=c(100, quant), oob=scales::squish,
                       breaks = c(seq(0, 1000, 200)),
                       name = "Total counts") +
  scale_y_reverse() +
  coord_fixed() +
  theme_void(base_size=8) +
  theme(legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"))
p_scatterbar
ggsave(p_scatterbar,
       filename = paste0(plot_path, "spatialplot_ncounts_", bin_size_str, ".pdf"),
       width = 11, height = 8, bg = "white")


# Boxplot - I don't like this anymore lol
# ori_meta_df <- visium_obj_ori@meta.data %>% 
#   mutate(
#     outlier_lwr = nCount_Spatial.008um < quantile(nCount_Spatial.008um, probs = 0.25) - IQR(nCount_Spatial.008um) * 1.5,
#     outlier_upr = nCount_Spatial.008um > quantile(nCount_Spatial.008um, probs = 0.75) + IQR(nCount_Spatial.008um) * 1.5
#   )
# p_boxplot <- ori_meta_df %>% 
#   ggplot(aes(y=nCount_Spatial.008um, x=0)) +
#   geom_violin(linewidth=0.15) +
#   geom_boxplot(linewidth=0.15, outlier.shape=NA, width=0.1) +
#   geom_point(data = function(x) subset(x, outlier_lwr | outlier_upr), 
#              position = 'jitter', shape=16, stroke=0, size=0.1, alpha=0.2) +
#   # Plot mean as a dot
#   stat_summary(fun = "mean", geom = "point", shape = 21, size = 0.75,
#                stroke=0.25, fill = "blue", color="blue",
#                show.legend = FALSE) +
#   theme_minimal(base_size=7) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(panel.grid.minor.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
#         axis.ticks.x = element_blank(),
#         axis.ticks.y = element_line(linewidth=0.25),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size=6),
#         axis.title = element_blank())
# p_boxplot
# ggsave(p_boxplot,
#        filename = paste0(plot_path, "boxplot_ncounts_", bin_size_str, ".pdf"),
#        width = 2, height = 4, bg = "white")

# Let's just go for a histogram
ncounts_features_df_long <- 
  ncounts_features_df %>% select(cell, contains("Spatial.008um")) %>%
  pivot_longer(cols = contains("Spatial.008um"), 
               values_to = "x") %>% 
  mutate(type = gsub("_Spatial.008um", "", name)) %>% 
  mutate(type = gsub("^n", "", type))

means_df <- 
  ncounts_features_df_long %>%
  group_by(type) %>%
  mutate(max_x = max(x)) %>%
  mutate(text_y = max(table(cut_width(x, width = unique(max_x)/50)))) %>%
  summarise(mean_x = prettyNum(round(mean(x)), big.mark=",", scientific=FALSE),
            n = prettyNum(n(), big.mark=",", scientific=FALSE),
            text_x = mean(x) + unique(max_x)/50*4, text_y = max(text_y))

set1_colors <- RColorBrewer::brewer.pal(8, "Set1")[2]

# Plot counts and features
p_histogram <- ggplot(ncounts_features_df_long, aes(x = x)) +
  geom_histogram(color="white", fill=set1_colors,
                 alpha=0.5, bins=50, position="identity", size=0.1) +
  geom_text(aes(x = text_x, y = text_y, color=set1_colors,
                label = paste0("Mean: ", mean_x, "\nn = ", n)),
            data = means_df, vjust = 1, hjust=0, size=2,
            show.legend = FALSE
  ) +
  facet_wrap(~ type, scales = "free",
             labeller = labeller( .multi_line=FALSE)) +
  scale_fill_manual(values = set1_colors) +
  scale_color_manual(values = set1_colors) +
  # Use scientific notation for y
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        # Add right spacing to legend title
        legend.title = element_text(margin = margin(r = 15), size=8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
  )

ggsave(p_histogram,
       filename = paste0(plot_path, "histogram_ncounts_features_", bin_size_str, ".pdf"),
       width = 6, height = 3, bg = "white")

design <- "AAAA
           AAAA
           #BB#"
first_row <- p_scatterbar + p_histogram + 
  plot_layout(design=design)
first_row
ggsave(first_row,
       filename = paste0(plot_path, "supplementary_spatial_deconv_overview_", bin_size_str, "_a.pdf"),
       width = 10, height = 6, bg = "white")


design <- "AB
           CC"
second_row <- p_umap_celltype + theme(legend.position = "none") +
  p_barplot + 
  p_heatmap + plot_layout(design = design)

all_rows <- first_row / second_row +
  plot_layout(heights = c(1.3, 2))
ggsave(all_rows,
       filename = paste0(plot_path, "supplementary_spatial_deconv_overview_", bin_size_str, ".pdf"),
       width = 8, height = 11, bg = "white")

## OLD IDEA 1: Heatmap of deviation ##
ggplot(fibro_plasma %>% mutate(deviation = ifelse(n() == 2,proportion - lead(proportion), proportion)) %>%
         ungroup() %>% distinct(spot, deviation, x, y, square_size), aes(x = y, y = x)) +
  geom_tile(aes(fill = deviation), color = "white", height=square_size, width=square_size) +
  scale_fill_continuous(trans='reverse') +
  theme_void() +
  scale_y_reverse() +
  coord_fixed() +
  ggtitle(paste0(bin_size, "\u00b5m: Plasma cells vs Fibroblasts")) +
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Prop. Deviation")

## OLD IDEA 2: Blended color map ##
# (Code from Seurat)
cols <- c('lightgrey', '#ff0000', '#00ff00')
color.matrix <- Seurat:::BlendMatrix(
  two.colors = cols[2:3],
  col.threshold = 0,
  negative.color = cols[1]
)
cols <- cols[2:3]
colors <- list(
  color.matrix[, 1],
  color.matrix[1, ],
  as.vector(x = color.matrix)
)
blended_expr <- Seurat:::BlendExpression(data = fibro_plasma %>% ungroup %>% 
                  tidyr::complete(spot, celltype, fill=list(proportion=0)) %>% 
                  select(spot, celltype, proportion) %>%
                  pivot_wider(names_from = celltype, values_from = proportion) %>% 
                  column_to_rownames("spot"))
cols.use <- as.numeric(x = as.character(x = blended_expr[, 3])) +1
cols.use <- colors[[3]][sort(x = unique(x = cols.use))] %>% 
  setNames(levels(blended_expr$Fibroblasts_Plasma_cells))

fibro_plasma_color <- fibro_plasma %>% select(spot, x, y, square_size) %>% ungroup() %>%  
  inner_join(blended_expr %>% rownames_to_column("spot") %>%
               select(-Fibroblasts, -Plasma_cells),
             by = "spot")

ggplot(fibro_plasma_color , aes(x = y, y = x)) +
  geom_tile(aes(fill = Fibroblasts_Plasma_cells), color = "white", height=square_size, width=square_size) +
  scale_fill_manual(values=cols.use) +
  theme_void() +
  scale_y_reverse() +
  coord_fixed() +
  ggtitle(paste0(bin_size, "\u00b5m: Plasma cells vs Fibroblasts")) +
  theme(legend.position = "none",
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5))

# Legend
Seurat:::BlendMap(color.matrix = color.matrix)

# Installing Giotto
# (Had to update Rcpp and install extra packages on system for terra)
# remotes::install_github('drieslab/GiottoVisuals', ref='R4.1.0')
# remotes::install_github('drieslab/GiottoClass', ref='R4.1.0')
# remotes::install_github('drieslab/GiottoUtils', ref='R4.1.0')
# devtools::install_github('drieslab/Giotto', ref='R4.1.0')
# remotes::install_github('drieslab/GiottoData')

g <- GiottoData::loadGiottoMini("visium")
enr <- runDWLSDeconv(gobject = g, sign_matrix = sign_matrix, return_gobject = FALSE)
