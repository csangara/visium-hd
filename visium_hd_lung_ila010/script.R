library(Seurat)
library(tidyverse)
library(patchwork)
library(sf)
library(ggnewscale)
use_color_palette <- "_k"

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

bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)
ext <- "_filtered100UMI"
plot_path <- paste0("visium_hd_lung_ila010/plots/")
first_run <- FALSE


if (first_run){
    visium_obj <- Load10X_Spatial(data.dir = "data/Visium_HD_Lung_ILA010/", bin.size = 8)
    bins_to_keep <- visium_obj@meta.data[paste0("nCount_Spatial.", bin_size_str)] > 100
    visium_obj <- visium_obj[, bins_to_keep]
    
    # Normalize data & cluster data
    visium_obj <- NormalizeData(visium_obj) %>% FindVariableFeatures() %>% ScaleData() %>% 
      RunPCA() %>% RunUMAP(dims = 1:30)
    visium_obj <- FindNeighbors(visium_obj, dims = 1:30) %>% FindClusters(resolution = 0.5)

    saveRDS(visium_obj, paste0("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_", bin_size_str, ext, "_pp.rds"))
} else {
    visium_obj <- readRDS(paste0("data/Visium_HD_Lung_ILA010/Visium_HD_Lung_ILA010_", bin_size_str, ext, "_pp.rds"))
}

#### CO-OCCURRENCE MATRIX ####
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

# Full heatmap
pheatmap::pheatmap(cooccurrence, color = RColorBrewer::brewer.pal(9, "Blues"),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   filename = paste0(plot_path, "cooccurrence_heatmap_", bin_size_str, ".pdf"))

# Subset
cooccurrence[c("B cells", "GC B cells", "ASCs"),] %>% 
  pheatmap::pheatmap(color = RColorBrewer::brewer.pal(9, "Blues"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     filename = paste0(plot_path, "cooccurrence_heatmap_subset_", bin_size_str, ".pdf"),
                     height = 2)

#### DECONVOLUTION BARPLOT ####
proportions_path <- paste0("visium_hd_lung_ila010/Visium_HD_Lung_ILA010_")

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
  scale_y_discrete(labels = function(x) proper_celltype_names[x]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Mean proportion across all bins") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
p_barplot
ggsave(paste0(plot_path, "deconv_barplot_", bin_size_str,  use_color_palette, ".pdf"),
       plot = p_barplot, width = 4, height = 4, bg = "white")

# Assign cell type with max proportion to each spot
visium_obj$celltype <- factor(colnames(deconv_props)[max.col(deconv_props)],
                              levels = names(color_palette))
visium_obj$celltype %>% table %>% sort

#### UMAP PLOTS #####

umap_meta_df <- Embeddings(visium_obj, "umap") %>% as.data.frame() %>%
  rownames_to_column("spot") %>% 
  inner_join(visium_obj@meta.data %>% rownames_to_column("spot"),
             by = "spot")

p_umap_celltype <- ggplot(umap_meta_df) +
  geom_point(aes(x = umap_1, y = umap_2, color = celltype), 
             size=0.4, stroke=0, shape=16) +
  scale_color_manual(values = color_palette,
                     labels = function(x) proper_celltype_names[x]) +
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
ggsave(paste0(plot_path, "umap/umap_celltype_", bin_size_str, use_color_palette, ".pdf"),
       p_umap_celltype, 
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


#### PROPORTIONS ON UMAP ####
cts_oi <- c("Fibroblasts", "Plasma_cells", "Macrophages", "Dendriticcells", "AT2_cells")
ct_pairs_oi <- list(
  c("Fibroblasts", "Plasma_cells"),
  c("AT2_cells", "Plasma_cells"),
  c("Macrophages", "Plasma_cells"),
  c("Dendriticcells", "Plasma_cells"),
  c("Macrophages", "Fibroblasts"),
  c("Dendriticcells", "Fibroblasts")
)
for (ct in cts_oi){
  visium_obj <- AddMetaData(visium_obj, deconv_props[, ct], col.name = paste0(ct, "_props"))
}

umap_meta_df <- Embeddings(visium_obj, "umap") %>% as.data.frame() %>%
  rownames_to_column("spot") %>% 
  inner_join(visium_obj@meta.data %>% rownames_to_column("spot"),
             by = "spot")

props_vec <- cts_oi %>% setNames(paste0(cts_oi, "_props"))

p_feature_props <- lapply(names(props_vec), function(props_name){
  ggplot(umap_meta_df %>% arrange(!!sym(props_name)) %>% 
           mutate(prop = !!sym(props_name))) +
    geom_point(aes(x = umap_1, y = umap_2, color = prop), size=0.4, stroke=0, shape=16) +
    # scale_color_viridis_c(limits=c(0,1)) +
    scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")[2:10]),
                          limits=c(0,1)) +
    # scale_color_gradientn(colors = rev(c("#FF0000", "#CC0033", "#990066", "#660099", "#3300CC", "#0000FF")),
    #                         limits=c(0,1)) +
    # scale_color_gradient2(low="blue", high="red",midpoint=0.5,
    #                       mid="#660099", limits=c(0,1))+
    ggtitle(paste0(proper_celltype_names[props_vec[props_name]], " proportions")) +
    theme_classic(base_size=8) +
    theme(plot.title = element_text(hjust = 0.5, size=8, face="bold"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 6))
  
}) %>% setNames(props_vec)


for (ct in cts_oi){
  ggsave(paste0(plot_path, "umap/umap_", ct, "_props_",  bin_size_str, ".pdf"),
         p_feature_props[[ct]], 
         width = 8, height = 8, bg = "white")
  
}

for (ct_pair in ct_pairs_oi){
  p1 <- p_feature_props[[ct_pair[1]]]
  p2 <- p_feature_props[[ct_pair[2]]]
  p_combined <- p1 + p2 + plot_layout(ncol=2, guides="collect") &
    theme(legend.position = "bottom")
  ggsave(paste0(plot_path, "umap/umap_props_", ct_pair[1], "_", ct_pair[2], "_", bin_size_str, ".pdf"),
         p_combined,
         width = 8, height = 5, bg = "white")
}

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

#### SPATIAL SCATTERBARPLOTS ####
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
roi_ext <- "inset_" # "inset_" or ""

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
  scale_fill_manual(values = color_palette,
                    labels = function(x) proper_celltype_names[x]) +
  theme_void(base_size=8) +
  scale_y_reverse() +
  coord_fixed() +
  guides(fill = guide_legend(ncol=1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"))

ggsave(paste0(plot_path, "spatial_proportions/scatterbarplot/spatial_scatterbarplot_", roi_ext, bin_size_str,
              use_color_palette, ".pdf"),
       p_scatterbar,
       width = width, height = height, bg = "white")

#### SPATIAL CO-OCCURRENCE PLOT: Plasma_cells and Fibroblasts ####
ct_pairs_oi <- list(
  c("Fibroblasts", "Plasma_cells"),
  c("AT2_cells", "Plasma_cells"),
  c("Macrophages", "Plasma_cells"),
  c("Dendriticcells", "Plasma_cells"),
  c("Fibroblasts", "Macrophages"),
  c("Fibroblasts", "Dendriticcells")
)

for (ct_pair in ct_pairs_oi){
  ct1_ct2 <- deconv_props_roi %>%
    filter(celltype %in% c(ct_pair[1], ct_pair[2])) %>% 
    # Group by spot and arrange by descending proportion
    group_by(spot) %>% arrange(spot, desc(proportion))
  
  ct1_ct2_scaled <- ct1_ct2 %>% filter(n() == 2) %>% 
    # get the higher proportion cell type
    slice_max(proportion, n=1) %>% ungroup() %>% 
    # Scale one cell type from 0 to 1, otherwise from 0 to -1
    mutate(scaled_proportion = ifelse(celltype == ct_pair[1],
                                      scales::rescale(proportion, to = c(0, 1), from=c(0.5, 1)),
                                      scales::rescale(proportion, to = c(0, -1), from=c(0.5, 1))))
  
  p_spatial_cooccurrence <- ggplot() +
    geom_tile(aes(fill = proportion, x=y, y=x), color = "white", height=square_size, width=square_size,
              data = ct1_ct2 %>% filter(celltype == ct_pair[1], n() == 1)) +
    scale_fill_gradientn(name = proper_celltype_names[ct_pair[1]], colors = RColorBrewer::brewer.pal(9, "Greens")[3:8],
                         guide = guide_colourbar(order = 1)) +
    new_scale_fill() +
    geom_tile(aes(fill = proportion, x=y, y=x), color = "white", height=square_size, width=square_size,
              data = ct1_ct2 %>% filter(celltype == ct_pair[2], n() == 1)) +
    scale_fill_gradientn(name = proper_celltype_names[ct_pair[2]], colors = RColorBrewer::brewer.pal(9, "Purples")[3:8],
                         guide = guide_colourbar(order = 3)) +
    new_scale_fill() +
    geom_tile(aes(fill = scaled_proportion, x=y, y=x), color = "white", height=square_size, width=square_size,
              data = ct1_ct2_scaled) +
    scale_fill_gradientn(name = "Co-occurring", colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")[2:10]),
                         limits = c(-1, 1), breaks=c(-1, 0, 1),
                         labels=c(proper_celltype_names[ct_pair[2]], "Equal", proper_celltype_names[ct_pair[1]]),
                         guide = guide_colourbar(order = 2)) +
    theme_void() +
    scale_y_reverse() +
    coord_fixed() +
    theme(legend.key.size = unit(0.3, "cm"),
          legend.title = element_text(size=7),
          legend.text = element_text(size=6))
  ggsave(paste0(plot_path, "spatial_proportions/cooccurrence/spatial_cooccurrence_", ct_pair[1], "_", ct_pair[2], "_", roi_ext, bin_size_str, ".pdf"),
         p_spatial_cooccurrence,
         width = width, height = height, bg = "white")
  
}

#### DISTANCE TO NEAREST PLASMA CELL ####

square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

interactions_oi <- list(
  list(sender = c("Fibroblasts"), receiver = "Plasma_cells",
       genes_oi = c("Tnfsf13b", "Cxcl12", "Il6", "Il1r1", "Il1rap")),
  list(sender = c("Macrophages", "Dendriticcells"), receiver = "Fibroblasts",
       genes_oi = c("Il1a", "Il1b"))
)

for(i in 1:length(interactions_oi)){
  sender <- interactions_oi[[i]]$sender
  receiver <- interactions_oi[[i]]$receiver
  genes_oi <- interactions_oi[[i]]$genes_oi
  
  # Calculate distances between fibroblasts and plasma cells
  sender_points <- GetTissueCoordinates(visium_obj) %>%
    right_join(deconv_props_df %>% filter(celltype %in% sender & proportion > 0),
               by = c("cell" = "spot")) %>%
    st_as_sf(coords = c("x", "y"))
  
  receiver_points <- GetTissueCoordinates(visium_obj) %>%
    right_join(deconv_props_df %>% filter(celltype == receiver & proportion > 0),
               by = c("cell" = "spot")) %>%
    st_as_sf(coords = c("x", "y"))
  
  # For each sender, find the nearest receiver
  sender_to_receiver <- st_nearest_feature(sender_points, receiver_points)
  
  # Calculate distance to nearest receiver
  sender_points$dist_to_receiver <- st_distance(sender_points, receiver_points[sender_to_receiver, ], by_element = TRUE)
  
  genes_expr <- GetAssayData(visium_obj, layer="counts")[genes_oi,] %>% t() %>%
    as.data.frame() %>% rownames_to_column("spot")
  
  sender_points_df <- sender_points %>%
    left_join(genes_expr, by = c("cell" = "spot")) %>% 
    pivot_longer(cols = all_of(genes_oi), names_to = "gene", values_to = "expression") %>% 
    # Convert distance to micron (8um = 18.08939 pixels)
    mutate(dist_to_receiver_um = as.numeric(dist_to_receiver) * (bin_size / square_size))
  
  ggplot(sender_points_df %>%
           mutate(gene = factor(gene, levels = genes_oi),
                  expression = as.factor(expression)),
         aes(x = dist_to_receiver_um, y = expression)) +
    facet_wrap(~gene, scales = "free_y") +
    geom_jitter(size=0.4, stroke=0, shape=16) +
    labs(y=paste0("Raw counts in bins with ",
                  paste0(proper_celltype_names[interactions_oi[[i]]$sender], collapse="/")),
         x = paste0("Distance to nearest ", proper_celltype_names[receiver], " (\u00b5m)")) +
    theme_classic(base_size=8)
  
  # Calculate width and height based on number of genes
  n_genes <- length(genes_oi)
  ncol <- ifelse(n_genes <= 3, n_genes, ceiling(sqrt(n_genes)))
  nrow <- ceiling(n_genes / ncol)
  width <- ncol * 3
  height <- nrow * 3
  
  ggsave(paste0(plot_path, "distance/",  paste0(sender, collapse="_"), "_expr_vs_dist_to_", receiver, "_", bin_size_str, ".pdf"),
         width = width, height = height, bg = "white")
}
