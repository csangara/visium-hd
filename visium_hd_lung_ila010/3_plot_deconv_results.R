use_color_palette <- "" # Use Karen's or my color palette
source("visium_hd_lung_ila010/0_utils.R")
source("utils/scatterbarplot_function.R")
bin_size <- 8
bin_size_str <- sprintf("%03dum", bin_size)

#### PREPROCESS DATA ####
# First run
# Normalize data & cluster data
# visium_obj <- NormalizeData(visium_obj) %>% FindVariableFeatures() %>% ScaleData() %>% 
#   RunPCA() %>% RunUMAP(dims = 1:30)
# visium_obj <- FindNeighbors(visium_obj, dims = 1:30) %>% FindClusters(resolution = 0.5)

# Save
# saveRDS(visium_obj, paste0(data_path, bin_size_str, ext, "_pp.rds"))

visium_obj <- readRDS(paste0(data_path, bin_size_str, ext, "_pp.rds"))

#### SIMPLE COOCCURRENCE HEATMAP ####
save_cooccurrence_plot <- FALSE

visium_obj <- readRDS(paste0(data_path, bin_size_str, ext, "_pp.rds"))
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

if (save_cooccurrence_plot){
  pheatmap::pheatmap(cooccurrence %>% `row.names<-`(proper_celltype_names) %>% 
                       `colnames<-`(proper_celltype_names),
                     color = RColorBrewer::brewer.pal(9, "Blues"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     filename = paste0(plot_path, "cooccurrence_heatmap_", bin_size_str, ".pdf"))
} else {
  p_heatmap <- ggplotify::as.ggplot(
    pheatmap::pheatmap(cooccurrence %>% `row.names<-`(proper_celltype_names) %>% 
                         `colnames<-`(proper_celltype_names),
                       color = RColorBrewer::brewer.pal(9, "Blues"),
                       cluster_rows = FALSE, cluster_cols = FALSE)
)
}

# Subset
cts_oi <- c("Bcells", "GC_Bcells", "Plasma_cells")
cooccurrence[cts_oi,] %>% 
  `row.names<-`(proper_celltype_names[cts_oi]) %>% 
  `colnames<-`(proper_celltype_names) %>% 
  pheatmap::pheatmap(color = RColorBrewer::brewer.pal(9, "Blues"),
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     #filename = paste0(plot_path, "cooccurrence_heatmap_subset_", bin_size_str, ".pdf"),
                     height = 2)

# Cooccurrence based on expected ratio (lift)
# Similar to Fisher's test
celltype_abundance <- colSums(deconv_props > 0)[names(color_palette)]
celltype_abundance_frac <- colSums(deconv_props > 0)[names(color_palette)]/nrow(deconv_props)
cooccurence_frac <- cooccurrence/nrow(deconv_props)
cooccurrence_scaled <- sweep(cooccurrence, 2, celltype_abundance, FUN = '/')

lift_matrix <- matrix(, nrow = length(color_palette), ncol = length(color_palette),
       dimnames = list(names(color_palette), names(color_palette)))
for (celltype1 in names(color_palette)){
  for (celltype2 in names(color_palette)){
    if (celltype1 == celltype2) next;
    
    p1 <- celltype_abundance[celltype1]
    p2 <- celltype_abundance[celltype2]
    p12 <- cooccurence_frac[celltype1,celltype2]
    
    lift_matrix[celltype1, celltype2] <- p12/(p1*p2)
  }
}

pheatmap::pheatmap(lift_matrix[1:3,], color = RColorBrewer::brewer.pal(9, "Blues"),
                   cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap::pheatmap(cooccurrence_scaled[1:3,], color = RColorBrewer::brewer.pal(9, "Blues"),
                   cluster_rows = FALSE, cluster_cols = FALSE)


#### DECONVOLUTION BARPLOT ####
bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)

visium_obj <- readRDS(paste0(data_path, bin_size_str, ext, "_pp.rds"))
deconv_props <- read.table(paste0(proportions_path, bin_size_str, ext,
                                  "/proportions_rctd_Visium_HD_Lung_ILA010_",
                                  bin_size_str, ext),
                           header = TRUE)

# Assign cell type with max proportion to each spot
visium_obj$celltype <- factor(colnames(deconv_props)[max.col(deconv_props)],
                              levels = names(color_palette))
visium_obj$celltype %>% table %>% sort

# Just for general use
deconv_props_df <- deconv_props %>%
  `rownames<-`(colnames(visium_obj)) %>% 
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

#### UMAP PLOTS ####
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
  umap_theme +
  theme(legend.key.size = unit(0.2, "cm"))
ggsave(paste0(plot_path, "umap/umap_celltype_", bin_size_str, use_color_palette, ".pdf"),
       p_umap_celltype, 
       width = 8, height = 8, bg = "white")

# For reference
umap_median <- umap_meta_df %>% 
  group_by(celltype) %>% 
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))
p_umap_celltype_text <- p_umap_celltype + geom_label(data = umap_median,
                                                     aes(x = umap_1, y = umap_2, label = celltype, fill = celltype),
                                                     color = "black", size=3) +
  scale_fill_manual(values = color_palette)
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
  umap_theme +
  theme(legend.position = "none")
ggsave(paste0(plot_path, "umap/umap_leiden_", bin_size_str, ".pdf"),
       p_umap_leiden, 
       width = 8, height = 8, bg = "white")

# Combind plots for Supplementary Figure 4.15
design <- "AABB
           CCCC"
p <- p_umap_celltype +
  theme(legend.position = "none") +
  p_barplot + 
  p_heatmap + plot_layout(design = design)
p
ggsave(p,
       filename = paste0(plot_path, "supplementary_spatial_deconv_overview_",
                         bin_size_str, "_bc2.pdf"),
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

color_scales <- list(
  "viridis" = scale_color_viridis_c(limits=c(0,1)),
  "RdYlBu" = scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")[2:10]),
                        limits=c(0,1)),
  "RdPuBu1" = scale_color_gradientn(colors = rev(c("#FF0000", "#CC0033", "#990066", "#660099", "#3300CC", "#0000FF")),
                        limits=c(0,1)),
  "RdPuBu2" =  scale_color_gradient2(low="blue", high="red",midpoint=0.5,
                          mid="#660099", limits=c(0,1))
)

color_scale_use <- "RdYlBu"
p_feature_props <- lapply(names(props_vec), function(props_name){
  ggplot(umap_meta_df %>% arrange(!!sym(props_name)) %>% 
           mutate(prop = !!sym(props_name))) +
    geom_point(aes(x = umap_1, y = umap_2, color = prop),
               size=0.4, stroke=0, shape=16) +
    color_scales[[color_scale_use]] +
    ggtitle(paste0(props_vec[props_name], " proportions")) +
    umap_theme
}) %>% setNames(props_vec)

props_feature_plot_path <- paste0(plot_path, "umap/", color_scale_use, "/")
if (!dir.exists(props_feature_plot_path)){
  dir.create(props_feature_plot_path)
}

for (ct in cts_oi){
  ggsave(paste0(props_feature_plot_path, "umap_", ct, "_props_",  bin_size_str, ".pdf"),
         p_feature_props[[ct]], 
         width = 8, height = 8, bg = "white")
}

for (ct_pair in ct_pairs_oi){
  p1 <- p_feature_props[[ct_pair[1]]]
  p2 <- p_feature_props[[ct_pair[2]]]
  p_combined <- p1 + p2 + plot_layout(ncol=2, guides="collect") &
    theme(legend.position = "bottom")
  ggsave(paste0(props_feature_plot_path, "umap_props_", ct_pair[1], "_", ct_pair[2], "_", bin_size_str, ".pdf"),
         p_combined,
         width = 8, height = 5, bg = "white")
}


# Ligand receptor expression on UMAP
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

for (ligrec in ligrec_pairs[4]){
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
      ggtitle(lig)
  })
  
  feature_plots_wrapped <- wrap_plots(feature_plots, nrow=1) &
    (umap_theme +
    theme(legend.key.size = unit(0.4, "cm")))
  ggsave(paste0(plot_path, "umap/ligand_receptor_expr/umap_expr_", paste0(ligrec, collapse="_"), "_", bin_size_str, ".pdf"),
         feature_plots_wrapped,
         width = 8, height = 4, bg = "white")
}

#### SPATIAL SCATTERBARPLOTS ####
roi_ext <- "inset_"
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

# To get coordinates
# ISpatialDimPlot(visium_obj)

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

plot_scatterbar(deconv_props, square_size=square_size,
                visium_obj_roi = visium_obj_roi) +
  scale_fill_manual(values = color_palette,
                    labels = proper_celltype_names,
                    breaks = names(color_palette)) +
  theme(legend.text = element_text(size=5),
        legend.key.size = unit(0.3, "cm"))
ggsave(paste0(plot_path, "spatial_proportions/spatial_scatterbarplot_",
              roi_ext, bin_size_str, ".pdf"),
       p_scatterbar,
       width = width, height = height, bg = "white")

#### SPATIAL CO-OCCURRENCE PLOT ####
ct_pairs_oi <- list(
  c("Fibroblasts", "Plasma_cells"),
  c("AT2_cells", "Plasma_cells"),
  c("Macrophages", "Plasma_cells"),
  c("Dendriticcells", "Plasma_cells"),
  c("Fibroblasts", "Macrophages"),
  c("Fibroblasts", "Dendriticcells")
)

deconv_props_roi <- deconv_props %>%
  `rownames<-`(colnames(visium_obj)) %>% 
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
         square_size = square_size)

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
  ggsave(paste0(plot_path, "spatial_proportions/cooccurrence/spatial_cooccurrence_",
                ct_pair[1], "_", ct_pair[2], "_", roi_ext, bin_size_str, ".pdf"),
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


#### (Rudimentary) niche detection based on KNN ####
# We want to see the fibroblast niche
# Try to cluster based on deconvolution proportions to see if there is a 'niche'
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

# Plot average fibroblast proportions on cluster
fibroblast_props_df <- visium_obj@meta.data %>% 
  mutate(fibroblast_props = deconv_props[, "Fibroblasts"]) 
cluster_order <- fibroblast_props_df %>% 
  group_by(Spatial.008um_snn_res.0.5) %>%
  summarise(mean_fibro = mean(fibroblast_props)) %>% 
  arrange(desc(mean_fibro)) %>%
  pull(Spatial.008um_snn_res.0.5) %>%
  as.character()

ggplot(fibroblast_props_df %>% 
         mutate(Spatial.008um_snn_res.0.5 = factor(Spatial.008um_snn_res.0.5,
                                                   levels = cluster_order)),
       aes(x=Spatial.008um_snn_res.0.5, y=fibroblast_props)) +
  geom_boxplot() +
  theme_classic()

visium_obj_top5_fib_clusters <- visium_obj %>% 
  subset(Spatial.008um_snn_res.0.5 %in% c(cluster_order[1:5], "4"))
SpatialDimPlot(visium_obj_top5_fib_clusters, group.by="Spatial.008um_snn_res.0.5",
               image.alpha = 0, stroke=NA, pt.size.factor = 6)
