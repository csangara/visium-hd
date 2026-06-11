use_color_palette <- ""
source("visium_hd_lung_ila010/0_utils.R")

# Which cell types express Pigr in the single-cell data?
seurat_obj <- readRDS("data/scref_Lung_UBla/lung_combined_macs_DCs_seurat_annot_ident.rds")

# Mostly AT2 cells, and some secretory cells
DimPlot(seurat_obj, label = TRUE) + FeaturePlot(seurat_obj, "Pigr")
VlnPlot(seurat_obj, "Pigr")

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

#### Visualize genes on UMAP ####
genes_oi <- c("Pigr", "Igha", "Ighm")
lapply(genes_oi, function(gene){
    p <- umap_meta_df %>% 
      inner_join(GetAssayData(visium_obj)[gene,,drop=FALSE] %>% t() %>% as.data.frame() %>%
                   rownames_to_column("spot"),
                 by = "spot") %>%
      pivot_longer(cols = all_of(gene), names_to = "gene", values_to = "expression") %>% 
      arrange(expression) %>%
      ggplot(aes(x=umap_1, y=umap_2, color=expression)) +
      geom_point(size=0.4, stroke=0, shape=16) +
      scale_color_gradient(low = "lightgrey", high = "blue") +
      ggtitle(gene) +
      umap_theme
  
  ggsave(paste0("visium_hd_lung_ila010/plots/pigr/featureplot_", gene, ".pdf"), p,
         width = 8, height = 8, bg = "white")
  ggsave(paste0("visium_hd_lung_ila010/plots/pigr/featureplot_", gene, ".png"), p,
         width = 8, height = 8, bg = "white")
})


#### Contingency tables ####
# How many AT2 cells expresses Pigr?
pigr_df <- data.frame(at2_props = deconv_props$AT2_cells,
                      pigr_expr = GetAssayData(visium_obj)["Pigr",])
fisher_pigr_df <- pigr_df %>% group_by(at2_presence = at2_props > 0,
                                       pigr_expression = pigr_expr > 0) %>% 
  summarise(n = n()) %>% 
  mutate(at2_presence = factor(at2_presence, levels = c(TRUE, FALSE)),
         pigr_expression = factor(pigr_expression, levels = c(FALSE, TRUE)))

fisher_pigr_df
ggplot(fisher_pigr_df,
       aes(x=at2_presence, y=pigr_expression, fill=n)) +
  geom_tile() +
  geom_text(aes(label=n))+
  theme_minimal() +
  scale_x_discrete(position="top") +
  scale_fill_gradient(low = "#F7FBFF", high = "#4292C6") +
  labs(fill="Count", x="AT2 cell predicted in bin", y = "Expresses Pigr") +
  coord_fixed() +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size=11))
ggsave("visium_hd_lung_ila010/plots/pigr/contingency_table_AT2_Pigr.pdf",
       width=4, height=3)

# How many ASCs expressed Igha or Ighm?
ig_df <-  data.frame(asc_props = deconv_props$Plasma_cells,
                     iga_expr = GetAssayData(visium_obj)["Igha",],
                     igm_expr = GetAssayData(visium_obj)["Ighm",])
fisher_ig_df <- bind_rows(ig_df %>% group_by(asc_props = asc_props > 0, ig_expr = iga_expr > 0) %>% 
                            summarise(n=n()) %>% mutate(gene_oi = "Igha"),
                          ig_df %>% group_by(asc_props =asc_props > 0, ig_expr = igm_expr > 0) %>% 
                            summarise(n=n()) %>% mutate(gene_oi = "Ighm"),
                          ig_df %>% group_by(asc_props =asc_props > 0, ig_expr = iga_expr > 0 | igm_expr > 0) %>% 
                            summarise(n=n()) %>% mutate(gene_oi = "Igha/Ighm")) %>% 
  mutate(asc_props = factor(asc_props, levels = c(TRUE, FALSE)),
         ig_expr = factor(ig_expr, levels = c(FALSE, TRUE)),
         gene_oi = factor(gene_oi, levels = c("Igha", "Ighm", "Igha/Ighm")))
ggplot(fisher_ig_df,
       aes(x=asc_props, y=ig_expr, fill=n)) +
  geom_tile() +
  geom_text(aes(label=n))+
  facet_wrap(~gene_oi, ncol = 1,
             strip.position = "left") +
  theme_minimal() +
  scale_x_discrete(position="top") +
  scale_fill_gradient(low = "#F7FBFF", high = "#4292C6") +
  labs(fill="Count", x="ASC cell predicted in bin", y = "Expresses gene") +
  coord_fixed() +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size=11),
        strip.placement = "outside")
ggsave("visium_hd_lung_ila010/plots/pigr/contingency_table_ASC_Ig.pdf",
       width=4, height=7)

# Bins with both Pigr+ AT2 cells AND Ig+ ASCs
at2_ascs <- deconv_props_df %>%
  filter(celltype %in% c("AT2_cells", "Plasma_cells"), proportion > 0) %>% 
  filter(n() == 2, .by = "spot") %>% 
  distinct(spot) %>% 
  left_join(GetAssayData(visium_obj, layer="counts")[c("Pigr", "Igha", "Ighm"),,drop=FALSE] %>% t() %>%
            as.data.frame() %>% rownames_to_column("spot"), by = "spot") %>% 
  group_by(pigr_expr = Pigr > 0, ig_expr = Igha > 0 | Ighm > 0) %>% 
  summarise(n = n()) %>% 
  mutate(pigr_expr = factor(pigr_expr, levels = c(TRUE, FALSE)),
         ig_expr = factor(ig_expr, levels = c(FALSE, TRUE)))
ggplot(at2_ascs,
       aes(x=pigr_expr, y=ig_expr, fill=n)) +
  geom_tile() +
  geom_text(aes(label=n))+
  theme_minimal() +
  scale_x_discrete(position="top", labels = c("Pigr+", "Pigr-")) +
  scale_y_discrete(labels = c("Ig-", "Ig+")) +
  scale_fill_gradient(low = "#F7FBFF", high = "#4292C6") +
  labs(fill="Count", x="Bins with AT2 cells", y = "Bins with ASCs") +
  coord_fixed() +
  theme(panel.grid = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size=11))
ggsave("visium_hd_lung_ila010/plots/pigr/contingency_table_AT2_ASC.pdf",
       width=4, height=3)

# Fisher test = no significance
fisher.test(data.frame(ig_TRUE = c(22, 846),
                       ig_FALSE = c(3, 72)))


#### DISTANCE CALCULATION ####
# Calculate distances between pigr+ AT2 and iga or igm plasma cells
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
at2_points <- GetTissueCoordinates(visium_obj) %>%
  right_join(deconv_props_df %>% filter(celltype == "AT2_cells" & proportion > 0),
             by = c("cell" = "spot")) %>%
  st_as_sf(coords = c("x", "y"))

plasma_points <- GetTissueCoordinates(visium_obj) %>%
  right_join(deconv_props_df %>% filter(celltype == "Plasma_cells" & proportion > 0),
             by = c("cell" = "spot")) %>%
  st_as_sf(coords = c("x", "y"))

# For each AT2, find the nearest plasma cell
at2_to_plasma <- st_nearest_feature(at2_points, plasma_points)

# Calculate distance to nearest plasma cell
at2_points$dist_to_plasma <- st_distance(at2_points, plasma_points[at2_to_plasma, ], by_element = TRUE)
ig_expr <-  GetAssayData(visium_obj, layer="counts")[c("Igha", "Ighm"),,drop=FALSE] %>% t() %>%
  as.data.frame() %>% rownames_to_column("spot")
at2_points <- bind_cols(at2_points,
          plasma_points[at2_to_plasma, ] %>% left_join(ig_expr, by = c("cell" = "spot")) %>% 
            rename(nearest_plasma_cell = cell, plasma_props = proportion) %>%
            select(-celltype) %>% st_drop_geometry())

genes_oi <- c("Pigr")
ligands_expr <- GetAssayData(visium_obj, layer="counts")[genes_oi,,drop=FALSE] %>% t() %>%
  as.data.frame() %>% rownames_to_column("spot")

at2_points_df <- at2_points %>%
  left_join(ligands_expr, by = c("cell" = "spot")) %>% 
  pivot_longer(cols = all_of(genes_oi), names_to = "gene", values_to = "expression") %>% 
  # Convert distance to micron (8um = 18.08939 pixels)
  mutate(dist_to_plasma_um = as.numeric(dist_to_plasma) * (bin_size / square_size))

at2_points_df 
ggplot(at2_points_df %>%
         mutate(gene = factor(gene, levels = genes_oi),
                Ig_binary = Igha > 0 | Ighm > 0),
       aes(x = dist_to_plasma_um, y = expression, color = Ig_binary)) +
  geom_jitter(size=0.7, stroke=0, shape=16) +
  scale_color_manual(values = c("lightgrey","#6A3D9A"),
                     labels = c("FALSE", "TRUE")) +
  scale_y_continuous(breaks=0:5) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(y="Pigr expression (raw counts) in bins with AT2 cells",
       x = "Distance to nearest ASC (\u00b5m)",
       color = "Igha or Ighm\nexpressed") +
  theme_classic(base_size=8) +
  theme(legend.title = element_text(hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9),
        legend.background = element_rect(color="black", linewidth=0.2))

ggsave(paste0("visium_hd_lung_ila010/plots/pigr/Pigr_expression_vs_dist_to_ASC.pdf"),
       width = 7, height = 6, bg = "white")

##### SPATIAL CO-OCCURRENCE PLOT ######
roi_ext <- "_inset"
if (roi_ext == "_inset"){
  roi <- c(xmin = 5700, xmax = 8000, ymin = 16800, ymax = 18300)
  width <- 6; height <- 7; rel_width <- 0.3
} else {
  roi <- c(xmin = 5200, xmax = 9100, ymin = 12200, ymax = 18800)
  width <- 11; height <- 8; rel_width <- 0.1
}

cells_roi <- GetTissueCoordinates(visium_obj) %>%
  filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
visium_obj_roi <- visium_obj %>% .[, colnames(.) %in% cells_roi]

scatter_at2_df <- deconv_props_df %>% filter(celltype == "AT2_cells", proportion > 0) %>% 
  left_join(GetAssayData(visium_obj, layer="counts")["Pigr",,drop=FALSE] %>% t() %>%
  as.data.frame() %>% rownames_to_column("spot"), by = "spot") %>% 
  mutate(Pigr_binary = ifelse(Pigr > 0, "Pigr_TRUE", "Pigr_FALSE")) %>% 
  inner_join(GetTissueCoordinates(visium_obj_roi),
            by = c("spot" = "cell"))

scatter_asc_df <- deconv_props_df %>% filter(celltype == "Plasma_cells", proportion > 0) %>% 
  left_join(GetAssayData(visium_obj, layer="counts")[c("Igha", "Ighm"),,drop=FALSE] %>% t() %>%
              as.data.frame() %>% rownames_to_column("spot"), by = "spot") %>% 
  mutate(Ig_binary = ifelse(Igha > 0 | Ighm > 0, "Ig_TRUE", "Ig_FALSE")) %>% 
  inner_join(GetTissueCoordinates(visium_obj_roi),
             by = c("spot" = "cell"))

both_df <- inner_join(scatter_at2_df %>% filter(Pigr_binary == "Pigr_TRUE"),
                      scatter_asc_df %>% filter(Ig_binary == "Ig_TRUE") %>% select(-x,-y),
                      by = "spot")
both_df

legend_order <- c("Pigr- AT2 cells", "Pigr+ AT2 cells", "Ig- ASCs", "Ig+ ASCs", "Both") %>%
  setNames(c("Pigr_FALSE", "Pigr_TRUE", "Ig_FALSE", "Ig_TRUE", "Pigr_TRUE,Ig_TRUE"))
color_palette <- c("#CAB2D6","#6A3D9A", "#FB9A99", "#E31A1C", "#33A02C") %>%
  setNames(names(legend_order))
p_scatterbar <- ggplot() +
  # Main barplot
  geom_tile(data = scatter_at2_df,
            aes(fill = Pigr_binary, x = y - (square_size/4), y = x),
            height = square_size, width = square_size/2,
            show.legend = TRUE) +
  geom_tile(data = scatter_asc_df,
            aes(fill = Ig_binary, x = y + (square_size/4), y = x),
            height = square_size, width = square_size/2,
            show.legend = TRUE) +
  geom_tile(data = both_df,
            aes(fill = paste0(Pigr_binary,",",Ig_binary), x = y, y = x),
            height = square_size, width = square_size,
            show.legend = TRUE) +
  # White borders
  geom_tile(data = GetTissueCoordinates(visium_obj_roi),
            aes(x = y, y = x), height = square_size, width = square_size,
            fill = NA, color = "gray90", inherit.aes = FALSE) +
  theme_void(base_size = 7) +
  scale_y_reverse() +
  coord_fixed() +
  scale_fill_manual(values = color_palette,
                    limits = names(legend_order),
                    labels = legend_order) +
  guides(fill = guide_legend(ncol=1)) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.3, "cm"))  

# Create legend showing AT2 cells at the left of the square
# and ASCs at the right side of the square
legend_square_size <- 1
spacing <- 0.1
p_legend_text_size <- 2.5

legend_df <- data.frame(categories = factor(names(legend_order), levels = names(legend_order)),
                        labels = legend_order) %>%
  mutate(x = case_when(grepl("Pigr", labels) ~ legend_square_size/4,
                       grepl("Ig", labels) ~ legend_square_size/4*3,
                       TRUE ~ legend_square_size/2),
         width =  case_when(labels == "Both" ~ legend_square_size,
                            TRUE ~ legend_square_size/2),
         y = as.numeric(categories) + (as.numeric(categories)-1)*spacing)


p_legend <- ggplot(legend_df, aes(x = x, y = y)) +
  # Cell type legend
  geom_tile(aes(fill=categories, width = width), height = legend_square_size, show.legend = FALSE) +
  geom_tile(data = legend_df, aes(x=0.5, y=y),
            width=legend_square_size, height=legend_square_size,
            color="gray70", fill=NA, linewidth=0.2) +
  geom_text(data = legend_df, 
            aes(label=labels, x=1.5), size=p_legend_text_size, hjust=0) +
  scale_fill_manual(values = color_palette) +
  xlim(0, 10) +
  theme_void() +
  scale_y_reverse() +
  coord_fixed()

# Combine plot and legend
cowplot::plot_grid(p_scatterbar, p_legend, rel_widths = c(1, rel_width))
ggsave(paste0("visium_hd_lung_ila010/plots/pigr/spatialgrid_pigr_ig_celltypes", roi_ext, ".png"),
       width=width, height = height, bg = "white")
ggsave(paste0("visium_hd_lung_ila010/plots/pigr/spatialgrid_pigr_ig_celltypes", roi_ext, ".pdf"),
       width=width, height = height, bg = "white")

##### SPATIAL CO-OCCURRENCE PLOT - LOG EXPRESSION ####
# Similar plot, but use raw data (log expr) but without cell type information
library(ggnewscale)
roi_ext <- "_inset"

# Can put in one or more Ig variants to check
# If >1, the log1p values of the expression will be added to each other
# ig_oi <- c("Ighg1")
ig_oi_list <- list("Igha",
                   "Ighm",
                   "Ighg1",
                   c("Igha", "Ighm"))

for (roi_ext in c("_inset", "")){
  if (roi_ext == "_inset"){
    roi <- c(xmin = 5700, xmax = 8000, ymin = 16800, ymax = 18300)
    width <- 6; height <- 7; rel_width <- 0.3
  } else {
    roi <- c(xmin = 5200, xmax = 9100, ymin = 12200, ymax = 18800)
    width <- 11;  height <- 8; rel_width <- 0.1
  }
  
  cells_roi <- GetTissueCoordinates(visium_obj) %>%
    filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)
  square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
  visium_obj_roi <- visium_obj %>% .[, colnames(.) %in% cells_roi]
  
  scatter_pigr_df <- GetAssayData(visium_obj, layer="counts")["Pigr",,drop=FALSE] %>% t() %>%
    as.data.frame() %>% rownames_to_column("spot") %>% 
    mutate(Pigr_binary = Pigr > 0,
           Pigr = log1p(Pigr)) %>% 
    inner_join(GetTissueCoordinates(visium_obj_roi),
               by = c("spot" = "cell"))
  
  for (ig_oi in ig_oi_list){
    scatter_ig_df <- GetAssayData(visium_obj, layer="counts")[ig_oi,,drop=FALSE] %>% t() %>%
      as.data.frame() %>% rownames_to_column("spot") %>% 
      mutate(Ig_binary = if_any(starts_with("Igh"), ~ . > 0)) %>% 
      pivot_longer(cols = starts_with("Igh"), names_to = "Igh") %>%
      mutate(value = log1p(value)) %>% 
      summarise(Ig = sum(value), .by = c(spot, Ig_binary)) %>% 
      inner_join(GetTissueCoordinates(visium_obj_roi),
                 by = c("spot" = "cell")) 
    
    both_df <- inner_join(scatter_pigr_df %>% filter(Pigr_binary),
                          scatter_ig_df %>% filter(Ig_binary) %>% select(-x,-y),
                          by = "spot") %>% 
      mutate(Pigr_Ig = Pigr + Ig)
    
    legend_ig <- paste0(paste0("log(", ig_oi, ")"), collapse = "+\n")
    p_scatterbar <- ggplot() +
      # Main barplot
      geom_tile(data = scatter_pigr_df,
                aes(fill = Pigr, x = y - (square_size/4), y = x),
                height = square_size, width = square_size/2,
                show.legend = TRUE) +
      scale_fill_gradient(low = "white", high =  "#6A3D9A", name = "log(Pigr)",
                          guide = guide_colourbar(order = 1)) +
      new_scale_fill() +
      geom_tile(data = scatter_ig_df,
                aes(fill = Ig, x = y + (square_size/4), y = x),
                height = square_size, width = square_size/2,
                show.legend = TRUE) +
      scale_fill_gradient(low = "white", high = "#E31A1C", name = legend_ig,
                          guide = guide_colourbar(order = 2)) +
      new_scale_fill() +
      geom_tile(data = both_df,
                aes(fill = Pigr_Ig, x = y, y = x),
                height = square_size, width = square_size,
                show.legend = TRUE) +
      scale_fill_gradient(low = "#C7E9C0", high = "#238B45", name = "All",
                          guide = guide_colourbar(order = 3)) +
      new_scale_fill() +
      # White borders
      geom_tile(data = GetTissueCoordinates(visium_obj_roi),
                aes(x = y, y = x), height = square_size, width = square_size,
                fill = NA, color = "gray90", inherit.aes = FALSE) +
      theme_void(base_size = 7) +
      scale_y_reverse() +
      coord_fixed() +
      theme(legend.text = element_text(size=5),
            legend.key.size = unit(0.3, "cm"))  
    #p_scatterbar
    
    ig_filename <- paste0(tolower(ig_oi), collapse = "_")
    ggsave(paste0("visium_hd_lung_ila010/plots/pigr/spatialgrid_pigr_", ig_filename, roi_ext, ".png"),
           width=width, height = height, bg = "white")
    ggsave(paste0("visium_hd_lung_ila010/plots/pigr/spatialgrid_pigr_", ig_filename, roi_ext, ".pdf"),
           width=width, height = height, bg = "white")
    
  }
}


#### Contingency tables - no cell type information ####
# Pigr vs Ig variants
ig_oi_names <- c("Igha", "Ighm", "Ighg1") %>% setNames(rep("Ig", 3))
expr_df <- GetAssayData(visium_obj, layer="counts")[c("Pigr", ig_oi_names),,drop=FALSE] %>%
  t %>% as.data.frame() %>% rownames_to_column("spot")
ig_oi <- "Igha"
for (ig_oi in ig_oi_names){
  print(ig_oi)
  pigr_ig_df <- expr_df %>% select(spot, Pigr, all_of(ig_oi)) %>% 
    rename(any_of(ig_oi_names)) %>% 
    group_by(Pigr_expr = Pigr > 0,
             Ig_expr = Ig > 0) %>% 
    summarise(n = n(), .groups = "drop") %>% 
    mutate(Pigr_expr = factor(Pigr_expr, levels = c(TRUE, FALSE)),
           Ig_expr = factor(Ig_expr, levels = c(FALSE, TRUE)))
  
  ggplot(pigr_ig_df,
         aes(x=Pigr_expr, y=Ig_expr, fill=n)) +
    geom_tile() +
    geom_text(aes(label=n))+
    theme_minimal() +
    scale_x_discrete(position="top") +
    scale_fill_gradient(low = "#F7FBFF", high = "#4292C6") +
    labs(fill="Count", x="Expresses Pigr", y = paste0("Expresses ", ig_oi)) +
    coord_fixed() +
    theme(panel.grid = element_blank(),
          legend.key.size = unit(0.5, "cm"),
          legend.title = element_text(size=11))
  ggsave(paste0("visium_hd_lung_ila010/plots/pigr/contingency_table_Pigr_", ig_oi, ".pdf"),
         width=4, height=3)
}