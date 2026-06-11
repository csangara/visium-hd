# Before the spatial co-occurence plot with three color scales,
# We also experimented with some other visualization ideas
use_color_palette <- ""
source("visium_hd_lung_ila010/0_utils.R")
bin_size <- 8
bin_size_str <- sprintf("%03dum", bin_size)

roi_ext <- "inset"
if (roi_ext == "inset"){
  roi <- c(xmin = 5700, xmax = 8000, ymin = 16800, ymax = 18300)
  width <- 6
  height <- 7
} else {
  roi <- c(xmin = 5200, xmax = 9100, ymin = 12200, ymax = 18800)
  width <- 11
  height <- 8
}

visium_obj <- readRDS(paste0(data_path, bin_size_str, ext, "_pp.rds"))
deconv_props <- read.table(paste0(proportions_path, bin_size_str, ext,
                                  "/proportions_rctd_Visium_HD_Lung_ILA010_",
                                  bin_size_str, ext),
                           header = TRUE)
rownames(deconv_props) <- colnames(visium_obj)

cells_roi <- GetTissueCoordinates(visium_obj) %>%
  filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)
visium_obj_roi <- visium_obj %>% .[, colnames(.) %in% cells_roi]
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

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

ct_pairs_oi <- list(
  c("Fibroblasts", "Plasma_cells"),
  c("AT2_cells", "Plasma_cells"),
  c("Macrophages", "Plasma_cells"),
  c("Dendriticcells", "Plasma_cells"),
  c("Fibroblasts", "Macrophages"),
  c("Fibroblasts", "Dendriticcells")
)
ct_pair <- ct_pairs_oi[[1]]
ct1_ct2 <- deconv_props_roi %>%
  filter(celltype %in% c(ct_pair[1], ct_pair[2])) %>% 
  # Group by spot and arrange by descending proportion
  group_by(spot) %>% arrange(spot, desc(proportion))

#### OLD IDEA 1: Heatmap of deviation ####
ggplot(ct1_ct2 %>% mutate(deviation = ifelse(n() == 2,proportion - lead(proportion), proportion)) %>%
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

blended_expr <- Seurat:::BlendExpression(data = ct1_ct2 %>% ungroup %>% 
                                           tidyr::complete(spot, celltype, fill=list(proportion=0)) %>% 
                                           select(spot, celltype, proportion) %>%
                                           pivot_wider(names_from = celltype, values_from = proportion) %>% 
                                           column_to_rownames("spot"))
cols.use <- as.numeric(x = as.character(x = blended_expr[, 3])) +1
cols.use <- colors[[3]][sort(x = unique(x = cols.use))] %>% 
  setNames(levels(blended_expr$Fibroblasts_Plasma_cells))

ct1_ct2_color <- ct1_ct2 %>% select(spot, x, y, square_size) %>% ungroup() %>%  
  inner_join(blended_expr %>% rownames_to_column("spot") %>%
               select(-!!sym(ct_pair[1]), -!!sym(ct_pair[2])),
             by = "spot")

ggplot(ct1_ct2_color , aes(x = y, y = x)) +
  geom_tile(aes(fill = Fibroblasts_Plasma_cells), color = "white",
            height=square_size, width=square_size) +
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
