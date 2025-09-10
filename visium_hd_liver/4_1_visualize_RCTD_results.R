library(Seurat)
library(tidyverse)

dataset <- "" #or "_caw009" or ""
data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
plot_path <- paste0("visium_hd_liver", dataset, "/plots/")


color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "CentralVeinEndothelialcells" = "#FED8B1FF",
                   "LSECs" = "#FBB05FFF",
                   "PortalVeinEndothelialcells" = "#CC7722FF",
                   "LymphaticEndothelialcells" = "#8F4716FF",
                   "Cholangiocytes" = "#C61B84FF",
                   "HsPCs" = "#F19FC3FF",
                   "Stellatecells" = "#A31A2AFF",
                   "Mesothelialcells" = "#D0110BFF",
                   "Fibroblasts" = "#E45466FF",
                   "CapsularFibroblasts" = "#D46F6CFF",
                   "Kupffercells" = "#5DA6DBFF",
                   "MonocytesMonocytederivedcells" = "#a3daf3",
                   "cDC1s" = "#893A86FF",
                   "cDC2s" = "#893A86FF",
                   "pDCs" = "#893A86FF",
                   "MigcDCs" = "#893A86FF",
                   "Bcells" = "#9C7EBAFF",
                   "NKcells" = "#4A6E34FF",
                   "Tcells" = "#3AB04AFF",
                   "ILC1s" = "#A3D7BAFF",
                   "Basophils" = "#191919",
                   "Neutrophils" = "#727272")

celltype_order <- names(color_palette)
                   
bin_size <- 8 # 8, 16, or 32

# Fill the text with 0
bin_size_str <- sprintf("%03dum", bin_size)

visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                             bin_size_str, ".rds"))
dim(visium_obj) # 19059 genes x 462269 spots

ext <- "_annot_cd45"

deconv_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                  "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                  bin_size_str, ext),
                           header = TRUE)
dim(deconv_props) #460942 spots

removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,paste0("nCount_Spatial.", bin_size_str)] < 100)]

# Check if removed rows + leftover rows == total rows (yes)
length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2]

# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]

dim(visium_obj_subset) # 19059 genes x 460942 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

# Assign barcode to most abundant cell type per spot
all(rownames(deconv_props) == colnames(visium_obj_subset))
visium_obj_subset$celltype <- colnames(deconv_props)[max.col(deconv_props)] %>% 
  factor(levels = celltype_order)
visium_obj_subset$celltype %>% table %>% sort

p_celltype <- SpatialDimPlot(visium_obj_subset, group.by = "celltype",
               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = color_palette) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"))
p_celltype
ggsave(paste0(plot_path, "spatialdimplot_celltype_", bin_size_str, ".png"), p_celltype,
       width = 10, height = 12, bg = "white")

# Let's get a region of interest
roi <- c(6000, 10000, 30000, 34000)
roi_name <- paste0("_", paste0(sprintf("%02d", roi/1000), collapse=""))

cells_roi <- GetTissueCoordinates(visium_obj_subset) %>%
  filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)

SpatialDimPlot(visium_obj_subset,
               cells.highlight = cells_roi)
ggsave(paste0(plot_path, "spatialdimplot_ROI_", roi_name, "_", bin_size_str, ".png"),
       width = 10, height = 8, bg = "white")

visium_obj_roi <- visium_obj_subset %>% .[, colnames(.) %in% cells_roi]

SpatialFeaturePlot(visium_obj_roi, "Glul", slot="counts",
                   stroke=NA, image.alpha=0, pt.size.factor=12)
ggsave(paste0(plot_path, "spatialfeatureplot_Glul_ROI", roi_name, "_",
              bin_size_str, ".png"),
       width = 8, height = 8, bg = "white")

SpatialFeaturePlot(visium_obj_roi, "Hal", slot="counts",
                   stroke=NA, image.alpha=0, pt.size.factor=12)
ggsave(paste0(plot_path, "spatialfeatureplot_Hal_ROI", roi_name, "_",
              bin_size_str, ".png"),
       width = 8, height = 8, bg = "white")


square_size <- visium_obj_subset@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
deconv_props_roi <- deconv_props[rownames(deconv_props) %in% cells_roi, ]

deconv_props_df <- deconv_props_roi %>% 
  rownames_to_column("spot") %>%
  pivot_longer(cols = -spot, names_to = "celltype",
               values_to = "proportion") %>% 
  filter(proportion > 0) %>% 
  # merge with coordinates
  left_join(GetTissueCoordinates(visium_obj_roi), by = c("spot" = "cell")) %>% 
  # the coords are x and y centroid, so we want x1, y1, x2, y2 as corners of squares
  # get whether there are one or two cell types
  mutate(n_celltypes = n(), .by = "spot")

# Pie plot of cell type proportions in the ROI
library(scatterpie)

deconv_props_scatterpie <- deconv_props_roi  %>%
  rownames_to_column("spot") %>%
  left_join(GetTissueCoordinates(visium_obj_roi), by = c("spot" = "cell")) %>% 
  mutate(radius = square_size / 2)

p_scatterpie <- ggplot() +
  geom_scatterpie(aes(r=radius, x = y, y = x, group=spot),
                  data = deconv_props_scatterpie,
                  cols=unique(deconv_props_df$celltype), color=NA) +
  scale_fill_manual(values = color_palette) +
  coord_equal() +
  scale_y_reverse() +
  theme_void()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(paste0(plot_path, "spatialscatterpie_ROI", roi_name, "_", bin_size_str, ".png"),
       p_scatterpie, width = 10, height = 8, bg = "white")

deconv_props_df_square <- deconv_props_df %>%
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

deconv_props_df_tri <- deconv_props_df %>% 
  filter(n_celltypes == 2) %>% 
  group_by(spot) %>% arrange(spot, desc(proportion)) %>% 
  mutate(rank = row_number(),
         group = paste0(spot, "_", rank)) %>% 
  # Select only lower proportion
  # Draw half triangles for each
  # Higher proportion in the lower half (since we will flip it later)
  mutate(x1 = y + square_size / 2,
         y1 = x - square_size / 2,
         x2 = y - square_size / 2,
         y2 = x + square_size / 2,
         x3 = case_when(rank == 1 ~ y - square_size / 2,
                        rank == 2 ~ y + square_size / 2),
         y3 = case_when(rank == 1 ~ x - square_size / 2,
                        rank == 2 ~ x + square_size / 2)
         
         ) %>% 
  rename(coord_x = x, coord_y = y) %>%
  pivot_longer(cols = c(x1, y1, x2, y2, x3, y3),
               names_to = c(".value", "corner"),
               names_pattern = "(x|y)([1-3])")


deconv_props_df_shapes <- bind_rows(deconv_props_df_square, deconv_props_df_tri) %>% 
  mutate(celltype = factor(celltype, levels = celltype_order))

# Create the ggplot
ggplot(deconv_props_df_shapes, aes(x = x, y = y, group = group)) +
  geom_polygon(aes(fill = celltype)) +
  scale_fill_manual(values = color_palette) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))

ggsave(paste0(plot_path, "spatialtriplot_doublet_celltype_ROI", roi_name, "_",
              bin_size_str, ".png"),
       width = 10, height = 8, bg = "white")

# Create barplot in squares
deconv_props_df_barplot <- deconv_props_df %>% 
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
  mutate(celltype = factor(celltype, levels = celltype_order))

# Create the ggplot
ggplot(deconv_props_df_shapes, aes(x = x, y = y)) +
  geom_polygon(aes(fill = celltype, group=group)) +
  # White border
  geom_tile(data = deconv_props_df %>% distinct(x, y),
            aes(x = y, y = x), height = square_size, width = square_size,
            fill = NA, color = "white") +
  scale_fill_manual(values = color_palette) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))

ggsave(paste0(plot_path, "spatialbarplot_doublet_celltype_ROI", roi_name, "_",
              bin_size_str, ".png"),
       width = 10, height = 8, bg = "white")

SpatialDimPlot(visium_obj_roi,
               group.by = "celltype",
               image.alpha = 0, pt.size.factor = 14, shape=22, stroke=0.1) +
  scale_fill_manual(values = color_palette) +
  guides(fill = guide_legend(ncol=1, override.aes = list(size = 3))) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
  
#p_celltype
ggsave(paste0(plot_path, "spatialdimplot_celltype_ROI", roi_name, "_",
              bin_size_str, ".png"),
       width = 10, height = 8, bg = "white") 

deconv_props_df_all <- data.frame(deconv_props) %>% 
  pivot_longer(cols = everything(), names_to = "celltype",
               values_to = "proportion")

# Cell type proportions barplot
deconv_props_summ <- deconv_props_df_all %>% 
  group_by(celltype) %>% 
  summarise(agg_proportion = mean(proportion))

# Barplot ordered by total abundance
p <- ggplot(deconv_props_summ, aes(x=reorder(celltype, agg_proportion), y=agg_proportion)) +
  geom_bar(stat="identity") +
  # Add text of value, rounded to two digits, only if value is > 0
  geom_text(aes(label=ifelse(agg_proportion > 0.001, round(agg_proportion, 3), "")), nudge_y = 0.03, size=2) +
  coord_flip() +
  theme_minimal(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  ggtitle("Average cell type proportions across tissue")

p
ggsave(paste0(plot_path, "barplot_avg_across_tissue_", bin_size_str, ".png"), p,
       width = 8, height = 6, bg = "white")


# What is the distribution of cell types across the tissue?
p_boxplot_all <- ggplot(deconv_props_df %>% filter(proportion > 0.0001),
       aes(y=reorder(celltype, proportion, median), x=proportion)) +
  coord_flip() +
  geom_boxplot(alpha=0.1) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = 8) +
  theme(axis.title.x = element_blank(),
        # angled text
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),) +
  ggtitle("Cell type proportions across tissue (prop > 0.0001)")
p_boxplot_all
ggsave(paste0(plot_path, "boxplot_all_", bin_size_str, ".png"), p_boxplot_all,
       width = 15, height = 8, bg = "white")


# DOUBLET MODE
doublet_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                   "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                   bin_size_str, ext),
                            header = TRUE)

doublet_info <- read.table(paste0(proportions_path, "_", bin_size_str,
                                  "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                  bin_size_str, ext, "_doublet_info.tsv"),
                           header = TRUE)

# Add rownames to deconv_props
rownames(doublet_props) <- colnames(visium_obj_subset)
all(rownames(doublet_props) == doublet_info$spot) # Check

table(doublet_info$spot_class)

# List colors from RColorBrewer paired palette
# RColorBrewer::brewer.pal(n = 6, name = "Paired")

classes_colors <- c("singlet"="#377EB8", "doublet_certain"="#33A02C",
                    "doublet_uncertain"="#B2DF8A", "reject"="#E41A1C")
spotclass_order <- names(classes_colors)

visium_obj_subset$spot_class <- factor(doublet_info$spot_class, levels = spotclass_order)
p_spot_class <- SpatialDimPlot(visium_obj_subset, group.by = "spot_class",
                               image.alpha = 0, stroke = NA) +
  scale_fill_manual(values = classes_colors) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
p_spot_class
ggsave(paste0(plot_path, "spatialdimplot_spotclass_", bin_size_str, ".png"), p_spot_class,
       width = 10, height = 12, bg = "white")

# Let's just plot the ROIs
doublet_props_roi <- doublet_props[rownames(doublet_props) %in% cells_roi, ]
doublet_info_roi <- doublet_info[doublet_info$spot %in% cells_roi, ]

all(colnames(visium_obj_roi) == doublet_info_roi$spot)
visium_obj_roi$spot_class <- factor(doublet_info_roi$spot_class, levels = spotclass_order)

p_spot_class_roi <- SpatialDimPlot(visium_obj_roi, group.by = "spot_class",
               image.alpha = 0, pt.size.factor = 14, shape=22, stroke=0.1) +
  scale_fill_manual(values = classes_colors) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right")
p_spot_class_roi

ggsave(paste0(plot_path, "spatialdimplot_spotclass_ROI", roi_name, 
              "_", bin_size_str, ".png"), p_spot_class_roi,
       width = 8, height = 6, bg = "white")

# Bar plot of spot classes over the whole tissue
ggplot(doublet_info, aes(y=spot_class)) +
  geom_bar(aes(fill=spot_class), width = 0.6, position="dodge") +
  scale_fill_manual(values = classes_colors) +
  scale_y_discrete(limits = rev(spotclass_order)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic(base_size = 8) +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  ggtitle(paste0("Spot classes across ", bin_size_str, " tissue"))
ggsave(paste0(plot_path, "barplot_spotclasses_", bin_size_str, ".png"),
       width = 5, height = 5, bg = "white")
