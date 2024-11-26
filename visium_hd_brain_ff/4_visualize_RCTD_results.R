library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(grid)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

rotate_image <- function(p, rot_angle) {
  gt <- ggplot_gtable(ggplot_build(p))
  panel_idx <- which(gt$layout$name == "panel")
  rot_vp <- viewport(angle = rot_angle)
  gt[["grobs"]][[panel_idx]] <- editGrob(gt[["grobs"]][[panel_idx]], vp = rot_vp)
  p_rot <- ggdraw() + draw_grob(gt)
  
  return(p_rot)
}
##########################################

visium_obj <- readRDS("data/Visium_HD_MouseBrain_FF/Visium_HD_MouseBrain_FF_008um.rds")
dim(visium_obj) # 19059 genes x 453820 spots

ext <- "" # "", "_converted_doublet", "_converted_full"

deconv_props <- read.table(paste0("visium_hd_brain_ff/Visium_HD_MouseBrain_FF_008um/proportions_rctd_Visium_HD_MouseBrain_FF_008um", ext),
                           header = TRUE)
dim(deconv_props) #451697 spots

# Get removed rows from the other file
removed_rows <- scan("visium_hd_brain_ff/Visium_HD_MouseBrain_FF_008um/proportions_rctd_rows_removed", what="character") %>% 
  .[grepl("s_008um", .)]

length(removed_rows) # 2123

# Check if removed rows + leftover rows == total rows (yes)
length(removed_rows) + dim(deconv_props)[1] == dim(visium_obj)[2]

# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_rows)]

dim(visium_obj_subset) # 19059 genes x 451697 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

# Check counts of removed rows
visium_obj[, colnames(visium_obj) %in% removed_rows]$nCount_Spatial.008um %>% hist()

p_ncount <- SpatialFeaturePlot(visium_obj, "nCount_Spatial.008um", image.alpha=0) +
  labs(fill = "nCount") +
  theme(legend.position = "right")
ggsave("visium_hd_brain_ff/plots/spatialfeatureplot_nCount.png",
       rotate_image(p_ncount, -90),
       width = 8, height = 6, bg = "white") 

orig_xlim <- ggplot_build(p_ncount)$layout$panel_scales_x[[1]]$range$range
orig_ylim <- ggplot_build(p_ncount)$layout$panel_scales_y[[1]]$range$range

p_ncount_subset <- SpatialFeaturePlot(visium_obj_subset, "nCount_Spatial.008um", image.alpha=0) +
  xlim(orig_xlim) + ylim(orig_ylim) +
  labs(fill = "nCount") + 
  theme(legend.position = "right") 
ggsave("visium_hd_brain_ff/plots/spatialfeatureplot_nCount_subset.png",
       rotate_image(p_ncount_subset, -90),
       width = 8, height = 6, bg = "white") 

# Assign barcode to most abundant cell type per spot
all(rownames(deconv_props) == colnames(visium_obj_subset))
visium_obj_subset$celltype <- colnames(deconv_props)[max.col(deconv_props)]

p_celltype <- SpatialDimPlot(visium_obj_subset, group.by = "celltype",
               image.alpha = 0, stroke=NA) +
  xlim(orig_xlim) + ylim(orig_ylim) +
  scale_fill_manual(values = col_vector) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right")
ggsave(paste0("visium_hd_brain_ff/plots/spatialdimplot_celltype", ext, ".png"),
       rotate_image(p_celltype, -90),
       width = 8, height = 6, bg = "white")

# Get region annotations
region_annotations <- readRDS("data/Visium_HD_MouseBrain/tissue_positions_with_annotations_008um.rds")

# Filter to only spots that were deconvolved
region_annotations <- region_annotations %>% filter(barcode %in% rownames(deconv_props),
                                                    !is.na(acronym_lvl6)) %>% 
  mutate(region_broad = case_when(
    acronym_lvl6 %in% c("VIS", "PTLp", "SS", "AUD", "TEa", "RHP") ~ "Cerebral cortex",
    acronym_lvl6 %in% c("STRd", "sAMY") ~ "Cerebral nuclei",
    acronym_lvl6 %in% c("DORpm", "DORsm") ~ "Thalamus",
    acronym_lvl6 %in% c("MEZ", "LZ", "PVR", "PVZ") ~ "Hypothalamus",
    acronym_lvl6 %in% c("HIP") ~ "Hippocampus",
    TRUE ~ acronym_lvl6
  ))

# Label regions in the plot, calculate median x and y coordinates of each region
region_annot_label <- region_annotations %>% 
  group_by(acronym_lvl6) %>% 
  summarise(median_x = median(pxl_col_in_lowres),
            median_y = median(pxl_row_in_lowres))

p_regions <- ggplot(region_annotations,aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres, fill=acronym_lvl6)) +
  geom_bin2d(bins=500) + scale_y_reverse() + coord_fixed(ratio=1) +
  geom_label(data=region_annot_label, aes(x=median_x, y=median_y, label=acronym_lvl6), size=5) +
  coord_fixed() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank())
ggsave("visium_hd_brain/plots/annotated_regions.png", p_regions,
       width = 8, height = 6, bg = "white")

# Filter deconvolved proportions to match region_annotations
deconv_props_df <- deconv_props %>% tibble::rownames_to_column("barcode") %>% 
  filter(barcode %in% region_annotations$barcode) %>% 
  pivot_longer(-barcode, names_to = "celltype", values_to = "proportion") %>% 
  inner_join(region_annotations %>% select(barcode, acronym_lvl6, pxl_col_in_lowres, pxl_row_in_lowres, region_broad),
             by = "barcode")

# Summarise per region
deconv_props_summ <- deconv_props_df %>% group_by(region_broad, acronym_lvl6, celltype) %>%
  summarise(proportion = mean(proportion))

# Stacked barplot per region
p <- ggplot(deconv_props_summ, aes(x = acronym_lvl6, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width=0.6) +
  theme_minimal(base_size = 8) +
  scale_fill_manual(values = col_vector) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  facet_grid(~region_broad, scales = "free_x", space='free') +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
p
ggsave(paste0("visium_hd_brain/plots/barplot_avg_by_region", ext, ".png"), p,
       width = 8, height = 6, bg = "white")


# What is the distribution of cell types across the tissue?
p_boxplot_all <- ggplot(deconv_props_df %>% filter(proportion > 0.0001),
       aes(y=celltype, x=proportion)) +
  coord_flip() +
  geom_boxplot(alpha=0.1) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = 8) +
  theme(axis.title.x = element_blank(),
        # angled text
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),) +
  ggtitle("Cell type proportions across tissue (prop > 0.0001)")
#p_boxplot_all
ggsave(paste0("visium_hd_brain/plots/boxplot_all", ext, ".png"), p_boxplot_all,
       width = 15, height = 8, bg = "white")


# DOUBLET MODE
ext <- "" # "", "_converted_doublet"
doublet_props <- read.table(paste0("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um", ext),
                            header = TRUE)
doublet_info <- read.table(paste0("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um", ext, "_doublet_info.tsv"),
                           header = TRUE)

# Add rownames to deconv_props
rownames(doublet_props) <- colnames(visium_obj_subset)
all(rownames(doublet_props) == doublet_info$spot) # Check

table(doublet_info$spot_class)

classes_colors <- c("singlet"="forestgreen", "doublet_certain"="navyblue",
                    "doublet_uncertain"="orange", "reject"="red")
visium_obj_subset$spot_class <- doublet_info$spot_class
p_spot_class <- SpatialDimPlot(visium_obj_subset, group.by = "spot_class",
               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = classes_colors) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right")
ggsave(paste0("visium_hd_brain/plots/spatialdimplot_spot_class", ext, ".png"), p_spot_class,
       width = 8, height = 6, bg = "white")

# Plot by region
doublet_info_region <- doublet_info %>%
  filter(spot %in% region_annotations$barcode) %>% 
  inner_join(region_annotations %>% select(barcode, acronym_lvl6, region_broad),
             by = c("spot" = "barcode"))

p_region_spotclass <- ggplot(doublet_info_region, aes(x = acronym_lvl6,
                              fill = factor(spot_class, levels=names(classes_colors)))) +
  geom_bar(position = "fill", width=0.6) +
  theme_minimal(base_size = 8) +
  scale_fill_manual(values = classes_colors, name="Spot class") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  facet_grid(~region_broad, scales = "free_x", space='free') +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank())
# p_region_spotclass
ggsave(paste0("visium_hd_brain/plots/barplot_spot_class_by_region", ext, ".png"),
       p_region_spotclass,
       width = 8, height = 6, bg = "white")

# Make confusion matrix
celltypes <- celltype_position$celltype %>% unique %>% sort
conf_matrix <- caret::confusionMatrix(colnames(deconv_props)[max.col(deconv_props)] %>% factor(levels=celltypes),
                                      colnames(doublet_props)[max.col(doublet_props)] %>% factor(levels=celltypes))
pheatmap::pheatmap(conf_matrix$table, scale="row",
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100),
                   cluster_rows = F, cluster_cols = FALSE,
                   main = "Full mode (rows) vs Doublet mode (cols) predictions",
                   filename = "visium_hd_brain/plots/conf_matrix_full_vs_doublet.png")

# Get correlation of proportions between full and doublet mode
mean(diag(cor(deconv_props, doublet_props)))

# Overlap of most abundant cell type
max_celltypes <- data.frame(celltype_full = colnames(deconv_props)[max.col(deconv_props)],
                                celltype_doublet = colnames(doublet_props)[max.col(doublet_props)])
max_celltypes %>% filter(celltype_full == celltype_doublet) %>% nrow / nrow(max_celltypes)
