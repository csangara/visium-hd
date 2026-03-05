library(Seurat)
library(tidyverse)

celltype_colors_df <- read.csv("data/scref_MouseBrain_ABA/cell_metadata_with_cluster_annotation_downsampled.csv") %>%
  distinct(class, class_color) %>% 
  mutate(class = str_replace_all(class, "[- ]", "")) %>% 
  mutate(class = paste0("X", class)) %>% arrange(class)
celltype_colors <- celltype_colors_df$class_color %>% setNames(celltype_colors_df$class)  

visium_obj <- readRDS("data/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds")
dim(visium_obj) # 19059 genes x 393543 spots

ext <- "_full" # "" or "_full"

deconv_props <- read.table(paste0("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um", ext),
                           header = TRUE)
dim(deconv_props) #278323 spots

removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,"nCount_Spatial.008um"] < 100)]

# Check if removed rows + leftover rows == total rows (yes)
length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2]

# Print number of removed spots, total spots, and percentage
cat("Number of removed spots:", length(removed_spots), "\n")
cat("Total spots:", dim(visium_obj)[2], "\n")
cat("Percentage of removed spots:", 
    round(length(removed_spots) / dim(visium_obj)[2] * 100, 2), "%\n")

# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]

dim(visium_obj_subset) # 19059 genes x 278323 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

p_ncount <- SpatialFeaturePlot(visium_obj, "nCount_Spatial.008um") +
  labs(fill = "nCount") +
  theme(legend.position = "right")
ggsave("visium_hd_brain/plots/spatialfeatureplot_nCount.png", p_ncount,
       width = 8, height = 6, bg = "white") 

p_ncount_subset <- SpatialFeaturePlot(visium_obj_subset, "nCount_Spatial.008um") +
  labs(fill = "nCount") +
  theme(legend.position = "right")
ggsave("visium_hd_brain/plots/spatialfeatureplot_nCount_subset.png", p_ncount_subset,
       width = 8, height = 6, bg = "white") 

# Assign barcode to most abundant cell type per spot
visium_obj_subset$celltype <- factor(colnames(deconv_props)[max.col(deconv_props)],
                                     levels = names(celltype_colors))

p_celltype <- SpatialDimPlot(visium_obj_subset, group.by = "celltype",
               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = celltype_colors,
                    labels = function(x) gsub("^X", "", x),
                    name="Cell type") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right",
        legend.text = element_text(size=5),
        legend.title = element_text(size=6),
        legend.key.size = unit(0.3, "cm")
        )
p_celltype

# ggsave(paste0("visium_hd_brain/plots/spatialdimplot_celltype", ext, ".pdf"), p_celltype,
#        width = 8, height = 6, bg = "white")
ggsave(paste0("visium_hd_brain/plots/spatialdimplot_celltype", ext, ".png"), p_celltype,
       width = 8, height = 6, bg = "white", dpi=300)

# Get region annotations
region_annotations <- readRDS("data/Visium_HD_MouseBrain/tissue_positions_with_annotations_008um.rds")

# Filter to only spots that were deconvolved
region_annotations <- region_annotations %>% filter(barcode %in% rownames(deconv_props),
                                                    !grepl("unassigned", division))

# Label regions in the plot, calculate median x and y coordinates of each region
region_annot_label <- region_annotations %>% 
  group_by(division) %>% 
  summarise(median_x = median(pxl_col_in_lowres),
            median_y = median(pxl_row_in_lowres))

p_regions <- ggplot(region_annotations,aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres, fill=division)) +
  geom_bin2d(bins=500) + scale_y_reverse() + coord_fixed(ratio=1) +
  geom_label(data=region_annot_label, aes(x=median_x, y=median_y, label=division), size=5) +
  coord_fixed() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank())
p_regions
ggsave("visium_hd_brain/plots/annotated_regions.png", p_regions,
       width = 8, height = 6, bg = "white")

# Filter deconvolved proportions to match region_annotations
deconv_props_df <- deconv_props %>% tibble::rownames_to_column("barcode") %>% 
  filter(barcode %in% region_annotations$barcode) %>% 
  pivot_longer(-barcode, names_to = "celltype", values_to = "proportion") %>% 
  inner_join(region_annotations %>% select(barcode, division, pxl_col_in_lowres, pxl_row_in_lowres),
             by = "barcode")

# Summarise per region
deconv_props_summ <- deconv_props_df %>% group_by(division, celltype) %>%
  summarise(proportion = mean(proportion))

# Stacked barplot per region
p <- ggplot(deconv_props_summ, aes(x = division, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width=0.6) +
  theme_minimal(base_size = 8) +
  scale_fill_manual(values = celltype_colors,
                    labels = function(x) gsub("^X", "", x),
                    name="Cell type") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank())
p
# ggsave(paste0("visium_hd_brain/plots/barplot_avg_by_region", ext, ".pdf"), p,
#        width = 8, height = 6, bg = "white")
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
p_boxplot_all
ggsave(paste0("visium_hd_brain/plots/boxplot_all", ext, ".png"), p_boxplot_all,
       width = 15, height = 8, bg = "white")


#### Only for Doublet mode
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
  inner_join(region_annotations %>% select(barcode, division),
             by = c("spot" = "barcode"))

p_region_spotclass <- ggplot(doublet_info_region, aes(x = division,
                              fill = factor(spot_class, levels=names(classes_colors)))) +
  geom_bar(position = "fill", width=0.6) +
  theme_minimal(base_size = 8) +
  scale_fill_manual(values = classes_colors, name="Spot class") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank())
# p_region_spotclass
ggsave(paste0("visium_hd_brain/plots/barplot_spot_class_by_region", ext, ".png"),
       p_region_spotclass,
       width = 8, height = 6, bg = "white")