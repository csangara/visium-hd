library(Seurat)
library(tidyverse)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

visium_obj <- readRDS("data/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds")
dim(visium_obj) # 19059 genes x 393543 spots

ext <- "" #"doublet_hd" or ""

deconv_props <- read.table(paste0("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um", ext),
                           header = TRUE)
dim(deconv_props) #278323 spots

# Get removed rows from the other file
removed_rows <- scan("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_rows_removed", what="character") %>% 
  .[grepl("s_008um", .)]

length(removed_rows) # 115220

# Check if removed rows + leftover rows == total rows (yes)
length(removed_rows) + dim(deconv_props)[1] == dim(visium_obj)[2]

# Subset visium_obj to only include spots that were not removedqual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_rows)]

dim(visium_obj_subset) # 19059 genes x 278323 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

# Visualize spots in visium_obj
# Read parquet file
tissue.positions <- arrow::read_parquet("data/Visium_HD_MouseBrain/binned_outputs/square_008um/spatial/tissue_positions.parquet")

# Check that all "barcode" of tissue.positions are in visium_obj and visium_obj_subset
all(colnames(visium_obj) %in% tissue.positions$barcode)
all(colnames(visium_obj_subset) %in% tissue.positions$barcode)

visium_obj$nCount_RNA <- GetAssayData(visium_obj, layer = "counts") %>% colSums()
visium_obj_subset$nCount_RNA <- GetAssayData(visium_obj_subset, layer = "counts") %>% colSums()

# Check counts of removed rows
visium_obj[, colnames(visium_obj) %in% removed_rows]$nCount_RNA %>% hist()

# Visualize which spots were removed
visium_obj_meta <- visium_obj@meta.data %>% 
  # Add x and y coordinates
  left_join(tissue.positions, by = c("name" = "barcode"))

# Plot spots
p_visium_bin2d <- ggplot(visium_obj_meta, aes(y = pxl_row_in_fullres, x = pxl_col_in_fullres, fill=nCount_RNA, group=nCount_RNA)) +
  geom_bin2d(bins=500) +
  coord_fixed() +
  scale_y_reverse() +
  theme_minimal(base_size = 8)
#p_visium_bin2d
ggsave("visium_hd_brain/plots/visium_bin2d.png", p_visium_bin2d,
       width = 8, height = 6, bg = "white")

# Visum subset
visium_obj_subset_meta <- visium_obj_subset@meta.data %>% 
  # Add x and y coordinates
  left_join(tissue.positions, by = c("name" = "barcode"))

# Plot spots
p_visium_subset_bin2d <- ggplot(visium_obj_subset_meta, aes(y = pxl_row_in_fullres, x = pxl_col_in_fullres, fill=nCount_RNA, group=nCount_RNA)) +
  geom_bin2d(bins=500) +
  coord_fixed() +
  scale_y_reverse() +
  theme_minimal(base_size = 8)
# p_visium_subset_bin2d
ggsave("visium_hd_brain/plots/visium_subset_bin2d.png", p_visium_subset_bin2d,
       width = 8, height = 6, bg = "white")


deconv_props_df <- deconv_props %>% pivot_longer(everything(), names_to = "celltype", values_to = "proportion") %>% 
  group_by(celltype) %>% summarise(proportion = mean(proportion))

# Plot proportions as barplot
ggplot(deconv_props_df, aes(y = celltype, x = proportion)) +
  geom_bar(stat = "identity")

# Assign barcode to most abundant cell type per spot
celltype_position <- data.frame(barcode = rownames(deconv_props),
           celltype = colnames(deconv_props)[max.col(deconv_props)]) %>% 
  left_join(tissue.positions, by = c("barcode" = "barcode"))

# Plot spots colored by most abundant cell type
p_celltype <- ggplot(celltype_position, aes(y = pxl_row_in_fullres, x = pxl_col_in_fullres, color=celltype)) +
  #geom_bin2d(bins=500) +
  geom_point(size=0.1) +
  coord_fixed() +
  scale_y_reverse() +
  # Increase size of legend
  scale_color_manual(values = col_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal(base_size = 8)

# Save plot
ggsave(paste0("visium_hd_brain/plots/spot_celltype", ext, ".png"), p_celltype,
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

# Filter deconvolved proportions to match region_annotations
deconv_props_df <- deconv_props %>% tibble::rownames_to_column("barcode") %>% 
  filter(barcode %in% region_annotations$barcode) %>% 
  pivot_longer(-barcode, names_to = "celltype", values_to = "proportion") %>% 
  inner_join(region_annotations %>% select(barcode, acronym_lvl6, pxl_col_in_lowres, pxl_row_in_lowres, region_broad),
             by = "barcode")

# Summarise per region
deconv_props_summ <- deconv_props_df %>% group_by(region_broad, acronym_lvl6, celltype) %>% summarise(proportion = mean(proportion))

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
ggsave(paste0("visium_hd_brain/plots/rctd_avg_by_region", ext, ".png"), p,
       width = 8, height = 6, bg = "white")


# Label regions in the plot, calculate median x and y coordinates of each region
region_annot_label <- region_annotations %>% 
  group_by(acronym_lvl6) %>% 
  summarise(median_x = median(pxl_col_in_lowres), median_y = median(pxl_row_in_lowres))

p_regions <- ggplot(region_annotations,aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres, fill=acronym_lvl6)) +
  geom_bin2d(bins=500) + scale_y_reverse() + coord_fixed(ratio=1) +
  geom_label(data=region_annot_label, aes(x=median_x, y=median_y, label=acronym_lvl6), size=5) +
  coord_fixed() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank())
ggsave("visium_hd_brain/plots/annotated_regions.png", p_regions,
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

# In the hippocampal region, what is the distribution of cell types?
p_boxplot_hip <- ggplot(deconv_props_df %>% filter(acronym_lvl6 == "HIP", proportion > 0.0001),
       aes(y=celltype, x=proportion)) +
  geom_boxplot(alpha=0.1) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  ggtitle("Hippocampal cell type proportions (prop > 0.0001)")
ggsave(paste0("visium_hd_brain/plots/boxplot_hip", ext, ".png"), p_boxplot_hip,
       width = 8, height = 10, bg = "white")

# DOUBLET MODE
doublet_props <- read.table("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um_doublet_hd",
                            header = TRUE)
doublet_info <- read.table("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um_doublet_hd_doublet_info.tsv",
                           header = TRUE)

# Add rownames to deconv_props
rownames(doublet_props) <- colnames(visium_obj_subset)
all(rownames(doublet_props) == doublet_info$spot) # Check

table(doublet_info$spot_class)

doublet_info <- doublet_info %>% 
  left_join(tissue.positions, by = c("spot" = "barcode")) 
classes_colors <- c("singlet"="forestgreen", "doublet_certain"="navyblue",
                    "doublet_uncertain"="orange", "reject"="red")
# Plot spots by prediction type
p_spot_class <- ggplot(doublet_info, aes(y = pxl_row_in_fullres, x = pxl_col_in_fullres,
                                         color=factor(spot_class, levels=names(classes_colors)))) +
  #geom_bin2d(bins=500) +
  geom_point(size=0.1) +
  coord_fixed() +
  scale_y_reverse() +
  # Increase size of legend
  scale_color_manual(values = classes_colors) +
  guides(color = guide_legend("Spot class", override.aes = list(size = 3))) +
  theme_minimal(base_size = 8)
# p_spot_class
ggsave("visium_hd_brain/plots/spot_class.png", p_spot_class,
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
ggsave("visium_hd_brain/plots/spot_class_by_region.png", p_region_spotclass,
       width = 8, height = 6, bg = "white")

# Assign barcode to most abundant cell type per spot
celltype_position_doublet <- data.frame(barcode = rownames(doublet_props),
                                celltype = colnames(doublet_props)[max.col(doublet_props)]) %>% 
  left_join(tissue.positions, by = c("barcode" = "barcode"))

# Make confusion matrix
celltypes <- celltype_position$celltype %>% unique %>% sort
conf_matrix <- caret::confusionMatrix(celltype_position$celltype %>% factor(levels=celltypes),
                       celltype_position_doublet$celltype %>% factor(levels=celltypes))
pheatmap::pheatmap(conf_matrix$table, scale="row",
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100),
                   cluster_rows = F, cluster_cols = FALSE,
                   main = "Full mode (rows) vs Doublet mode (cols) predictions",
                   filename = "visium_hd_brain/plots/conf_matrix_full_vs_doublet.png")

# Barplot 
doublet_props_df <- doublet_props %>% tibble::rownames_to_column("barcode") %>% 
  filter(barcode %in% region_annotations$barcode) %>% 
  pivot_longer(-barcode, names_to = "celltype", values_to = "proportion") %>% 
  inner_join(region_annotations %>% select(barcode, acronym_lvl6, pxl_col_in_lowres, pxl_row_in_lowres, region_broad),
             by = "barcode")

# Summarise per region
doublet_props_summ <- doublet_props_df %>% group_by(region_broad, acronym_lvl6, celltype) %>% summarise(proportion = mean(proportion))

# Get correlation of proportions between full and doublet mode
mean(diag(cor(deconv_props, doublet_props)))

# Overlap of most abundant cell type
max_celltypes <- data.frame(celltype_full = colnames(deconv_props)[max.col(deconv_props)],
                                celltype_doublet = colnames(doublet_props)[max.col(doublet_props)])
max_celltypes %>% filter(celltype_full == celltype_doublet) %>% nrow / nrow(max_celltypes)