library(tidyverse)
celltype_colors_df <- read.csv("data/ABA_metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_color.csv") %>% 
  inner_join(read.csv("data/ABA_metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_pivoted.csv"),
             by = "cluster_alias") %>% 
  mutate(class = str_replace_all(class, "[- ]", "")) %>% 
  distinct(class, class_color)
celltype_colors <-  celltype_colors_df$class_color %>% setNames(celltype_colors_df$class)

parcellations <- inner_join(read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_color.csv"),
                            read.csv("data/ABA_metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv"),
                            by="parcellation_index")
division_colors_df <- parcellations %>% distinct(division, division_color) %>% 
  # If color is #FFFFFF or #CCCCCC, change division to "fiber tracts/VS"
  mutate(division = ifelse(division_color %in% c("#FFFFFF", "#CCCCCC"),
                           "fiber tracts/VS", division)) %>% 
  distinct(division, division_color)
division_colors <- setNames(division_colors_df$division_color, division_colors_df$division)


sections_oi <- c("C57BL6J-638850.40","C57BL6J-638850.50")




merfish_cells <- read.csv("data/ABA_metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/cell_metadata_with_parcellation_annotation.csv")
merfish_cells %>% head

merfish_cells$brain_section_label %>% table

merfish_cells_section <- merfish_cells %>% 
  filter(brain_section_label %in% sections_oi) %>% 
  mutate(class = str_replace_all(class, "[- ]", ""))

ggplot(merfish_cells_section %>% filter(brain_section_label == sections_oi[1]),
       aes(y=y_section, x=x_section, color=class)) +
  geom_point(size=0.1) + scale_y_reverse() + coord_fixed(ratio=1) +
  #geom_label(data=region_annot_label, aes(x=median_x, y=median_y, label=acronym_lvl5), size=5) +
  scale_color_manual(values = celltype_colors) +
  coord_fixed() +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_blank())

# Get the region annotations I made
region_annotations <- readRDS("data/Visium_HD_MouseBrain/tissue_positions_with_annotations_008um.rds")

region_annot_label <- region_annotations %>% 
  group_by(division) %>% 
  summarise(median_x = median(pxl_col_in_lowres),
            median_y = median(pxl_row_in_lowres))

p_regions <- ggplot(region_annotations, aes(y=pxl_row_in_lowres, x=pxl_col_in_lowres, fill=division)) +
  geom_bin2d(bins=500) + scale_y_reverse() + coord_fixed(ratio=1) +
  geom_label(data=region_annot_label, aes(x=median_x, y=median_y, label=division), size=5) +
  coord_fixed() +
  theme_classic() +
  # scale_fill_manual(values = division_colors) +
  theme(legend.position = "none",
        axis.title = element_blank())
p_regions
ggsave("visium_hd_brain_combined/plots/annotated_regions_ffpe.png", p_regions,
       width = 8, height = 6, bg = "white")

# Filter to only spots that were deconvolved
region_annotations <- region_annotations %>% filter(barcode %in% deconv_props_df %>% filter(source == "ffpe") %>% pull(spot) %>% unique,
                                                    !is.na(division))

region_annotations$division %>% table

common_division <- lapply(5:8, function(levels) {
  intersect(unique(merfish_cells$parcellation_division), 
            unique(region_annotations[[paste0("acronym_lvl", levels)]]))
}) %>% unlist

common_structure <- lapply(5:8, function(levels) {
  intersect(unique(merfish_cells$parcellation_structure), 
            unique(region_annotations[[paste0("acronym_lvl", levels)]]))
}) %>% unlist

common_substructure <- lapply(5:8, function(levels) {
  intersect(unique(merfish_cells$parcellation_substructure), 
            unique(region_annotations[[paste0("acronym_lvl", levels)]]))
}) %>% unlist

# Cells with common division (2369220 out of 3739961)
merfish_cells %>% filter(parcellation_division %in% common_division) %>% nrow
# Cells with common structure (909472)
merfish_cells %>% filter(parcellation_structure %in% common_structure) %>% nrow
# Cells with common substructure (485184)
merfish_cells %>% filter(parcellation_substructure %in% common_substructure) %>% nrow

# Check intersection
 intersect(merfish_cells %>% filter(parcellation_structure %in% common_structure) %>% pull(cell_label),
           merfish_cells %>% filter(parcellation_substructure %in% common_substructure) %>% pull(cell_label)) %>% 
   length
# = 485184
 
# Let's just use division level for now - but only with acronym lvl 5 (2321423)
merfish_cells %>% filter(parcellation_division %in% unique(region_annotations$division)) %>% nrow
 
merfish_cells_filtered <- merfish_cells %>% 
  filter(parcellation_division %in% unique(region_annotations$division)) %>%
  rename(region = parcellation_division)

# Create barplot per region
merfish_cells_summ <- merfish_cells_filtered %>%
  group_by(region, class, class_color) %>%
  summarise(count = n()) %>%
  group_by(region) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>% 
  mutate(class = str_replace_all(class, "[- ]", "")) 

ggplot(merfish_cells_summ %>%  filter(region != "PAL"),
       aes(x = region, y = prop, fill = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MERFISH Cell Type Composition by Brain Region",
       x = "Brain Region",
       y = "Proportion of Cells",
       fill = "Cell Class")

# Read RCTD
visium_obj <- readRDS("data/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds")

ext <- "_converted_doublet" # "", "_converted_doublet", "_converted_full"
deconv_props <- read.table(paste0("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um", ext),
                           header = TRUE)
dim(deconv_props) #278323 spots

# Get removed rows from the other file
removed_rows <- scan("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_rows_removed", what="character") %>% 
  .[grepl("s_008um", .)]

length(removed_rows) # 115220

# Check if removed rows + leftover rows == total rows (yes)
length(removed_rows) + dim(deconv_props)[1] == dim(visium_obj)[2]

# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_rows)]

dim(visium_obj_subset) # 19059 genes x 278323 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

deconv_props_df <- deconv_props %>% tibble::rownames_to_column("barcode") %>% 
  filter(barcode %in% region_annotations$barcode) %>% 
  pivot_longer(-barcode, names_to = "celltype", values_to = "proportion") %>% 
  inner_join(region_annotations,
             by = "barcode")


deconv_props_df$acronym_lvl5 %>% table

# Summarise per region
deconv_props_summ <- deconv_props_df %>% group_by(acronym_lvl5, celltype) %>%
  summarise(proportion = mean(proportion)) %>% 
  # remove X in front of cell type names
  ungroup() %>% 
  mutate(celltype = str_replace_all(celltype, "^X", "")) %>% 
  # remove NA and PAL
  filter(acronym_lvl5 != "PAL") %>% 
  filter(!is.na(acronym_lvl5))

# Stacked barplot per region
ggplot(deconv_props_summ, aes(x = acronym_lvl5, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width=0.6) +
  theme_minimal(base_size = 8) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  # facet_grid(~region_broad, scales = "free_x", space='free') +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank())

# Combined barplot
combined_summ_df <- bind_rows(deconv_props_summ %>% rename(region = acronym_lvl5) %>%
            mutate(source = "VisiumHD RCTD"),
          merfish_cells_summ %>% select(region, class, prop) %>%
            rename(celltype = class,
                   proportion = prop) %>%
            mutate(source = "MERFISH")) %>%
  filter(region != "PAL")

ggplot(combined_summ_df,
       aes(x = source, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width=0.6) +
  theme_minimal(base_size = 8) +
  facet_wrap(~region) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(panel.grid.major.x = element_blank())


# Calculate stats between the two
# Do spatial scatterplot of the cortex
square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

deconv_props_roi <- deconv_props %>%
  rownames_to_column("spot") %>%
  pivot_longer(cols = -spot, names_to = "celltype",
               values_to = "proportion") %>% 
  filter(proportion > 0) %>% 
  # merge with coordinates
  left_join(GetTissueCoordinates(visium_obj_subset), by = c("spot" = "cell")) %>% 
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
  # Remove X from celltype 
  mutate(celltype = gsub("^X", "", celltype)) %>% 
  mutate(celltype = factor(celltype, levels = names(celltype_colors)))

p_scatterbar <- ggplot(deconv_props_df_shapes,
                       aes(x = x, y = y)) +
  geom_polygon(aes(fill = celltype, group = group), show.legend = TRUE) +
  # White border
  # geom_tile(data = deconv_props_roi %>% distinct(x, y),
  #           aes(x = y, y = x), height = square_size, width = square_size,
  #           fill = NA, color = "white", inherit.aes = FALSE) +
  scale_fill_manual(values = celltype_colors) +
  theme_void(base_size=8) +
  scale_y_reverse() +
  coord_fixed() +
  guides(fill = guide_legend(ncol=1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"))
p_scatterbar

## Visium FF

visium_obj <- readRDS("data/Visium_HD_MouseBrain_FF/Visium_HD_MouseBrain_FF_008um.rds")
ext <- "" # "", "_converted_doublet", "_converted_full"
deconv_props <- read.table(paste0("visium_hd_brain_ff/Visium_HD_MouseBrain_FF_008um/proportions_rctd_Visium_HD_MouseBrain_FF_008um", ext),
                           header = TRUE)
removed_rows <- scan("visium_hd_brain_ff/Visium_HD_MouseBrain_FF_008um/proportions_rctd_rows_removed", what="character") %>% 
  .[grepl("s_008um", .)]
# Check if removed rows + leftover rows == total rows (yes)
length(removed_rows) + dim(deconv_props)[1] == dim(visium_obj)[2]
# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_rows)]
# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

# TODO: download and get proportions from Zhuang (X.Z) 1 and 2
