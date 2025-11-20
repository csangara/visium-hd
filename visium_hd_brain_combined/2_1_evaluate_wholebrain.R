library(tidyverse)
library(Seurat)
celltype_colors_df <- read.csv("data/ABA_metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_color.csv") %>% 
  inner_join(read.csv("data/ABA_metadata/WMB-taxonomy/20231215/views/cluster_to_cluster_annotation_membership_pivoted.csv"),
             by = "cluster_alias") %>% 
  mutate(class = str_replace_all(class, "[- ]", "")) %>% 
  distinct(class, class_color)
celltype_colors <-  celltype_colors_df$class_color %>% setNames(celltype_colors_df$class)
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

ggplot(merfish_cells_section %>% filter(brain_section_label == sections_oi[2]),
       aes(y=y_section, x=x_section, color=class)) +
  geom_point(size=0.1) + scale_y_reverse() + coord_fixed(ratio=1) +
  #geom_label(data=region_annot_label, aes(x=median_x, y=median_y, label=acronym_lvl5), size=5) +
  scale_color_manual(values = celltype_colors) +
  coord_fixed() +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_blank())


# Create barplot of cell types
merfish_cells_summ <- merfish_cells_section %>%
  rename(section = brain_section_label) %>% 
  group_by(class, class_color, section) %>%
  summarise(count = n()) %>%
  group_by(section) %>% 
  mutate(prop = count / sum(count))

ggplot(merfish_cells_summ,
       aes(x = section, y = prop, fill = class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  facet_wrap(~class) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MERFISH Cell Type Composition by Brain Region",
       x = "Brain Region",
       y = "Proportion of Cells",
       fill = "Cell Class")

first_run <- TRUE

# Read RCTD files
deconv_props_list <- list()
visium_objs_list <- list()
for (ds in c("ffpe", "fresh_frozen")){
  ext <- ifelse(ds == "ffpe", "_converted_doublet", "")
  ds_ext <- ifelse(ds == "ffpe", "", "_FF")
  
  visium_obj <- readRDS(paste0("data/Visium_HD_MouseBrain", ds_ext, "/Visium_HD_MouseBrain", ds_ext, "_008um.rds"))
  deconv_props <- read.table(paste0("visium_hd_brain", tolower(ds_ext), "/Visium_HD_MouseBrain", ds_ext,
                                    "_008um/proportions_rctd_Visium_HD_MouseBrain", ds_ext, "_008um", ext),
                             header = TRUE)
  removed_rows <- scan(paste0("visium_hd_brain", tolower(ds_ext), "/Visium_HD_MouseBrain", ds_ext,
                              "_008um/proportions_rctd_rows_removed"), what="character") %>% 
    .[grepl("s_008um", .)]
  
  stopifnot(length(removed_rows) + dim(deconv_props)[1] == dim(visium_obj)[2])
  visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_rows)]
  rownames(deconv_props) <- colnames(visium_obj_subset)
  
  visium_objs_list[[ds]] <- visium_obj_subset
  deconv_props_list[[ds]] <- deconv_props
  
  rm(visium_obj, visium_obj_subset, deconv_props, removed_rows)
}

deconv_props_df <- lapply(1:2, function(i){
  deconv_props_list[[i]] %>% tibble::rownames_to_column("spot") %>% 
    pivot_longer(-spot, names_to = "celltype", values_to = "proportion") %>% 
    mutate(celltype = str_replace_all(celltype, "^X", ""),
           source = names(deconv_props_list)[i])
}) %>% bind_rows()
  

# Summarise per region
deconv_props_summ <- deconv_props_df %>% group_by(celltype, source) %>%
  summarise(proportion = mean(proportion))

merfish_both_summ <- merfish_cells_summ %>% group_by(class) %>% 
  summarise(proportion = mean(prop)) %>% 
  mutate(source = "MERFISH_two_sections")

# Combined barplot
combined_summ_df <- bind_rows(deconv_props_summ,
                              merfish_cells_summ %>% select(class, prop, section) %>%
                                rename(celltype = class,
                                       proportion = prop,
                                       source = section),
                              merfish_both_summ %>% rename(celltype = class))

# Stacked barplot per region
ggplot(combined_summ_df,
       aes(x = source, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", width=0.6) +
  theme_minimal(base_size = 8) +
  facet_wrap(~celltype) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none")


first_run <- FALSE
deconv_props_df_shapes <- lapply(unique(deconv_props_summ$source), function(src) {
  src <- "fresh_frozen"
  square_size <- visium_objs_list[[src]]@images[["slice1.008um"]]@scale.factors$spot
  
  if (first_run){
    deconv_props_roi <- deconv_props_df %>%
      filter(source == src) %>%
      filter(proportion > 0) %>% 
      # merge with coordinates
      left_join(GetTissueCoordinates(visium_objs_list[[src]]), by = c("spot" = "cell")) %>% 
      mutate(n_celltypes = n(), .by = "spot",
             square_size = square_size)
    
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
    
    saveRDS(deconv_props_df_shapes, paste0("visium_hd_brain_combined/rds/deconv_props_shapes_",
                                             tolower(src), "_008um.rds"))
    
  } else {
    deconv_props_df_shapes <- readRDS(paste0("visium_hd_brain_combined/rds/deconv_props_shapes_",
                                             tolower(src), "_008um.rds"))
  }
  deconv_props_df_shapes$source <- src
  deconv_props_df_shapes
}) %>% bind_rows()

p_scatterbar <- ggplot(deconv_props_df_shapes %>% filter(source == "fresh_frozen"),
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

p_scatterplots[[1]]


# TODO: download and get proportions from Zhuang (X.Z) 1 and 2
