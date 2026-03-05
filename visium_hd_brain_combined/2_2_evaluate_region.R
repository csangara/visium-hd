library(tidyverse)
library(Seurat)
source("visium_hd_brain_combined/0_brain_init.R")

# Read RCTD files 
deconv_props_list <- list()
visium_objs_list <- list()
for (ds in c("ffpe", "fresh_frozen")){
  ds_ext <- ifelse(ds == "ffpe", "", "_FF")
  
  visium_obj <- readRDS(paste0("data/Visium_HD_MouseBrain", ds_ext, "/Visium_HD_MouseBrain", ds_ext, "_008um.rds"))
  deconv_props <- read.table(paste0("visium_hd_brain", tolower(ds_ext), "/Visium_HD_MouseBrain", ds_ext,
                                    "_008um/proportions_rctd_Visium_HD_MouseBrain", ds_ext, "_008um"),
                             header = TRUE)
  
  removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,"nCount_Spatial.008um"] < 100)]
  
  stopifnot(length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2])
  visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]
  rownames(deconv_props) <- colnames(visium_obj_subset)
  
  visium_objs_list[[ds]] <- visium_obj_subset
  deconv_props_list[[ds]] <- deconv_props
  
  rm(visium_obj, visium_obj_subset, deconv_props, removed_spots)
}

deconv_props_df <- lapply(1:2, function(i){
  deconv_props_list[[i]] %>% tibble::rownames_to_column("spot") %>% 
    pivot_longer(-spot, names_to = "celltype", values_to = "proportion") %>% 
    mutate(celltype = str_replace_all(celltype, "^X", ""),
           source = names(deconv_props_list)[i])
}) %>% bind_rows()

# Get the region annotations I made for the Visium HD mouse brain FFPE data
region_annotations_df <- lapply(c("ffpe", "fresh_frozen"), function(ds){
  ds_ext <- ifelse(ds == "ffpe", "", "_FF")
  region_annotations <- readRDS(paste0("data/Visium_HD_MouseBrain", ds_ext, "/tissue_positions_with_annotations_008um.rds")) %>% 
    filter(!grepl("unassigned", division))
  
  # Filter to only spots that were deconvolved
  region_annotations <- region_annotations %>% filter(barcode %in% colnames(visium_objs_list[[ds]]))
  
  region_annotations %>% 
    mutate(source = ds)
  
}) %>% bind_rows()
                       
# Read merfish_cells from 1_3
merfish_cells <- readRDS("visium_hd_brain_combined/rds/merfish_cells_3datasets_combined.rds")

common_regions <- intersect(unique(region_annotations_df$division),
                            merfish_cells$division) %>% 
  # remove any that is 'unassigned'
  .[!grepl("unassigned", .)]

p_regions <- lapply(c("ffpe", "fresh_frozen"), function(ds){
  # Filter to only spots that were deconvolved
  ds_ext <- ifelse(ds == "ffpe", "", "_FF")
  region_annotations <- readRDS(paste0("data/Visium_HD_MouseBrain", ds_ext, "/tissue_positions_with_annotations_008um.rds")) %>% 
    filter(!grepl("unassigned", division))
  
  if (ds == "fresh_frozen"){
    # Swap x and y
    region_annotations_rename <- region_annotations %>% 
      rename(x_axis = pxl_row_in_lowres,
             y_axis = pxl_col_in_lowres) %>%
      mutate(x_axis = -x_axis, y_axis=-y_axis) %>% 
      filter(x_axis < -40)
  } else {
    region_annotations_rename <- region_annotations %>% 
      rename(x_axis = pxl_col_in_lowres,
             y_axis = pxl_row_in_lowres) %>% 
      mutate(y_axis = -y_axis) 
  }
  
  region_annot_label <- region_annotations_rename %>% 
    group_by(division) %>% 
    summarise(median_x = median(x_axis),
              median_y = median(y_axis)) %>% 
    # rename fiber tracts/VS to ""
    mutate(division = ifelse(division == "fiber tracts/VS", "", division))
    
  
  ggplot(region_annotations_rename, aes(x=x_axis, y=y_axis, fill=division)) +
    geom_bin2d(bins=100, linewidth=0,
               show.legend = ifelse(ds=="ffpe", TRUE, FALSE)) +
    coord_fixed(ratio=1) +
    #geom_text(data=region_annot_label, aes(x=median_x, y=median_y, label=division), size=1) +
    theme_void(base_size=7) +
    scale_fill_manual(values = division_colors, labels=proper_region_names,
                      breaks=names(proper_region_names), name = "Region") +
    theme(legend.position = "none",
          axis.title = element_blank())
  })

vs <- c("VL", "V3", "AQ", "V4")
fiber_tracts <- c("cm", "eps", "lfbs", "mfbs", "scwm")

merfish_cells_section <- merfish_cells %>% 
  filter(brain_section_label %in% sections_oi[1], x < 5.5,
         !division %in% c("brain-unassigned", "unassigned")) %>% 
  mutate(class = str_replace_all(class, "[- ]", ""),
         # rename any divisions in vs or fiber tracts to fiber tracts/VS
         division = ifelse(division %in% fiber_tracts | division %in% vs | division == "fiber tracts-unassigned",
                           "fiber tracts/VS", division))
         
merfish_cells_label <- merfish_cells_section %>% 
  group_by(division) %>% 
  summarise(median_x = median(x),
            median_y = median(y)) %>% 
  # rename fiber tracts/VS to ""
  mutate(division = ifelse(division == "fiber tracts/VS", "", division))

p_merfish <- ggplot(merfish_cells_section,
       aes(x=x, y=y, color=division)) +
  geom_point(shape=16, size=0.4, stroke=0, show.legend = FALSE) +
  coord_fixed(ratio=1) +
  scale_y_reverse() +
  # geom_text(data=merfish_cells_label, aes(x=median_x, y=median_y, label=division),
  #           size = 1, color="black") +
  theme_void(base_size=7) +
  scale_color_manual(values = division_colors) +
  theme(legend.position = "none",
        axis.title = element_blank())

p_region_ffpe <- p_regions[[1]] +
  guides(fill=guide_legend(nrow=1)) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.key.spacing.x = unit(0.15, "cm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=5, margin = margin(l=3)),)
p_regions_all <- p_merfish + p_region_ffpe + p_regions[[2]]
p_regions_all
ggsave("visium_hd_brain_combined/plots/merfish_vs_visiumhd_region_annotations.pdf",
       plot = p_regions_all,
       width = 7, height = 5)

ggplot(merfish_cells_section %>% filter(brain_section_label == sections_oi[1]),
       aes(y=y_section, x=x_section, color=class)) +
  geom_point(size=0.1) + scale_y_reverse() + coord_fixed(ratio=1) +
  #geom_label(data=region_annot_label, aes(x=median_x, y=median_y, label=acronym_lvl5), size=5) +
  scale_color_manual(values = celltype_colors) +
  coord_fixed() +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_blank())

merfish_cells_filtered <- merfish_cells %>% 
  filter(division %in% common_regions) %>%
  rename(region = division) %>% 
  mutate(brain_section_id = str_replace_all(brain_section_label, ".*\\.(\\d+)$", "\\1")) %>% 
  mutate(class = str_replace_all(class, "[- ]", "")) %>% 
  mutate(class = factor(class, levels = names(celltype_colors)),
         region = factor(region, levels = names(proper_region_names)))

# Create barplot per region
merfish_cells_summ <- merfish_cells_filtered %>% 
  group_by(brain_section_id, region, class, source) %>%
  summarise(count = n()) %>%
  group_by(brain_section_id, region, source) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup() %>% 
  mutate(class = factor(class, levels = names(celltype_colors)),
         region = factor(region, levels = names(proper_region_names)))

# Boxplot of proportions per region
deconv_props_region_summ <- deconv_props_df %>% 
  left_join(region_annotations_df %>% 
              select(barcode, division, source),
            by = c("spot" = "barcode", "source")) %>% 
  group_by(celltype, source, division) %>%
  summarise(mean_prop = mean(proportion)) %>% 
  mutate(celltype = str_replace_all(celltype, "^X", "")) %>% 
  rename(region = division, class=celltype) %>% 
  filter(region %in% common_regions) %>% 
  mutate(class = factor(class, levels = names(celltype_colors)),
         region = factor(region, levels = names(proper_region_names)))

merfish_cells_summ_top10  <- merfish_cells_summ %>% 
  group_by(region, class) %>% 
  summarise(median_prop = median(prop)) %>% 
  slice_max(order_by = median_prop, n = 10) %>% 
  select(-median_prop) %>% 
  inner_join(merfish_cells_summ,
             by = c("region", "class")) 
  
deconv_props_region_summ_top10 <- merfish_cells_summ_top10 %>% group_by(region) %>% 
  distinct(class) %>% 
  inner_join(deconv_props_region_summ,
             by = c("class", "region"))

return_boxplot <- function(merfish_dataset, deconv_props_df, regions_to_exclude, nrow=2) {
  ggplot(merfish_dataset %>% filter(!region %in% regions_to_exclude),
         aes(x = class, y = prop, color=class, group=interaction(class, source))) +
    geom_boxplot(linewidth=0.2, outlier.shape = 16, outlier.size = 0.5, outlier.stroke = 0, show.legend = FALSE, fill="white") +
    geom_point(data = deconv_props_df %>% filter(source == "ffpe", !region %in% regions_to_exclude, mean_prop > 0),
               aes(x=class, y=mean_prop, shape = source),
               position = position_nudge(x=-0.1), size=1, stroke=0.25, fill = "white", show.legend = TRUE) +
    geom_point(data = deconv_props_df %>% filter(source == "fresh_frozen", !region %in% regions_to_exclude, mean_prop > 0),
               aes(x=class, y=mean_prop, shape = source),
               position = position_nudge(x=0.1), size=1, stroke=0.25, fill="white") +
    ggh4x::facet_wrap2(~region, scales = "free_x",
                       labeller = as_labeller(proper_region_names),
                       nrow = nrow
    ) +
    scale_shape_manual(values = c(22, 23), labels = c("FFPE", "Fresh Frozen")) +
    scale_color_manual(values = celltype_colors, labels = proper_celltype_names, drop=FALSE) +
    #scale_fill_manual(values = c(celltype_colors), guide = "none") +
    scale_x_discrete(labels = function(x) str_replace_all(x, "^([0-9]+)([a-zA-Z0-9]+)$", "\\1")) +
    scale_y_continuous(breaks = seq(0, 1, by=0.2)) +
    # Override aes for legend
    guides(color = guide_legend(override.aes = list(shape = 15, size=3), nrow=6, order=1), 
           shape = guide_legend(override.aes = list(size=1.5), ncol=1, order=2)) +
    labs(y = "Cell type proportion", color = "Cell type", shape = "VisiumHD Dataset") +
    theme_bw(base_size=8) +
    theme(legend.position = "bottom",
          legend.key.size = unit(0.3, "cm"),
          legend.key.spacing.x = unit(0.15, "cm"),
          legend.key.spacing.y = unit(0, "cm"),
          legend.title.position = "top",
          legend.title = element_text(size=6),
          legend.text = element_text(size=5, margin = margin(l=3)),
          strip.background = element_blank(),
          strip.text.x.top = element_text(size=6),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_line(linewidth=0.1),
          panel.grid.major.y = element_line(linewidth=0.2),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.4),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=6),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.2),
          axis.text.y = element_text(size=5),
          axis.text.x = element_text(size=5))
}
                         
p_top10 <- return_boxplot(merfish_cells_summ_top10,
                           deconv_props_region_summ_top10,
                           regions_to_exclude = c("CTXsp", "PAL"))
p_top10

ggsave("visium_hd_brain_combined/plots/deconv_vs_merfish_boxplot_top10_celltypes_per_region.pdf",
       device = cairo_pdf,
       width = 7, height = 5.25)

# Supplementary -> remaining cell types
merfish_cells_summ_bottom <- merfish_cells_summ %>% 
  group_by(region, class) %>% 
  summarise(median_prop = median(prop)) %>% 
  ungroup() %>% 
  tidyr::complete(region, class, fill=list(median_prop=0)) %>% 
  group_by(region) %>% 
  arrange(-median_prop) %>%
  slice_tail(n=24) %>% 
  select(-median_prop) %>% 
  left_join(merfish_cells_summ,
             by = c("region", "class"))

deconv_props_summ_bottom <- merfish_cells_summ_bottom %>% 
  group_by(region) %>%  distinct(class) %>% 
  inner_join(deconv_props_region_summ,
             by = c("class", "region")) 
  

p_bottom10 <- return_boxplot(merfish_cells_summ_bottom %>% filter(!is.na(prop)),
                             deconv_props_summ_bottom,
                          regions_to_exclude = c("CTXsp", "PAL"),
                          nrow=3)
p_bottom10
ggsave(filename = "visium_hd_brain_combined/plots/deconv_vs_merfish_boxplot_bottom_celltypes_per_region.pdf",
       plot = p_bottom10,
       device = cairo_pdf,
       width = 7.25, height = 7)

p_ctxsp_pal <- return_boxplot(merfish_cells_summ %>% filter(region %in% c("CTXsp", "PAL")),
                              deconv_props_region_summ %>%
                                filter(region == "PAL" & source == "fresh_frozen" | region == "CTXsp" & source == "ffpe"),
                             regions_to_exclude = setdiff(common_regions, c("CTXsp", "PAL")),
                             nrow=2)

p_ctxsp_pal
ggsave("visium_hd_brain_combined/plots/deconv_vs_merfish_boxplot_celltypes_ctxpal.pdf",
       plot = p_ctxsp_pal,
       device = cairo_pdf,
       width = 7.25, height = 5)



deconv_props_region_df$division %>% table(useNA = "ifany")

# combined_props_summ <- bind_rows(deconv_props_summ %>% rename(region = division),
#           merfish_cells_summ %>%
#             rename(celltype = class,
#                    proportion = prop,
#                    source = brain_section_label))
# 
# # Stacked barplot per region
# ggplot(combined_props_summ, aes(x = source, y = proportion, fill = celltype)) +
#   geom_bar(stat = "identity", width=0.6) +
#   theme_minimal(base_size = 8) +
#   scale_fill_manual(values = celltype_colors) +
#   scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
#   facet_grid(~region, scales = "free_x", space='free') +
#   theme(axis.title = element_blank(),
#         panel.grid.major.x = element_blank())

# Calculate FPR for each region
fprs <- merfish_cells_summ %>% distinct(region, class, source) %>% 
  #tidyr::complete(region, class, source) %>%
  count(region, class) %>% 
  tidyr::complete(region, class, fill=list(n=0)) %>% 
  filter(n == 0) %>% select(-n) %>% 
  # Merge with number of spots in each region
  left_join(region_annotations_df %>%
              group_by(division, source) %>%
              summarise(n_spots = n()) %>%
              rename(region = division),
            by = c("region"), relationship = "many-to-many") %>% 
  left_join(deconv_props_count,
            by = c("region" = "division", "class" = "celltype", "source")) %>% 
  # Replace NA with 0
  mutate(n = ifelse(is.na(n), 0, n),
         tn = n_spots-n) %>% 
  rename(fp = n) %>% 
  mutate(fpr = fp / (fp + tn)) %>%
  group_by(region, source) %>%
  summarise(mean_fpr = mean(fpr)) %>% 
  pivot_wider(names_from = source, values_from = mean_fpr) %>% 
  ungroup() %>% 
  # PUT NA FOR region == PAL and ffpe
  mutate(ffpe = ifelse(region == "PAL", NA, ffpe),
         fresh_frozen = ifelse(region == "CTXsp", NA, fresh_frozen))
fprs
# Save as csv
write.csv(fprs,
          file = "visium_hd_brain_combined/rds/deconv_fpr_per_region.csv",
          row.names = FALSE, quote=FALSE)

fprs %>% select(region, ffpe) %>% filter(region != "PAL") %>% pull(ffpe) %>% mean
fprs %>% select(region, fresh_frozen) %>% filter(region != "CTXsp") %>% pull(fresh_frozen) %>% mean

deconv_props_count <- deconv_props_df %>% filter(proportion > 0) %>% 
  left_join(region_annotations_df %>% 
              select(barcode, division, source),
            by = c("spot" = "barcode", "source")) %>% 
  count(division, source, celltype)
