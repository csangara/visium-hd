library(Seurat)
library(tidyverse)
source("visium_hd_brain_ff/0_utils.R")

visium_obj <- readRDS("data/Visium_HD_MouseBrain_FF/Visium_HD_MouseBrain_FF_008um.rds")
dim(visium_obj) # 19059 genes x 453820 spots

doublet <- read.csv("visium_hd_brain_ff/Visium_HD_MouseBrain_FF_008um/scrublet_Visium_HD_MouseBrain_FF_008um.csv",
                    header = TRUE) %>% rename(spot_id = X)
dim(doublet)

# # Get region annotations
# region_annotations <- readRDS("data/Visium_HD_MouseBrain/tissue_positions_with_annotations_008um.rds")
# 
# # Filter to only spots that were deconvolved
# region_annotations <- region_annotations %>% filter(!is.na(acronym_lvl6)) %>% 
#   mutate(region_broad = case_when(
#     acronym_lvl6 %in% c("VIS", "PTLp", "SS", "AUD", "TEa", "RHP") ~ "Cerebral cortex",
#     acronym_lvl6 %in% c("STRd", "sAMY") ~ "Cerebral nuclei",
#     acronym_lvl6 %in% c("DORpm", "DORsm") ~ "Thalamus",
#     acronym_lvl6 %in% c("MEZ", "LZ", "PVR", "PVZ") ~ "Hypothalamus",
#     acronym_lvl6 %in% c("HIP") ~ "Hippocampus",
#     TRUE ~ acronym_lvl6
#   ))
# 
# # Filter deconvolved proportions to match region_annotations
# doublet_props_df <- doublet %>% 
#   filter(spot_id %in% region_annotations$barcode) %>% 
#   inner_join(region_annotations %>% select(barcode, acronym_lvl6, pxl_col_in_lowres, pxl_row_in_lowres, region_broad),
#              by = c("spot_id" = "barcode"))
# 
# # Summarise per region
# doublet_props_summ <- doublet_props_df %>% group_by(region_broad, acronym_lvl6) %>%
#   summarise(doublet_props = sum(predicted_doublet == "True")/n())
# 
# # Stacked barplot per region
# p <- ggplot(doublet_props_summ, aes(x = acronym_lvl6, y = doublet_props)) +
#   geom_bar(stat = "identity", width=0.6) +
#   theme_minimal(base_size = 8) +
#   #scale_fill_manual(values = col_vector) +
#   scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
#   facet_grid(~region_broad, scales = "free_x", space='free') +
#   theme(axis.title = element_blank(),
#         panel.grid.major.x = element_blank(),
#         legend.position = "none")
# p
# ggsave(paste0("visium_hd_brain/plots/barplot_avg_by_region", ext, ".png"), p,
#        width = 8, height = 6, bg = "white")

table(doublet$predicted_doublet) %>% prop.table()
# 0.99 False, 0.08 True

doublet$predicted_doublet %>% .[. != ""] %>% table # Only four doublets...

# About 25% of bins have two nuclei (see 4_3), so let's find the threshold for 75% of data
q75_scrublet <- doublet$doublet_score %>% .[!is.na(.)] %>% quantile(probs = 0.75)
doublet <- doublet %>% mutate(predicted_doublet = case_when(
  is.na(doublet_score) ~ "",
  doublet_score >= q75_scrublet ~ "True",
  TRUE ~ "False"
))

classes_colors <- c("False"="forestgreen", "True"="navyblue")
visium_obj$spot_class <- doublet$predicted_doublet

# # Highlight spots with no prediction
# no_prediction_cells <- doublet %>% filter(is.na(doublet_score)) %>% pull(spot_id)
# SpatialDimPlot(visium_obj,
#                cells.highlight = no_prediction_cells,
#                image.alpha=0)

p_spot_class <- SpatialDimPlot(visium_obj %>% .[, .$spot_class != ""],
                               group.by = "spot_class",
                               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = classes_colors, labels = c("Singlet", "Doublet")) +
  xlim(orig_xlim) + ylim(orig_ylim) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right",
        legend.title = element_blank())
ggsave(paste0("visium_hd_brain_ff/plots/spatialdimplot_scrublet_doublet.png"), rotate_image(p_spot_class, -90),
       width = 8, height = 6, bg = "white")
