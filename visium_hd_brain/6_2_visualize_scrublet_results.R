library(Seurat)
library(tidyverse)
library(RColorBrewer)

visium_obj <- readRDS("data/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds")
dim(visium_obj) # 19059 genes x 393543 spots

doublet <- read.csv("visium_hd_brain/Visium_HD_MouseBrain_008um/scrublet_Visium_HD_MouseBrain_008um.csv",
                    header = TRUE) %>% rename(spot_id = X)
dim(doublet)

# Get region annotations
region_annotations <- readRDS("data/Visium_HD_MouseBrain/tissue_positions_with_annotations_008um.rds")

# Filter to only spots that were deconvolved
region_annotations <- region_annotations %>% filter(!is.na(acronym_lvl6)) %>% 
  mutate(region_broad = case_when(
    acronym_lvl6 %in% c("VIS", "PTLp", "SS", "AUD", "TEa", "RHP") ~ "Cerebral cortex",
    acronym_lvl6 %in% c("STRd", "sAMY") ~ "Cerebral nuclei",
    acronym_lvl6 %in% c("DORpm", "DORsm") ~ "Thalamus",
    acronym_lvl6 %in% c("MEZ", "LZ", "PVR", "PVZ") ~ "Hypothalamus",
    acronym_lvl6 %in% c("HIP") ~ "Hippocampus",
    TRUE ~ acronym_lvl6
  ))

# Filter deconvolved proportions to match region_annotations
doublet_props_df <- doublet %>% 
  filter(spot_id %in% region_annotations$barcode) %>% 
  inner_join(region_annotations %>% select(barcode, acronym_lvl6, pxl_col_in_lowres, pxl_row_in_lowres, region_broad),
             by = c("spot_id" = "barcode"))

# Summarise per region
doublet_props_summ <- doublet_props_df %>% group_by(region_broad, acronym_lvl6) %>%
  summarise(doublet_props = sum(predicted_doublet == "True")/n())

# Stacked barplot per region
p <- ggplot(doublet_props_summ, aes(x = acronym_lvl6, y = doublet_props)) +
  geom_bar(stat = "identity", width=0.6) +
  theme_minimal(base_size = 8) +
  #scale_fill_manual(values = col_vector) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  facet_grid(~region_broad, scales = "free_x", space='free') +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
p
# ggsave(paste0("visium_hd_brain/plots/barplot_avg_by_region", ext, ".png"), p,
#        width = 8, height = 6, bg = "white")

table(doublet$predicted_doublet) %>% prop.table()
# 0.86 False, 0.14 True

# scale_fill_manual(values=c("FALSE"="forestgreen", "TRUE"="navyblue"),
#                   labels=c("Singlet", "Doublet"),
#                   breaks=c("FALSE", "TRUE")) +

classes_colors <- c("False"="forestgreen", "True"="navyblue")
visium_obj$spot_class <- doublet$predicted_doublet

p_spot_class <- SpatialDimPlot(visium_obj, group.by = "spot_class",
                               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = classes_colors, labels = c("Singlet", "Doublet")) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right",
        legend.title = element_blank())
ggsave(paste0("visium_hd_brain/plots/spatialdimplot_scrublet_doublet.png"), p_spot_class,
       width = 8, height = 6, bg = "white")
