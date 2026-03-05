library(Seurat)
library(tidyverse)

visium_obj <- readRDS("data/Visium_HD_MouseBrain/Visium_HD_MouseBrain_008um.rds")
dim(visium_obj) # 19059 genes x 393543 spots

doublet <- read.csv("visium_hd_brain/Visium_HD_MouseBrain_008um/scrublet_Visium_HD_MouseBrain_008um.csv",
                    header = TRUE) %>% rename(spot_id = X)
dim(doublet)

# Get region annotations
region_annotations <- readRDS("data/Visium_HD_MouseBrain/tissue_positions_with_annotations_008um.rds")
region_annotations <- region_annotations %>% filter(!grepl("unassigned", division))

# Filter deconvolved proportions to match region_annotations
doublet_props_df <- doublet %>% 
  filter(spot_id %in% region_annotations$barcode) %>% 
  inner_join(region_annotations %>% select(barcode, division, pxl_col_in_lowres, pxl_row_in_lowres),
             by = c("spot_id" = "barcode"))

# Summarise per region
doublet_props_summ <- doublet_props_df %>% group_by(division) %>%
  summarise(doublet_props = sum(predicted_doublet == "True")/n())

# Stacked barplot per region
p <- ggplot(doublet_props_summ, aes(x = division, y = doublet_props)) +
  geom_bar(stat = "identity", width=0.6) +
  theme_minimal(base_size = 8) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme(axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
p

table(doublet$predicted_doublet) %>% prop.table()
# 0.86 False, 0.14 True

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
