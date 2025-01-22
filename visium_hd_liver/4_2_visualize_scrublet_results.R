library(Seurat)
library(tidyverse)

visium_obj <- readRDS("data/Visium_HD_Liver/Visium_HD_Liver_008um.rds")
dim(visium_obj) # 19059 genes x 462269 spots

doublet <- read.csv("visium_hd_liver/Visium_HD_Liver_008um/scrublet_Visium_HD_Liver_008um.csv",
                    header = TRUE) %>% rename(spot_id = X)
dim(doublet)
table(doublet$predicted_doublet) # 37 doublets...

# For now, let's use the 0.75 quantile again
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
ggsave(paste0("visium_hd_liver/plots/spatialdimplot_scrublet_doublet.png"),
       width = 8, height = 6, bg = "white")
