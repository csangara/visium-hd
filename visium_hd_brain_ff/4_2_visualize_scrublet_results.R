library(Seurat)
library(tidyverse)
source("visium_hd_brain_ff/4_0_utils.R")

visium_obj <- readRDS("data/Visium_HD_MouseBrain_FF/Visium_HD_MouseBrain_FF_008um.rds")
dim(visium_obj) # 19059 genes x 453820 spots

doublet <- read.csv("visium_hd_brain_ff/Visium_HD_MouseBrain_FF_008um/scrublet_Visium_HD_MouseBrain_FF_008um.csv",
                    header = TRUE) %>% rename(spot_id = X)
dim(doublet)

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

# 153 spots have no prediction
table(doublet$predicted_doublet)

p_spot_class <- SpatialDimPlot(visium_obj %>% .[, .$spot_class != ""],
                               group.by = "spot_class",
                               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = classes_colors, labels = c("Singlet", "Doublet")) +
  xlim(orig_xlim) + ylim(orig_ylim) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right",
        legend.title = element_blank())
ggsave(paste0("visium_hd_brain_ff/plots/spatialdimplot_scrublet_doublet.png"),
       rotate_image(p_spot_class, -90),
       width = 8, height = 6, bg = "white")
