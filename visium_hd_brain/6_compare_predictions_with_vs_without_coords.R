library(tidyverse)

# Compare doublet mode with and without spatial coordinates
props_with_coord <- read.table("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um",
                               header = TRUE)
props_wo_coord <- read.table("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um_converted_doublet",
                             header = TRUE)

doublet_info_with_coord <- read.table("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um_doublet_info.tsv",
                                      header = TRUE)
doublet_info_wo_coord <- read.table("visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um_converted_doublet_doublet_info.tsv",
                                    header = TRUE)

# Calculate correlation of proportions
mean(diag(cor(props_with_coord, props_wo_coord))) # 0.97

# Concordance of most abundant cell type?
max_celltypes <- data.frame(spot = doublet_info_with_coord$spot,
                            celltype_with_coord = colnames(props_with_coord)[max.col(props_with_coord)],
                            celltype_wo_coord = colnames(props_wo_coord)[max.col(props_wo_coord)])
max_celltypes %>% filter(celltype_with_coord == celltype_wo_coord) %>%
  nrow / nrow(max_celltypes) # 0.985

# Confusion matrix
celltypes <- colnames(props_with_coord)
conf_matrix <- caret::confusionMatrix(max_celltypes$celltype_with_coord %>% factor(levels=celltypes),
                                      max_celltypes$celltype_wo_coord %>% factor(levels=celltypes))
pheatmap::pheatmap(conf_matrix$table, scale="row",
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100),
                   cluster_rows = F, cluster_cols = FALSE,
                   main = "Doublet mode with (rows) vs without (cols) coordinates predictions",
                   filename = "visium_hd_brain/plots/conf_matrix_doublet_with_vs_without_coords.png")
# Write table
conf_matrix$table %>% 
  write.table("visium_hd_brain/plots/conf_matrix_doublet_with_vs_without_coords.tsv",
              sep = "\t", quote = FALSE)

# Compare doublet class predictions
classes <- c("singlet", "doublet_certain", "doublet_uncertain", "reject")
spot_classes <- data.frame(spot = doublet_info_with_coord$spot,
                           spot_class_with_coord = doublet_info_with_coord$spot_class,
                           spot_class_wo_coord = doublet_info_wo_coord$spot_class)
spot_classes %>% filter(spot_class_with_coord == spot_class_wo_coord) %>%
  nrow / nrow(spot_classes) # 0.98
conf_matrix <- caret::confusionMatrix(spot_classes$spot_class_with_coord %>% factor(levels=classes),
                                      spot_classes$spot_class_wo_coord %>% factor(levels=classes))

pheatmap::pheatmap(conf_matrix$table, scale="row",
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100),
                   cluster_rows = F, cluster_cols = FALSE,
                   main = "Doublet classification with (rows) vs without (cols) coordinates class predictions",
                   filename = "visium_hd_brain/plots/conf_matrix_spotclass_with_vs_without_coords_class.png")
# Write table
conf_matrix$table %>% 
  write.table("visium_hd_brain/plots/conf_matrix_spotclass_with_vs_without_coords_class.tsv",
              sep = "\t", quote = FALSE)
classes_colors <- c("singlet"="forestgreen", "doublet_certain"="navyblue",
                    "doublet_uncertain"="orange", "reject"="red")
# For the spots that differ in maximum cell type prediction
# Are they more likely to be rejects?
max_celltypes %>% filter(celltype_with_coord != celltype_wo_coord) %>%
  inner_join(spot_classes, by = "spot") %>% 
  select(spot_class_with_coord, spot_class_wo_coord) %>%
  pivot_longer(cols = everything(),
               names_to = "coord",
               values_to = "spot_class") %>%
  mutate(spot_class = factor(spot_class, levels = classes),
         coord = str_replace(str_remove(coord, "spot_class_"), "wo", "without")) %>%
  # Stacked bar plot
  ggplot(aes(x = coord, fill = spot_class)) +
  geom_bar(width=0.6) +
  scale_fill_manual(values = classes_colors) +
  # Expand scale
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(axis.title.x = element_blank())
# Save plot
ggsave("visium_hd_brain/plots/barplot_spot_class_of_different_preds_with_vs_without_coords.png",
       width = 5, height = 4)
