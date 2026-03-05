library(tidyverse)

file_name <- "visium_hd_brain/Visium_HD_MouseBrain_008um/proportions_rctd_Visium_HD_MouseBrain_008um"

doublet_props <- read.table(file_name, header=TRUE)
full_props <- read.table(paste0(file_name, "_full"), header=TRUE)
                           
# Make confusion matrix
celltypes <- colnames(doublet_props)
conf_matrix <- caret::confusionMatrix(colnames(full_props)[max.col(full_props)] %>% factor(levels=celltypes),
                                      colnames(doublet_props)[max.col(doublet_props)] %>% factor(levels=celltypes))
pheatmap::pheatmap(conf_matrix$table, scale="row",
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdGy")))(100),
                   cluster_rows = F, cluster_cols = FALSE,
                   main = "Full mode (rows) vs Doublet mode (cols) predictions",
                   filename = "visium_hd_brain/plots/conf_matrix_full_vs_doublet.png")

# Get correlation of proportions between full and doublet mode
mean(diag(cor(full_props, doublet_props)))

# Overlap of most abundant cell type
max_celltypes <- data.frame(celltype_full = colnames(full_props)[max.col(full_props)],
                            celltype_doublet = colnames(doublet_props)[max.col(doublet_props)])
max_celltypes %>% filter(celltype_full == celltype_doublet) %>% nrow / nrow(max_celltypes)