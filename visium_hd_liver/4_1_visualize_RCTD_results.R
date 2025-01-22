library(Seurat)
library(tidyverse)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

visium_obj <- readRDS("data/Visium_HD_Liver/Visium_HD_Liver_008um.rds")
dim(visium_obj) # 19059 genes x 462269 spots

# Plot in sections -> label regions of visium_obj
# Subset each run into four quadrants
coords <- GetTissueCoordinates(visium_obj)

# Get quadrants
x_mid <- mean(c(max(coords$x), min(coords$x)))
y_mid <- mean(c(max(coords$y), min(coords$y)))

# Label each cell per quadrant
coords <- coords %>% 
  mutate(quadrant = case_when(x > x_mid & y > y_mid   ~ "Q1",
                              x <= x_mid & y > y_mid  ~ "Q2",
                              x <= x_mid & y <= y_mid ~ "Q3",
                              TRUE ~ "Q4"))

visium_obj$quadrant <- coords$quadrant

for (i in 1:4){
  p_ncount <- SpatialFeaturePlot(visium_obj %>% .[, .$quadrant == paste0("Q", i)],
                                 "nCount_Spatial.008um",
                                 pt.size.factor = 5, image.alpha=0) +
    labs(fill = "nCount") +
    theme(legend.position = "right")
  
  #p_ncount
  ggsave(paste0("visium_hd_liver/plots/spatialfeatureplot_nCount_Q", i, ".png"),
         p_ncount,
         width = 20, height = 15, bg = "white") 
}

p_ncount <- SpatialFeaturePlot(visium_obj, "nCount_Spatial.008um",
                               pt.size.factor = 2, image.alpha=0) +
  labs(fill = "nCount") +
  theme(legend.position = "right")
#p_ncount
ggsave(paste0("visium_hd_liver/plots/spatialfeatureplot_nCount.png"),
       p_ncount,
       width = 40, height = 30, bg = "white") 

ext <- ""

deconv_props <- read.table(paste0("visium_hd_liver/Visium_HD_Liver_008um/proportions_rctd_Visium_HD_Liver_008um", ext),
                           header = TRUE)
dim(deconv_props) #460942 spots

# Get removed rows from the other file
removed_rows <- scan("visium_hd_liver/Visium_HD_Liver_008um/proportions_rctd_rows_removed", what="character") %>% 
  .[grepl("s_008um", .)]

length(removed_rows) # 1327

# Check if removed rows + leftover rows == total rows (yes)
length(removed_rows) + dim(deconv_props)[1] == dim(visium_obj)[2]

# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_rows)]

dim(visium_obj_subset) # 19059 genes x 278323 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

# Check counts of removed rows
visium_obj[, colnames(visium_obj) %in% removed_rows]$nCount_Spatial.008um %>% hist()

# Assign barcode to most abundant cell type per spot
all(rownames(deconv_props) == colnames(visium_obj_subset))
visium_obj_subset$celltype <- colnames(deconv_props)[max.col(deconv_props)]
visium_obj_subset$celltype %>% table %>% prop.table() %>% sort

p_celltype <- SpatialDimPlot(visium_obj_subset, group.by = "celltype",
               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = col_vector) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right")
p_celltype
ggsave(paste0("visium_hd_liver/plots/spatialdimplot_celltype", ext, ".png"), p_celltype,
       width = 8, height = 6, bg = "white")


# Plot by quadrant
for (i in 1:4){
  p_celltype <- SpatialDimPlot(visium_obj_subset %>% .[, .$quadrant == paste0("Q", i)],
                               group.by = "celltype", interactive = TRUE,
                               image.alpha = 0, pt.size.factor = 5, shape=21, stroke=0.1) +
    scale_fill_manual(values = col_vector) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme(legend.position = "none")
  
  #p_celltype
  ggsave(paste0("visium_hd_liver/plots/spatialdimplot_celltype_Q", i, ext, ".png"),
         p_celltype,
         width = 10, height = 8, bg = "white") 
}
# Visualize ROI in Q2
# xmin, ymin, xmax, ymax = 8000, 34151, 10000, 35900
cells_roi <- GetTissueCoordinates(visium_obj_subset) %>% filter(x > 8000 & x < 10000 & y > 34151 & y < 35900) %>% pull(cell)

SpatialDimPlot(visium_obj_subset %>% .[, .$quadrant == "Q2"],
               cells.highlight = cells_roi)

visium_obj_roi <- visium_obj_subset %>% .[, .$quadrant == "Q2"] %>%
  .[, colnames(.) %in% cells_roi]

SpatialDimPlot(visium_obj_roi,
               group.by = "celltype",
               image.alpha = 0, pt.size.factor = 30, shape=22, stroke=0.1) +
  scale_fill_manual(values = col_vector) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

#p_celltype
ggsave(paste0("visium_hd_liver/plots/spatialdimplot_ROI.png"),
       width = 10, height = 5, bg = "white") 


deconv_props_df <- data.frame(deconv_props) %>% 
  pivot_longer(cols = everything(), names_to = "celltype",
               values_to = "proportion")

# Cell type proportions barplot
deconv_props_summ <- deconv_props_df %>% 
  group_by(celltype) %>% 
  summarise(agg_proportion = mean(proportion))

# Barplot ordered by total abundance
p <- ggplot(deconv_props_summ, aes(x=reorder(celltype, agg_proportion), y=agg_proportion)) +
  geom_bar(stat="identity") +
  # Add text of value, rounded to two digits, only if value is > 0
  geom_text(aes(label=ifelse(agg_proportion > 0.001, round(agg_proportion, 3), "")), nudge_y = 0.03, size=2) +
  coord_flip() +
  theme_minimal(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  ggtitle("Average cell type proportions across tissue")

p
ggsave(paste0("visium_hd_liver/plots/barplot_avg_across_tissue", ext, ".png"), p,
       width = 8, height = 6, bg = "white")


# What is the distribution of cell types across the tissue?
p_boxplot_all <- ggplot(deconv_props_df %>% filter(proportion > 0.0001),
       aes(y=reorder(celltype, proportion, median), x=proportion)) +
  coord_flip() +
  geom_boxplot(alpha=0.1) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = 8) +
  theme(axis.title.x = element_blank(),
        # angled text
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),) +
  ggtitle("Cell type proportions across tissue (prop > 0.0001)")
#p_boxplot_all
ggsave(paste0("visium_hd_liver/plots/boxplot_all", ext, ".png"), p_boxplot_all,
       width = 15, height = 8, bg = "white")


# DOUBLET MODE
ext <- "" # "", "_converted_doublet"
doublet_props <- read.table(paste0("visium_hd_liver/Visium_HD_Liver_008um/proportions_rctd_Visium_HD_Liver_008um", ext),
                            header = TRUE)
doublet_info <- read.table(paste0("visium_hd_liver/Visium_HD_Liver_008um/proportions_rctd_Visium_HD_Liver_008um", ext, "_doublet_info.tsv"),
                           header = TRUE)

# Add rownames to deconv_props
rownames(doublet_props) <- colnames(visium_obj_subset)
all(rownames(doublet_props) == doublet_info$spot) # Check

table(doublet_info$spot_class)

classes_colors <- c("singlet"="forestgreen", "doublet_certain"="navyblue",
                    "doublet_uncertain"="orange", "reject"="red")
visium_obj_subset$spot_class <- doublet_info$spot_class
p_spot_class <- SpatialDimPlot(visium_obj_subset, group.by = "spot_class",
               image.alpha = 0, stroke=NA) +
  scale_fill_manual(values = classes_colors) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right")
p_spot_class
ggsave(paste0("visium_hd_liver/plots/spatialdimplot_spot_class", ext, ".png"), p_spot_class,
       width = 8, height = 6, bg = "white")
