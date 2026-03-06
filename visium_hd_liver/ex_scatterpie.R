library(scatterpie)
library(Seurat)
library(tidyverse)

##### OTHER PLOTS - SCATTERPIE & SCATTERTRI #####
##### NOT USED #####

dataset <- "" #or "_caw009" or ""
data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
plot_path <- paste0("visium_hd_liver", dataset, "/plots/")

color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "CentralVeinEndothelialcells" = "#FED8B1FF",
                   "LSECs" = "#FBB05FFF",
                   "PortalVeinEndothelialcells" = "#CC7722FF",
                   "LymphaticEndothelialcells" = "#8F4716FF",
                   "Cholangiocytes" = "#C61B84FF",
                   "HsPCs" = "#F19FC3FF",
                   "Stellatecells" = "#A31A2AFF",
                   "Mesothelialcells" = "#D0110BFF",
                   "Fibroblasts" = "#E45466FF",
                   "CapsularFibroblasts" = "#D46F6CFF",
                   "Kupffercells" = "#5DA6DBFF",
                   "MonocytesMonocytederivedcells" = "#a3daf3",
                   "cDC1s" = "#893A86FF",
                   "cDC2s" = "#893A86FF",
                   "pDCs" = "#893A86FF",
                   "MigcDCs" = "#893A86FF",
                   "Bcells" = "#9C7EBAFF",
                   "NKcells" = "#4A6E34FF",
                   "Tcells" = "#3AB04AFF",
                   "ILC1s" = "#A3D7BAFF",
                   "Basophils" = "#191919",
                   "Neutrophils" = "#727272")

celltype_order <- names(color_palette)
bin_size <- 8 # 8, 16, or 32

# Fill the text with 0
bin_size_str <- sprintf("%03dum", bin_size)
visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                             bin_size_str, ".rds"))
dim(visium_obj) # 19059 genes x 462269 spots

ext <- "_annot_cd45"
deconv_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                  "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                  bin_size_str, ext),
                           header = TRUE)
dim(deconv_props) #460942 spots

removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,paste0("nCount_Spatial.", bin_size_str)] < 100)]

# Subset visium_obj to only include spots that were not removed
visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]

dim(visium_obj_subset) # 19059 genes x 460942 spots

# Add rownames to deconv_props
rownames(deconv_props) <- colnames(visium_obj_subset)

# Let's get a region of interest
roi <- c(6000, 10000, 30000, 34000)
roi_name <- paste0("_", paste0(sprintf("%02d", roi/1000), collapse=""))

cells_roi <- GetTissueCoordinates(visium_obj_subset) %>%
  filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)

visium_obj_roi <- visium_obj_subset %>% .[, colnames(.) %in% cells_roi]
square_size <- visium_obj_subset@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
deconv_props_roi <- deconv_props[rownames(deconv_props) %in% cells_roi, ]

deconv_props_scatterpie <- deconv_props_roi  %>%
  rownames_to_column("spot") %>%
  left_join(GetTissueCoordinates(visium_obj_roi), by = c("spot" = "cell")) %>% 
  mutate(radius = square_size / 2)

p_scatterpie <- ggplot() +
  geom_scatterpie(aes(r=radius, x = y, y = x, group=spot),
                  data = deconv_props_scatterpie,
                  cols=unique(deconv_props_df$celltype), color=NA) +
  scale_fill_manual(values = color_palette) +
  coord_equal() +
  scale_y_reverse() +
  theme_void()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))
ggsave(paste0(plot_path, "scatter_test/spatialscatterpie_ROI", roi_name, "_", bin_size_str, ".png"),
       p_scatterpie, width = 10, height = 8, bg = "white")

deconv_props_df_square <- deconv_props_df %>%
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

deconv_props_df_tri <- deconv_props_df %>% 
  filter(n_celltypes == 2) %>% 
  group_by(spot) %>% arrange(spot, desc(proportion)) %>% 
  mutate(rank = row_number(),
         group = paste0(spot, "_", rank)) %>% 
  # Select only lower proportion
  # Draw half triangles for each
  # Higher proportion in the lower half (since we will flip it later)
  mutate(x1 = y + square_size / 2,
         y1 = x - square_size / 2,
         x2 = y - square_size / 2,
         y2 = x + square_size / 2,
         x3 = case_when(rank == 1 ~ y - square_size / 2,
                        rank == 2 ~ y + square_size / 2),
         y3 = case_when(rank == 1 ~ x - square_size / 2,
                        rank == 2 ~ x + square_size / 2)
         
  ) %>% 
  rename(coord_x = x, coord_y = y) %>%
  pivot_longer(cols = c(x1, y1, x2, y2, x3, y3),
               names_to = c(".value", "corner"),
               names_pattern = "(x|y)([1-3])")


deconv_props_df_shapes <- bind_rows(deconv_props_df_square, deconv_props_df_tri) %>% 
  mutate(celltype = factor(celltype, levels = celltype_order))

# Create the ggplot
ggplot(deconv_props_df_shapes, aes(x = x, y = y, group = group)) +
  geom_polygon(aes(fill = celltype)) +
  scale_fill_manual(values = color_palette) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))

ggsave(paste0(plot_path, "scatter_test/spatialtriplot_doublet_celltype_ROI", roi_name, "_",
              bin_size_str, ".png"),
       width = 10, height = 8, bg = "white")
