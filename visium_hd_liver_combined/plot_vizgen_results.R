library(Seurat)
library(tidyverse)
library(patchwork)

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
bin_sizes <- c(8,16,32)
bin_size_strs <- sprintf("%03dum", bin_sizes)
bin_size <- 8

vizgen <- read.csv("visium_hd_liver_combined/vizgen/vizgen_grid_8um_celltype_proportions.csv")


vizgen_df <- vizgen %>%
  pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion")

vizgen_df_ranked <- vizgen %>%
  pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>% 
  filter(proportion > 0) %>% 
  group_by(grid_id) %>% 
  # Only get top two
  slice_max(order_by=proportion, n=2) %>%
  mutate(rank = data.table::frankv(proportion, order=-1),
         doublet_type = case_when(n() == 1 ~ "singlet",
                                  n() == 2 ~ "doublet"))

# Complete the dataframe
vizgen_df_ranked_complete <- vizgen_df_ranked %>% 
  # Recalculate proportions to become 1
  group_by(grid_id) %>%
  mutate(proportion = proportion / sum(proportion)) %>%
  ungroup %>% 
  tidyr::complete(grid_id, celltype, fill = list(proportion = 0, rank=NA, doublet_type=NA))

# Cell types present in each spot
rowSums(vizgen[,-1] > 0) %>% table

vizgen_summ <- vizgen_df_ranked_complete %>% 
  group_by(celltype) %>% 
  summarise(mean_prop=mean(proportion))

# Plot barplot
ggplot(vizgen_summ, aes(x=reorder(celltype, -mean_prop), y=mean_prop, fill=celltype)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Cell Type", y="Mean Proportion", title="Mean Cell Type Proportions in Vizgen Data") +
  scale_fill_brewer(palette="Set3") +
  guides(fill=FALSE)

# Plot barplot of total times a cell type is the top cell type
ggplot(vizgen_df %>% filter(rank == 1), aes(x=celltype, fill=celltype)) +
  geom_bar() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Cell Type", y="Count", title="Top Cell Types in Vizgen Data") +
  scale_fill_brewer(palette="Set3") +
  guides(fill=FALSE)

# Plot barplot total times a cell type is the second cell type
ggplot(vizgen_df %>% filter(rank == 2), aes(x=celltype, fill=celltype)) +
  geom_bar() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Cell Type", y="Count", title="Second Cell Types in Vizgen Data") +
  scale_fill_brewer(palette="Set3") +
  guides(fill=FALSE)


# Plot boxplot
ggplot(vizgen_df_ranked %>% filter(proportion > 0), aes(x=reorder(celltype, -proportion, FUN=median), y=proportion, fill=celltype)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Cell Type", y="Proportion", title="Cell Type Proportions in Vizgen Data") +
  scale_fill_brewer(palette="Set3") +
  guides(fill=FALSE)

vizgen_summ

# Plot the squares
centroids <- read.csv("visium_hd_liver_combined/vizgen/vizgen_grid8um_centroids.csv")
square_size <- bin_size*1000/108

# merge with coordinates
vizgen_df_coords <- vizgen_df_ranked %>% ungroup() %>% 
  left_join(centroids, by = "grid_id") %>% 
  rename(x = centroid_x, y = centroid_y) %>% 
  mutate(n_celltypes = n(), .by = "grid_id",
         square_size = square_size,
         bin_size = bin_size)

roi <- c(10000,15000,8000,12000)

vizgen_df_coords <- vizgen_df_coords %>% 
  filter(x >= roi[1], x <= roi[2],
         y >= roi[3], y <= roi[4])

vizgen_df_square <- vizgen_df_coords %>%
  filter(n_celltypes == 1) %>% 
  mutate(group = as.character(grid_id)) %>% 
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

# Create barplot in squares
vizgen_df_barplot <- vizgen_df_coords %>% 
  filter(n_celltypes == 2) %>% 
  group_by(grid_id) %>% arrange(grid_id, desc(proportion)) %>% 
  mutate(rank = row_number(),
         group = paste0(grid_id, "_", rank)) %>% 
  mutate(x1 = case_when(rank == 1 ~ y - (square_size / 2) + (proportion*square_size),
                        rank == 2 ~ y + (square_size / 2) - (proportion*square_size)),
         y1 = x - square_size / 2,
         x2 = case_when(rank == 1 ~ y - (square_size / 2) + (proportion*square_size),
                        rank == 2 ~ y + (square_size / 2) - (proportion*square_size)),
         y2 = x + square_size / 2,
         x3 = case_when(rank == 1 ~ y - square_size / 2,
                        rank == 2 ~ y + square_size / 2),
         y3 = x + square_size / 2,
         x4 = case_when(rank == 1 ~ y - square_size / 2,
                        rank == 2 ~ y + square_size / 2),
         y4 = x - square_size / 2
  ) %>% 
  rename(coord_x = x, coord_y = y) %>%
  pivot_longer(cols = c(x1, y1, x2, y2, x3, y3, x4, y4),
               names_to = c(".value", "corner"),
               names_pattern = "(x|y)([1-4])")

vizgen_df_shapes <- bind_rows(vizgen_df_square, vizgen_df_barplot)
  #mutate(celltype = factor(celltype, levels = celltype_order))

ggplot(vizgen_df_shapes, aes(x = x, y = y)) +
  geom_polygon(aes(fill = celltype, group = group), show.legend = TRUE) +
  # White border
  geom_tile(data = vizgen_df_coords %>% distinct(x, y),
            aes(x = y, y = x), height = square_size, width = square_size,
            fill = NA, color = "white", inherit.aes = FALSE) +
  #scale_fill_manual(values = color_palette) +
  theme_void() +
  scale_y_reverse() +
  coord_fixed() +
  #ggtitle(paste0(bs, "\u00b5m")) +
  guides(fill = guide_legend(ncol=1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust = 0.5))

# Try with another file
vizgen_area <- read.csv("visium_hd_liver_combined/vizgen/vizgen_grid_8um_celltype_proportions_by_area.csv")

vizgen_area_df <- vizgen_area %>%
  pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion")

vizgen_area_summ <- vizgen_area_df %>%
  group_by(celltype) %>% 
  summarise(mean_prop=mean(proportion))

# Plot barplot
ggplot(vizgen_area_summ, aes(x=reorder(celltype, -mean_prop), y=mean_prop, fill=celltype)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Cell Type", y="Mean Proportion", title="Mean Cell Type Proportions in Vizgen Data by Area") +
  scale_fill_brewer(palette="Set3") +
  guides(fill=FALSE)

# Plot boxplot
ggplot(vizgen_area_df %>% filter(proportion > 0), aes(x=reorder(celltype, -proportion, FUN=median), y=proportion, fill=celltype)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Cell Type", y="Proportion", title="Cell Type Proportions in Vizgen Data by Area") +
  scale_fill_brewer(palette="Set3") +
  guides(fill=FALSE)

vizgen_area_summ
vizgen_summ

# Calculate moran's I on KCs
vizgen_KCs <- vizgen[, c("grid_id", "Kuppfer.cells")]
  #inner_join(centroids, by="grid_id")
# KC_dists <- as.matrix(dist(cbind(vizgen_KCs$centroid_x, vizgen_KCs$centroid_y)))
# KC_dists.inv <- 1/KC_dists
# diag(KC_dists.inv) <- 0
# KC_dists.inv[1:5, 1:5]
# ape::Moran.I(vizgen_KCs$Kuppfer.cells, KC_dists.inv)

library(spdep)
library(sf)
grid_shp <- read_sf("visium_hd_liver_combined/vizgen/grid_8um.shp") %>% 
  mutate(FID = FID + 1) %>% 
  filter(FID %in% vizgen_KCs$grid_id)
grid_nb <- poly2nb(grid_shp, queen=TRUE, snap=1)
lw <- nb2listw(grid_nb, style="W", zero.policy=TRUE)

I <- moran(vizgen_KCs$Kuppfer.cells, lw, length(grid_nb), Szero(lw))[1]
I #0.32 for 8um
moran.test(vizgen_KCs$Kuppfer.cells,lw, alternative="greater") #p-value < 0.001


deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))

deconv_props <- deconv_props_rank %>% 
  filter(dataset == "caw009", bin_size == "008um") %>% 
  ungroup() %>% 
  select(spot, celltype, proportion) %>% 
  tidyr::complete(spot, celltype, fill = list(proportion = 0))
  
deconv_props_KC <- deconv_props %>%
  filter(celltype == "Kupffercells")

dataset <- "_caw009"
bin_size_str <- "008um"
bin_size <- 8

#visium_obj <- readRDS("data/Visium_HD_Liver_CAW009/Visium_HD_Liver_CAW009_008um.rds")
#square_size <- visium_obj@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

vishd_shapes <- read_sf("visium_hd_liver_combined/rds/shapes/CAW009_square_008um.shp")
vishd_location <- read.csv("visium_hd_liver_combined/rds/shapes/CAW009_square_008um_locationid.csv", row.names=1) %>% 
  rownames_to_column("spot")

vishd_shapes_sub <- vishd_shapes %>%
  inner_join(vishd_location %>% select(spot, location_id),
            by=c("location_i"="location_id")) %>% 
  filter(spot %in% deconv_props_KC$spot)


grid_vishd_nb <- poly2nb(vishd_shapes_sub, queen=TRUE, snap=0.1)
lw_vishd <- nb2listw(grid_vishd_nb, style="W", zero.policy=TRUE)

vishd_I <- moran(deconv_props_KC$proportion, lw_vishd, length(grid_vishd_nb), Szero(lw_vishd))[1]
vishd_I #0 for vishd
moran.test(deconv_props_KC$proportion,lw_vishd, alternative="greater") #p-value = 0.53


#### Cells per grid

ncells_df <- lapply(bin_sizes, function(bin_size) {

  read.csv(paste0("visium_hd_liver_combined/vizgen/intersection_result_", bin_size, "um.csv")) %>% 
    group_by(grid_id) %>% 
    summarise(n_cells = n(),
              n_celltypes = n_distinct(annotation_own_score_genes)) %>% 
    mutate(bin_size = bin_size)
}) %>% bind_rows()

ncells_df_long <- ncells_df %>%
  pivot_longer(cols=c(n_cells, n_celltypes), names_to="metric", values_to="value") 

# Plot of number of cells and number of cell types per bin size
ggplot(ncells_df_long, aes(x=value)) +
  geom_histogram(binwidth=1) +
  ggh4x::facet_grid2(metric~bin_size, scales="free") 

# Get the proportions plot

vizgen <- read.csv("visium_hd_liver_combined/vizgen/vizgen_grid_8um_celltype_proportions.csv")

props_df <- lapply(bin_sizes, function(bin_size) {
  read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size, "um_celltype_proportions.csv")) %>%
    pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>%
    mutate(bin_size = bin_size)
}) %>% bind_rows()

vizgen_ori_props <- read.csv("visium_hd_liver_combined/vizgen/sdata_polyT_obs.csv") %>% 
  select(cells, annotation_own_score_genes) %>% 
  group_by(annotation_own_score_genes) %>%
  summarise(n_cells = n()) %>% 
  mutate(proportion = n_cells / sum(n_cells),
         bin_size = 0)

props_df_summ <- props_df %>%
  group_by(bin_size, celltype) %>% 
  summarise(mean_prop=mean(proportion)) %>% 
  bind_rows(vizgen_ori_props %>% rename(celltype = annotation_own_score_genes) %>% select(bin_size, celltype, mean_prop=proportion)) %>% 
  # remove space and dot from name of celltype
  mutate(celltype = gsub("[ .]", "", celltype))

# Plot barplot
ggplot(props_df_summ, aes(x=reorder(celltype, -mean_prop), y=mean_prop, fill=celltype)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="Cell Type", y="Mean Proportion", title="Mean Cell Type Proportions in Vizgen Data by Bin Size") +
  scale_fill_brewer(palette="Set3") +
  guides(fill=FALSE) +
  facet_wrap(~bin_size, nrow=1)


