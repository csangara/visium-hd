library(Seurat)
library(tidyverse)
library(patchwork)
library(sf)

plot_path <- paste0("visium_hd_liver_combined/plots/")
bin_sizes <- c(8, 16)
n_bins <- 4

vizgen_df <- lapply(bin_sizes, function(bs) {
  vizgen <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_",
                            bs, "um_celltype_proportions.csv"))
  
  vizgen %>%
    pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>% 
    mutate(celltype = gsub("\\.", "", celltype)) %>% 
    mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype),
           bin_size = bs) %>% 
    filter(proportion > 0) %>% 
    ungroup()
}) %>% bind_rows()

segmented_cells <- read_sf("visium_hd_liver_combined/vizgen/segmentation_mask_boundaries.shp") 

intersection <- read.csv("visium_hd_liver_combined/vizgen/intersection_result_16um.csv")
grid_8um <- read_sf(paste0("visium_hd_liver_combined/vizgen/grid_", bin_sizes[1], "um.shp"))
grid_16um <- read_sf(paste0("visium_hd_liver_combined/vizgen/grid_", bin_sizes[2], "um.shp"))

grid_oi <- 8815
cells_oi <- intersection %>% filter(grid_id == grid_oi) 
  # left_join(vizgen_ori_annot %>% select(cells, celltype),
  #           by=c("index"="cells"))

segmented_cells_subset <- segmented_cells %>% mutate(index=as.integer(index)) %>% 
  right_join(cells_oi, by=c("index"="segmentation_id")) %>% 
  rename(celltype = annotation_own_score_genes)

bbox <- st_bbox(grid_16um %>% filter(FID+1 == grid_oi) %>% select(geometry))

# Subset grid_8um to the bbox
grid_8um_subset <- grid_8um %>% 
  filter(st_intersects(geometry, st_as_sfc(bbox), sparse = FALSE)[,1])

grid_16um_subset <- grid_16um %>% 
  filter(st_intersects(geometry, st_as_sfc(bbox), sparse = FALSE)[,1])

ggplot() +
  geom_sf(aes(fill=as.factor(celltype)), color="black", data=segmented_cells_subset) +
  geom_sf(color="black", fill=NA, data=grid_8um_subset) +
  geom_sf(color="red", fill=NA, size=1, data=grid_16um_subset) +
  lims(x=c(bbox$xmin, bbox$xmax), y=c(bbox$ymin, bbox$ymax)) +
  theme_void()

# Read parquet file
tissue_positions <- arrow::read_parquet(
  paste0("visium_hd_liver_combined/vizgen/tissue_positions_",
         bin_sizes[1], "um_", bin_sizes[2], "um.parquet")) %>% 
  select(barcode, containing_barcode) %>% 
  # Filter to those with 4 barcodes per containing barcode
  group_by(containing_barcode) %>%
  filter(n() == n_bins) %>% ungroup()
  
vizgen_df_subset <- vizgen_df %>%
  filter(bin_size %in% bin_sizes)

vizgen_df_celltype <- vizgen_df_subset %>%
  filter(celltype == "Kupffercells")

mother_bins <- tissue_positions %>% 
  select(barcode, containing_barcode) %>% 
  inner_join(vizgen_df_celltype %>% filter(bin_size == bin_sizes[1]) %>%
              ungroup() %>% select(grid_id, proportion, celltype) %>%
               rename(prop_child = proportion),
             by = c("barcode" = "grid_id")) %>% 
  inner_join(vizgen_df_celltype %>% filter(bin_size == bin_sizes[2]) %>%
              ungroup() %>% select(grid_id, proportion, celltype) %>% 
              rename(prop_mother = proportion),
            by = c("containing_barcode" = "grid_id")) %>% 
  # replace prop by 0 if it's NA
  mutate(prop_child = ifelse(is.na(prop_child), 0, prop_child))

mother_bins_range <- mother_bins %>% 
  group_by(containing_barcode) %>% 
  mutate(range = max(prop_child) - min(prop_child))

ggplot(mother_bins_range %>% ungroup() %>% 
         distinct(containing_barcode, prop_mother, range), 
       aes(x=prop_mother, y=range, text=containing_barcode)) +
  geom_abline(slope=c(1), intercept=0, color="gray70", linetype="dashed") +
  geom_point(size=0.5) +
  labs(x=paste0("Proportions of Kupffer cells in mother bins (", bin_sizes[2], "um)"),
       y=paste0("Range of proportions across ", n_bins, " child bins (", bin_sizes[1], "um)")) +
  # Geom ab line
  theme_classic()

# Bin props_mother into 20 bins
mother_bins_binned <- mother_bins %>% 
  mutate(prop_mother_bin = cut(prop_mother, breaks=seq(0, 1, by=0.05), include.lowest=TRUE)) %>% 
  group_by(prop_mother_bin) %>% 
  mutate(bin_size = n()) %>% ungroup()

p <- ggplot(mother_bins_binned,
       aes(x=prop_mother_bin, y=prop_child)) +
  geom_boxplot()
p

