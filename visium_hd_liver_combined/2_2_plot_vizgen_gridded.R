library(spdep)
library(sf)
source("scripts/scatterbarplot_function.R")
source("visium_hd_liver_combined/0_utils.R")

bin_sizes <- c(2, 4, 8, 16, 32)
make_combined_plot <- FALSE
first_run <- FALSE

#### PLOT SEGMENTED CELLS AND OVERLAY GRID ####
# Read in shape file
vizgen_ori_annot <- read.csv("visium_hd_liver_combined/vizgen/sdata_polyT_obs.csv") %>% 
  select(cells, annotation_own_score_genes) %>% 
  mutate(cells = as.character(cells)) %>% 
  rename(celltype = annotation_own_score_genes) %>%
  mutate(celltype = gsub(" ", "", celltype)) %>%
  mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype))

segmented_cells <- read_sf("visium_hd_liver_combined/vizgen/segmentation_mask_boundaries_roi.shp") %>% 
  left_join(vizgen_ori_annot %>% select(cells, celltype),
            by=c("index"="cells"))

if (make_combined_plot){
  #crd <- c(10500, 12500, 9700, 11000)
  crd <- c(10666, 12445, 9776, 11260) # from st_bbox
} else {
  crd <- c(10000,15000,8000,12000)
}

crd_poly <-  st_multipoint(
  matrix(c(crd[1], -crd[3],
           crd[2], -crd[3],
           crd[2], -crd[4],
           crd[1], -crd[4],
           crd[1], -crd[3]), ncol=2, byrow=TRUE)
) %>% st_cast("POLYGON")


p_grids <- lapply(bin_sizes, function(bin_size) {
  grid <- read_sf(paste0("visium_hd_liver_combined/vizgen/grid_", bin_size, "um_roi.shp"))
  grid_roi <- st_filter(grid, crd_poly, .predicate = st_within)
  lims <- st_bbox(grid_roi)
  
  ggplot() +
    geom_sf(aes(fill=celltype), color = "white", linewidth=0.1, data = segmented_cells) +
    geom_sf(fill=NA, color="gray40", linewidth=0.05, data = grid_roi) +
    theme_void() +
    coord_sf(xlim = c(lims$xmin, lims$xmax),
             ylim = c(lims$ymin, lims$ymax),
             expand = FALSE) +
    scale_fill_manual(values=color_palette_vizgen) +
    labs(title=paste0(bin_size, "\u03bcm")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5, size=8))
})

wrap_plots(p_grids, ncol=1) & theme(plot.title = element_blank())

if (!make_combined_plot){
  p_grids_wrapped <- wrap_plots(p_grids, nrow=2) &
    theme(plot.margin = margin(l=5, r=5))
  ggsave(paste0(plot_path, "segmented_cells_with_grid_overlay.pdf"),
         p_grids_wrapped,
         device = cairo_pdf,
         width=9, height=6, dpi=300)
}

#### COUNTS CELLS AND CELLTYPE PER GRID ####
if (first_run){
  ncells_df <- lapply(bin_sizes, function(bin_size) {
    read.csv(paste0("visium_hd_liver_combined/vizgen/intersection_result_", bin_size, "um.csv")) %>% 
      group_by(grid_id) %>% 
      summarise(n_cells = n(),
                n_celltypes = n_distinct(annotation_own_score_genes)) %>% 
      mutate(bin_size = bin_size)
  }) %>% bind_rows()
  saveRDS(ncells_df, "visium_hd_liver_combined/vizgen/ncells_ncelltypes_per_grid.rds")
} else {
  ncells_df <- readRDS("visium_hd_liver_combined/vizgen/ncells_ncelltypes_per_grid.rds")
}

ncells_df_long <- ncells_df %>%
  pivot_longer(cols=c(n_cells, n_celltypes), names_to="metric", values_to="value") 

ncells_df_long %>% group_by(bin_size, metric, value) %>% 
  summarise(n=n()) %>% 
  group_by(bin_size, metric) %>%
  mutate(prop_n = n/sum(n))

ncells_df_long_summ <- ncells_df_long %>%
  group_by(bin_size, metric, value) %>%
  mutate(n = n()) %>% 
  group_by(bin_size, metric) %>%
  summarise(median = median(value),
            # n of the median value
            n = n[which(value == median)][1])

# Plot of number of cells and number of cell types per bin size
p_cell_hist <- ggplot(ncells_df_long, aes(x=value)) +
  geom_histogram(binwidth=1, color="black", fill="gray70", linewidth=0.25) +
  geom_point(aes(x = median, y = n), size=1, color="black", shape=21, fill="white",
             data = ncells_df_long_summ, 
             inherit.aes = FALSE
  ) +
  ggh4x::facet_grid2(metric~bin_size, scales="free", independent="y",
                     labeller=labeller(metric=c("n_cells"="Cells", "n_celltypes"="Cell types"),
                                       bin_size=c("2" = "2\u03bcm","4" = "4\u03bcm",
                                         "8" = "8\u03bcm","16" = "16\u03bcm","32" = "32\u03bcm"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::comma,
                     name="Number of bins") +
  ggtitle("Distribution of cells and cell types per simulated bin size") +
  theme_classic(base_size=8) +
  theme_histogram

p_cell_hist

if (!make_combined_plot){
  ggsave(paste0(plot_path, "histogram_vizgen_ncells_ncelltypes_per_binsize.pdf"),
         p_cell_hist,
         device = cairo_pdf,
         width=8, height=5, dpi=300)
}
# Stats
ncells_df_long %>% group_by(bin_size, metric) %>% 
  summarise(mean=mean(value), median=median(value), sd=sd(value))

#### MAKE FIGURE 4.1 ####
# Figure with 8µm grid and histograms
p_grids_8um <- p_grids[[3]]

bin_size_oi <- c(8, 16, 32)
p_hist_oi <- ggplot(ncells_df_long %>% filter(bin_size %in% bin_size_oi), aes(x=value)) +
  geom_histogram(binwidth=1, color="black", fill="gray70", linewidth=0.25) +
  geom_point(aes(x = median, y = n), size=1, color="black", shape=21, fill="white",
             data = ncells_df_long_summ %>% filter(bin_size %in% bin_size_oi), 
             inherit.aes = FALSE
  ) +
  ggh4x::facet_grid2(metric~bin_size, scales="free", independent="y",
                     labeller=labeller(metric=c("n_cells"="Cells", "n_celltypes"="Cell types"),
                                       bin_size=paste0(bin_size_oi, "\u03bcm") %>% setNames(bin_size_oi))) +
  ggh4x::facetted_pos_scales(
    y = list(
      bin_size == 8 ~ scale_y_continuous(labels = scales::label_number(scale_cut =scales::cut_short_scale()),
                                         name = "Number of bins", expand = expansion(mult=c(0, 0.05))),
      bin_size == 16 ~ scale_y_continuous(labels = scales::label_number(scale_cut =scales::cut_short_scale()),
                                          expand = expansion(mult=c(0, 0.05))),
      bin_size == 32 ~ scale_y_continuous(labels = scales::comma, expand = expansion(mult=c(0, 0.05)))
    )
  ) +
  ggtitle("(b) Distribution of cells and cell types per simulated bin size") +
  theme_classic(base_size=7) +
  theme_histogram

p_hist_oi

p_combined <- p_grids_8um + ggtitle("(a) Artificial 8\u03bcm grid on Vizgen data") +
  theme(plot.title = element_text(size=7, face="bold", hjust=0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.spacing.y = unit(0.05, "cm"),
        legend.text = element_text(size=6)) +
  # Change cell type label names
  scale_fill_manual(values=color_palette_vizgen, 
                    labels=sort(c("Hepatocytes", "Endothelial cells", "Cholangiocytes",
                                  "Stromal cells", "Kupffer cells", "Other immune cells",
                                  "B cells", "Unknown"))) +
  p_hist_oi

bin_size_oi_strs <- paste0(bin_size_oi, collapse="_")
ggsave(paste0(plot_path, "vizgen_grid_8um_with_histogram_", bin_size_oi_strs, ".pdf"),
       device = cairo_pdf,
       p_combined,
       width=8, height=5, dpi=300)

# Extra figure for 2 and 4µm (not used)
p_grids_supp <- p_grids[c(1,2,4,5)] # 2, 4, 16, 32 um
p_grids_supp_wrap <- wrap_elements(wrap_plots(p_grids_supp, ncol=2, guides="collect") & 
                                     theme(plot.title = element_text(size=6, hjust=0.5),
                                           legend.position = "bottom",
                                           legend.title = element_blank(),
                                           legend.key.size = unit(0.3, "cm"),
                                           legend.key.spacing.x = unit(0.3, "cm"),
                                           legend.key.spacing.y = unit(0.05, "cm"),
                                           legend.text = element_text(size=6)) &
                                     # Change cell type label names
                                     scale_fill_manual(values=color_palette_vizgen, 
                                                       labels=sort(c("Hepatocytes", "Endothelial cells", "Cholangiocytes",
                                                                     "Stromal cells", "Kupffer cells", "Other immune cells",
                                                                     "B cells", "Unknown"))))

p_hist_supp <- ggplot(ncells_df_long %>% filter(!bin_size %in% bin_size_oi), aes(x=value)) +
  geom_histogram(binwidth=1, color="black", fill="gray70", linewidth=0.25) +
  geom_point(aes(x = median, y = n), size=1, color="black", shape=21, fill="white",
             data = ncells_df_long_summ %>% filter(!bin_size %in% bin_size_oi), 
             inherit.aes = FALSE) +
  ggh4x::facet_grid2(metric~bin_size, scales="free", independent="y",
                     labeller=labeller(metric=c("n_cells"="Cells", "n_celltypes"="Cell types"),
                                       bin_size=c("4" = "4\u03bcm",
                                                  "2" = "2\u03bcm"))) +
  scale_y_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale()),
                    name = "Number of bins", expand = expansion(mult=c(0, 0.05))) +
  ggtitle("(b) Distribution of cells and cell types per simulated bin size") +
  theme_classic(base_size=7) +
  theme_histogram

p_hist_supp

p_grids_supp_wrap + ggtitle("(a) Artificial grid on Vizgen data") +
  theme(plot.title = element_text(size=7, face="bold", hjust=0)) +
  p_hist_supp + plot_layout(widths=c(2, 1))

ggsave(paste0(plot_path, "vizgen_grid_supp_with_histogram.pdf"),
       device = cairo_pdf,
       width=8, height=4, dpi=300)


#### PLOT PROPORTIONS ####
bin_size_oi <- c(8, 16, 32)

# Compare proportion based on area vs transcripts
for (bin_size in bin_sizes){
  
  vizgen_transcript <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size, "um_celltype_proportions.csv"), row.names=1)
  vizgen_area <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size, "um_celltype_proportions_by_area.csv"), row.names=1)
  
  stopifnot("all row names must be the same" = all(rownames(vizgen_transcript) == rownames(vizgen_area)))
  stopifnot("all col names must be the same" = all(colnames(vizgen_transcript) == colnames(vizgen_area)))
  
  mean_cor <- mean(diag(cor(as.matrix(vizgen_transcript), as.matrix(vizgen_area), method="pearson")))
  print(paste0("Mean correlation at ", bin_size, "um: ", round(mean_cor, 3)))
}

# We will use proportion based on transcripts so it aligns with deconvolution
# Get the proportions plot
props_df <- lapply(bin_sizes, function(bin_size) {
  read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size, "um_celltype_proportions.csv")) %>%
    pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>%
    mutate(bin_size = bin_size)
}) %>% bind_rows()

# Calculate proportion of hepatocytes
read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size,
                "um_celltype_proportions.csv"),
         row.names = 1) %>% 
  # Get max proportion per grid
  apply(1, which.max) %>% 
  table %>% prop.table %>% 
  setNames(colnames(read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size,
                                    "um_celltype_proportions.csv"),
                                         row.names = 1)))

vizgen_ori_props <- read.csv("visium_hd_liver_combined/vizgen/sdata_polyT_obs.csv") %>% 
  select(cells, annotation_own_score_genes) %>% 
  group_by(annotation_own_score_genes) %>%
  summarise(n_cells = n()) %>% 
  mutate(proportion = n_cells / sum(n_cells),
         bin_size = 0)

props_df_summ <- props_df %>%
  group_by(bin_size, celltype) %>% 
  summarise(mean_prop=mean(proportion)) %>% 
  bind_rows(vizgen_ori_props %>% rename(celltype = annotation_own_score_genes) %>%
              select(bin_size, celltype, mean_prop=proportion)) %>% 
  # remove space and dot from name of celltype
  mutate(celltype = gsub("[ .]", "", celltype)) %>% 
  mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype))

celltype_order <- props_df_summ %>%
  group_by(celltype) %>%
  summarise(overall_mean = mean(mean_prop)) %>%
  arrange(-overall_mean) %>%
  pull(celltype)

props_df_summ <- props_df_summ %>% 
  # factor celltype based on mean proportion
  mutate(celltype = factor(celltype, levels = celltype_order),
         bin_size = factor(bin_size, levels=c(0, bin_sizes),
                           labels=c("Original Dataset", paste0(bin_sizes, "\u03bcm"))))

if (!make_combined_plot){
  # Plot barplot
  p_props_bar <- ggplot(props_df_summ, aes(x=celltype, y=mean_prop, fill=celltype)) +
    geom_bar(stat="identity") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x="Cell Type", y="Mean Proportion") +
    scale_fill_manual(values=color_palette_vizgen, name="Cell type", 
                      labels=c("Hepatocytes", "Endothelial cells", "Unknown",
                               "Stromal cells", "Kupffer cells", "Other immune cells",
                               "Cholangiocytes", "B cells")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(~bin_size, nrow=1) +
    # one legend legend
    guides(fill=guide_legend(nrow=1, title.position = "top")) +
    labs(title = "Mean Cell Type Proportions in Gridded Vizgen Data") +
    theme_minimal(base_size=8) +
    theme_barplot

  ggsave(paste0(plot_path, "barplot_vizgen_meanprops_per_binsize.pdf"),
         p_props_bar,
         device = cairo_pdf,
         width=10, height=3, dpi=300)
} else {
  p_props_bar <- ggplot(props_df_summ %>%
                          filter(bin_size %in% c(paste0(bin_size_oi, "\u03bcm"), "Original Dataset")),
                        aes(x=bin_size, y=mean_prop, fill=celltype)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x="Cell Type", y="Mean Proportion") +
    scale_fill_manual(values=color_palette_vizgen, name="Cell type", 
                      labels=c("Hepatocytes", "Endothelial cells", "Unknown",
                               "Stromal cells", "Kupffer cells", "Other immune cells",
                               "Cholangiocytes", "B cells")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    guides(fill=guide_legend(nrow=1, title.position = "top")) +
    labs(title = "Mean Cell Type Proportions in Gridded Vizgen Data") +
    theme_minimal(base_size=8) +
    theme_barplot +
    theme(axis.text.x = element_text(size=7))
}

#### SPATIAL BARPLOT ####
doubletize <- FALSE
vizgen_df <- lapply(bin_sizes, function(bin_size) {
  read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size, "um_celltype_proportions.csv")) %>%
    pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>%
    mutate(celltype = gsub("\\.", "", celltype)) %>% 
    mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype)) %>% 
    mutate(bin_size = bin_size)
}) %>% bind_rows() %>% filter(proportion > 0)

bin_sizes_oi <- c(8, 16, 32)
p_spatial_barplots <- lapply(bin_sizes_oi, function(bs) {
  centroids <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid",
                               bs, "um_centroids.csv")) %>% 
    rename(x = centroid_x, y = centroid_y) %>% 
    filter(x >= crd[1], x <= crd[2],
           y >= crd[3], y <= crd[4])
  
  square_size <- bs*1000/108
  vizgen_df_coords <- vizgen_df %>%
    filter(bin_size == bs) %>%
    ungroup() %>% 
    inner_join(centroids, by = "grid_id") %>% 
    rename(spot = grid_id, x=y, y=x)
  
  plot_scatterbar(vizgen_df_coords,
                  ncelltypes = ifelse(doubletize, 2, Inf),
                  square_size = square_size,
                  color_palette = color_palette_vizgen)
})
p_spatial_barplots

#### VIZGEN TRANSCRIPT COUNTS BOXPLOT ####

vizgen_metadata_df <- read.csv("visium_hd_liver_combined/vizgen/sdata_polyT_obs.csv") %>% 
  mutate(cells = as.character(cells)) %>% 
  rename(celltype = annotation_own_score_genes) %>%
  mutate(celltype = gsub(" ", "", celltype)) %>%
  mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype)) %>% 
  mutate(celltype = factor(celltype, levels=celltype_order))

# Calculate average transcripts
vizgen_metadata_df %>% group_by(celltype) %>% 
  summarise(mean_counts = mean(total_counts)) %>% 
  # How much higher is Hepatocytes than the rest
  mutate(ratio_to_hepato = mean_counts[celltype == "Hepatocytes"]/mean_counts) %>% 
  filter(!celltype %in% c("Hepatocytes", "Unknown")) %>% 
  pull(ratio_to_hepato) %>% 
  mean()

# Make two plots - size vs transcripts per cell
vizgen_meta_plots <- lapply(c("shapeSize", "total_counts"), function(metric){
  titles <- c("shapeSize" = "Total cell size (pixels)",
              "total_counts" = "Total transcripts per cell")
  ggplot(vizgen_metadata_df,
         aes(x=celltype, y=!!sym(metric), fill=celltype, color=celltype)) +
    geom_boxplot(show.legend = FALSE, outlier.shape=NA, linewidth=0.15) +
      scale_fill_manual(values = color_palette_vizgen, name = "Cell type",
                        labels=add_space_celltype) +
      scale_colour_manual(values = c(rep("black", 2), "white", rep("black", 5))) + 
    geom_boxplot(color="black", fill=NA, fatten=NULL, linewidth=0.15,
      outlier.size = 0.4, outlier.shape = 16, outlier.stroke = 0) +
    scale_y_continuous(labels = scales::comma) +
    guides(fill=guide_legend(nrow=1, title.position = "top")) +
    labs(y=titles[metric]) +
    theme_minimal(base_size=7) +
    theme_barplot
})
vizgen_meta_plots

# Finally make Figure 4.4a-c
p_grids_wrapped <- wrap_elements(p_grids[[3]] + p_grids[[4]] + p_grids[[5]] + plot_layout(ncol=1) &
                                   theme(plot.title = element_blank()))
p_spatialbarplots_wrapped <- wrap_elements(wrap_plots(p_spatial_barplots) + plot_layout(ncol=1) &
                                             theme(legend.position = "none"))
p_spatial <- wrap_elements(p_grids_wrapped + p_spatialbarplots_wrapped + plot_layout(ncol=2, widths=c(1,1.2)))

p_vizgen_meta <- wrap_elements(wrap_plots(vizgen_meta_plots) + plot_layout(nrow=1))
design <- "AABB
           AACD
           AAEE"

combined_plots <- #free(p_grids_wrapped) + (p_cell_hist / p_props_bar) +
  p_spatial +
  (p_props_bar + theme(legend.position = "none")) +
  vizgen_meta_plots[[2]] + vizgen_meta_plots[[1]] +
    plot_spacer() +
  plot_layout(design = design)
combined_plots
ggsave(paste0(plot_path, "vizgen_grid_with_props_UMIs",
              ifelse(doubletize, "_doubletize", ""), ".pdf"),
       combined_plots,
       device = cairo_pdf,
       width=8, height=6, dpi=300)