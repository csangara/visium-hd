library(Seurat)
library(tidyverse)
library(patchwork)
library(spdep)
library(sf)

color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "Endothelialcells" = "#dca754",
                   "Cholangiocytes" = "#C61B84FF",
                   "Stromalcells" = "#79151e",
                   "Kupffercells" = "#5DA6DBFF",
                   "Otherimmunecells" = "#893A86FF",
                   "Bcells" = "#9C7EBAFF",
                   "Unknown" = "#191919")
add_space_celltype <- function(celltype_strs) {
  celltype_strs %>%
  # Add space before the word "cells"
  stringr::str_replace_all("([A-Za-z]+)(cells)$", "\\1 \\2") %>%
  # Add space before the word "immune"
  stringr::str_replace_all("([A-Za-z]+)(immune)", "\\1 \\2")
}
  
bin_sizes <- c(2, 4, 8,16,32)
bin_size_strs <- sprintf("%03dum", bin_sizes)
bin_size <- 8
make_combined_plot <- TRUE
first_run <- FALSE

#### SEGMENTED CELLS AND GRID ####
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
  #crd <- c(11000, 14000, 9700, 11000)
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
    scale_fill_manual(values=color_palette) +
    labs(title=paste0(bin_size, "\u03bcm")) +
    theme(legend.position = "none",
          plot.title = element_text(hjust=0.5, size=8))
})

wrap_plots(p_grids, ncol=1) & theme(plot.title = element_blank())

if (!make_combined_plot){
  p_grids_wrapped <- wrap_plots(p_grids) &
    theme(plot.margin = margin(l=5, r=5))
  ggsave("visium_hd_liver_combined/plots/segmented_cells_with_grid_overlay.pdf",
         p_grids_wrapped,
         device = cairo_pdf,
         width=length(p_grids)*3, height=3, dpi=300)
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
                                       bin_size=c("8" = "8\u03bcm",
                                                  "16" = "16\u03bcm",
                                                  "32" = "32\u03bcm"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::comma,
                     name="Number of bins") +
  theme_classic(base_size=8) +
  ggtitle("Distribution of cells and cell types per simulated bin size") +
  theme(strip.background = element_rect(linewidth = 0.5),
        strip.text.y.right = element_text(size=7),
        strip.text.x.top = element_text(hjust=0.5, size=7),
        panel.spacing.y = unit(5, "mm"),
        panel.grid.major.x = element_line(color="gray90", linewidth=0.1),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, color="black",linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6, margin = margin(r=7)),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(size=5),
        axis.line.y.left = element_line(linewidth=0.25),
        axis.line.x.bottom = element_line(linewidth=0.1),
        axis.ticks = element_line(linewidth=0.1),
        plot.title = element_text(size=7, face="bold"))

if (!make_combined_plot){
  ggsave("visium_hd_liver_combined/plots/histogram_vizgen_ncells_ncelltypes_per_binsize.pdf",
         p_cell_hist,
         device = cairo_pdf,
         width=8, height=5, dpi=300)
}
# Stats
ncells_df_long %>% group_by(bin_size, metric) %>% 
  summarise(mean=mean(value), median=median(value), sd=sd(value))

#### PLOT PROPORTIONS ####
# Get the proportions plot
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
    theme_minimal(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x="Cell Type", y="Mean Proportion") +
    scale_fill_manual(values=color_palette, name="Cell type", 
                      labels=c("Hepatocytes", "Endothelial cells", "Unknown",
                               "Stromal cells", "Kupffer cells", "Other immune cells",
                               "Cholangiocytes", "B cells")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(~bin_size, nrow=1) +
    # one legend legend
    guides(fill=guide_legend(nrow=1, title.position = "top")) +
    labs(title = "Mean Cell Type Proportions in Gridded Vizgen Data") +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(hjust=0.5, size=6),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=6, margin = margin(r=7)),
          plot.title = element_text(size=7, face="bold"),
          legend.position = "bottom",
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
          legend.key.size = unit(0.4, "cm"),
          legend.spacing.x = unit(0.4, "cm"),
          legend.key.spacing.x = unit(0.3, "cm"),
          legend.title = element_text(size=6, hjust=0.5, margin=margin(b=5)),
          legend.text = element_text(size=6, margin = margin(l=3)))
  
  ggsave("visium_hd_liver_combined/plots/barplot_vizgen_meanprops_per_binsize.pdf",
         p_props_bar,
         device = cairo_pdf,
         width=8, height=3, dpi=300)
} else {
  p_props_bar <- ggplot(props_df_summ, aes(x=bin_size, y=mean_prop, fill=celltype)) +
    geom_bar(stat="identity", position="dodge") +
    theme_minimal(base_size=8) +
    labs(x="Cell Type", y="Mean Proportion") +
    scale_fill_manual(values=color_palette, name="Cell type", 
                      labels=c("Hepatocytes", "Endothelial cells", "Unknown",
                               "Stromal cells", "Kupffer cells", "Other immune cells",
                               "Cholangiocytes", "B cells")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    guides(fill=guide_legend(nrow=1, title.position = "top")) +
    labs(title = "Mean Cell Type Proportions in Gridded Vizgen Data") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.text.x = element_text(size=7),
          axis.text.y = element_text(size=5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=6),
          plot.title = element_text(size=7, face="bold"),
          legend.position = "bottom",
          legend.box.margin = margin(t = 0, r = 200, b = 0, l = 0),
          legend.key.size = unit(0.4, "cm"),
          legend.spacing.x = unit(0.4, "cm"),
          legend.key.spacing.x = unit(0.3, "cm"),
          legend.title = element_text(size=6, hjust=0.5, margin=margin(b=5)),
          legend.text = element_text(size=6, margin = margin(l=3)))
}

# Wrap element of p_grids
p_grids_wrapped <- wrap_elements(wrap_plots(p_grids) + plot_layout(ncol=1) &
                                   theme(plot.title = element_blank()))
                                   
p_spatialbarplots_wrapped <- wrap_elements(wrap_plots(p_spatial_barplots) + plot_layout(ncol=1))

p_spatial <- wrap_elements(p_grids_wrapped + p_spatialbarplots_wrapped + plot_layout(ncol=2, widths=c(1,1.2)))

design <- "AB
           AB
           AC"

combined_plots <- #free(p_grids_wrapped) + (p_cell_hist / p_props_bar) +
  p_spatial + p_cell_hist + p_props_bar +
  plot_layout(design = design) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')
combined_plots
ggsave("visium_hd_liver_combined/plots/vizgen_grid_analysis.pdf",
       combined_plots,
       device = cairo_pdf,
       width=8, height=6, dpi=300)

##### MAKE FIGURE WITH 8 MICRON GRID AND GRIDSIZE 2,8, 16 MICRON #####
# MAIN FIGURE
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
  # Works for bin_size_oi (2, 8, 16) 
  #scale_y_continuous(labels = scales::label_number(scale_cut=append(scales::cut_short_scale())),
   #                   name = "Number of bins", expand = expansion(mult=c(0, 0.05))) 
  ggh4x::facetted_pos_scales(
    y = list(
      bin_size == 8 ~ scale_y_continuous(labels = scales::label_number(scale_cut =scales::cut_short_scale()),
                                         name = "Number of bins", expand = expansion(mult=c(0, 0.05))),
      bin_size == 16 ~ scale_y_continuous(labels = scales::label_number(scale_cut =scales::cut_short_scale()),
                                          expand = expansion(mult=c(0, 0.05))),
      bin_size == 32 ~ scale_y_continuous(labels = scales::comma, expand = expansion(mult=c(0, 0.05)))
    )
  ) +
  theme_classic(base_size=7) +
  ggtitle("(b) Distribution of cells and cell types per simulated bin size") +
  theme(strip.background = element_rect(linewidth = 0.5),
        strip.text.y.right = element_text(size=7),
        strip.text.x.top = element_text(hjust=0.5, size=7),
        panel.spacing.y = unit(5, "mm"),
        panel.grid.major.x = element_line(color="gray90", linewidth=0.1),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, color="black",linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6, margin = margin(r=7)),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(size=5),
        axis.line.y.left = element_line(linewidth=0.25),
        axis.line.x.bottom = element_line(linewidth=0.1),
        axis.ticks = element_line(linewidth=0.1),
        plot.title = element_text(size=7, face="bold"))

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
  scale_fill_manual(values=color_palette, 
                    labels=sort(c("Hepatocytes", "Endothelial cells", "Cholangiocytes",
                             "Stromal cells", "Kupffer cells", "Other immune cells",
                             "B cells", "Unknown"))) +
  p_hist_oi

bin_size_oi_strs <- paste0(bin_size_oi, collapse="_")
ggsave(paste0("visium_hd_liver_combined/plots/vizgen_grid_8um_with_histogram_", bin_size_oi_strs, ".pdf"),
       device = cairo_pdf,
       p_combined,
       width=8, height=5, dpi=300)

# Supplementary figure
p_grids_supp <- p_grids[c(1,2,4,5)] # 2, 4, 16, 32 um +
p_grids_supp_wrap <- wrap_elements(wrap_plots(p_grids_supp, ncol=2, guides="collect") & 
  theme(plot.title = element_text(size=6, hjust=0.5),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.key.spacing.y = unit(0.05, "cm"),
        legend.text = element_text(size=6)) &
  # Change cell type label names
  scale_fill_manual(values=color_palette, 
                    labels=sort(c("Hepatocytes", "Endothelial cells", "Cholangiocytes",
                                  "Stromal cells", "Kupffer cells", "Other immune cells",
                                  "B cells", "Unknown"))))

p_hist_supp <- ggplot(ncells_df_long %>% filter(!bin_size %in% bin_size_oi), aes(x=value)) +
  geom_histogram(binwidth=1, color="black", fill="gray70", linewidth=0.25) +
  geom_point(aes(x = median, y = n), size=1, color="black", shape=21, fill="white",
             data = ncells_df_long_summ %>% filter(!bin_size %in% bin_size_oi), 
             inherit.aes = FALSE
  ) +
  ggh4x::facet_grid2(metric~bin_size, scales="free", independent="y",
                     labeller=labeller(metric=c("n_cells"="Cells", "n_celltypes"="Cell types"),
                                       bin_size=c("4" = "4\u03bcm",
                                                  "32" = "32\u03bcm")))  +
  ggh4x::facetted_pos_scales(
    y = list(
      bin_size == 4 ~ scale_y_continuous(labels = scales::label_number(scale_cut =scales::cut_short_scale()),
                                         name = "Number of bins", expand = expansion(mult=c(0, 0.05))),
      bin_size == 32 ~ scale_y_continuous(labels = scales::comma, expand = expansion(mult=c(0, 0.05)))
    )
  ) +
  theme_classic(base_size=7) +
  ggtitle("(b) Distribution of cells and cell types per simulated bin size") +
  theme(strip.background = element_rect(linewidth = 0.5),
        strip.text.y.right = element_text(size=7),
        strip.text.x.top = element_text(hjust=0.5, size=7),
        panel.spacing.y = unit(5, "mm"),
        panel.grid.major.x = element_line(color="gray90", linewidth=0.1),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, color="black",linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6, margin = margin(r=7)),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(size=5),
        axis.line.y.left = element_line(linewidth=0.25),
        axis.line.x.bottom = element_line(linewidth=0.1),
        axis.ticks = element_line(linewidth=0.1),
        plot.title = element_text(size=7, face="bold"))

p_hist_supp

p_grids_supp_wrap + ggtitle("(a) Artificial grid on Vizgen data") +
  theme(plot.title = element_text(size=7, face="bold", hjust=0)) +
  p_hist_supp + plot_layout(widths=c(2, 1))
ggsave("visium_hd_liver_combined/plots/vizgen_grid_supp_with_histogram.pdf",
       device = cairo_pdf,
       width=8, height=4, dpi=300)

#### SPATIAL BARPLOT ####
doubletize <- FALSE
if (doubletize){
  if (first_run) {
    vizgen_df <- lapply(bin_sizes, function(bs) {
      vizgen <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_",
                                bs, "um_celltype_proportions.csv"))
      
      vizgen %>%
        pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>% 
        mutate(celltype = gsub("\\.", "", celltype)) %>% 
        mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype)) %>% 
        filter(proportion > 0) %>%
        group_by(grid_id) %>%
        # Only get top two
        slice_max(order_by=proportion, n=2) %>%
        mutate(rank = data.table::frankv(proportion, order=-1),
               doublet_type = case_when(n() == 1 ~ "singlet",
                                        n() == 2 ~ "doublet")) %>% 
        # Recalculate proportions to become 1
        group_by(grid_id) %>%
        mutate(proportion = proportion / sum(proportion),
               bin_size = bs) %>%
        ungroup()
    }) %>% bind_rows()
    
    # Save rds
    saveRDS(vizgen_df, "visium_hd_liver_combined/vizgen/vizgen_props_doublet.rds")
  } else {
    vizgen_df <- readRDS("visium_hd_liver_combined/vizgen/vizgen_props_doublet.rds")
  }
} else {
  vizgen_df <- lapply(bin_sizes, function(bin_size) {
    read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size, "um_celltype_proportions.csv")) %>%
      pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>%
      mutate(celltype = gsub("\\.", "", celltype)) %>% 
      mutate(celltype = ifelse(celltype == "Kuppfercells", "Kupffercells", celltype)) %>% 
      mutate(bin_size = bin_size)
  }) %>% bind_rows() %>% 
    filter(proportion > 0)
}


bin_sizes <- c(8, 16, 32)
p_spatial_barplots <- lapply(bin_sizes, function(bs) {
  centroids <- read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid",
                               bs, "um_centroids.csv"))
  square_size <- bs*1000/108
  
  # merge with coordinates
  vizgen_df_coords <- vizgen_df %>% filter(bin_size==bs) %>% 
    ungroup() %>% 
    left_join(centroids, by = "grid_id") %>% 
    rename(x = centroid_x, y = centroid_y) %>% 
    mutate(n_celltypes = n(), .by = "grid_id",
           square_size = square_size,
           bin_size = bin_size)
  
  vizgen_df_coords <- vizgen_df_coords %>% 
    filter(x >= crd[1], x <= crd[2],
           y >= crd[3], y <= crd[4])
  
  vizgen_df_square <- vizgen_df_coords %>%
    filter(n_celltypes == 1) %>% 
    mutate(group = as.character(grid_id)) %>% 
    # Draw squares for each
    mutate(x1 = x - square_size / 2,
           y1 = y - square_size / 2,
           x2 = x + square_size / 2,
           y2 = y - square_size / 2,
           x3 = x + square_size / 2,
           y3 = y + square_size / 2,
           x4 = x - square_size / 2,
           y4 = y + square_size / 2) %>%
    rename(coord_x = x, coord_y = y) %>% 
    pivot_longer(cols = c(x1, y1, x2, y2, x3, y3, x4, y4),
                 names_to = c(".value", "corner"),
                 names_pattern = "(x|y)([1-4])")
  
  # Create barplot in squares
  vizgen_df_barplot <- vizgen_df_coords %>% 
    group_by(grid_id) %>% arrange(grid_id, desc(proportion)) %>% 
    mutate(rank = row_number(),
           group = paste0(grid_id, "_", rank),
           cum_prop = cumsum(proportion)) %>% 
    mutate(x1 = x - (square_size/2) + (cum_prop*square_size),
           y1 = y - square_size / 2,
           x2 = x - (square_size/2) + (cum_prop*square_size),
           y2 = y + square_size / 2,
           x3 = x - (square_size/2) + ((cum_prop - proportion)*square_size),
           y3 = y + square_size / 2,
           x4 = x - (square_size/2) + ((cum_prop - proportion)*square_size),
           y4 = y - square_size / 2
    ) %>% 
    rename(coord_x = x, coord_y = y) %>%
    pivot_longer(cols = c(x1, y1, x2, y2, x3, y3, x4, y4),
                 names_to = c(".value", "corner"),
                 names_pattern = "(x|y)([1-4])")
  
  vizgen_df_shapes <- bind_rows(vizgen_df_square, vizgen_df_barplot)
  
  ggplot(vizgen_df_shapes, aes(x = x, y = y)) +
    geom_polygon(aes(fill = celltype, group = group), show.legend = TRUE) +
    # White border
    geom_tile(data = vizgen_df_coords %>% distinct(x, y),
              aes(x = x, y = y), height = square_size, width = square_size,
              fill = NA, color = "white", inherit.aes = FALSE) +
    scale_fill_manual(values = color_palette) +
    theme_void() +
    coord_fixed() +
    scale_y_reverse() +
    theme(legend.position="none")
})

# Sanity check
wrap_elements(wrap_plots(p_grids, ncol=1) & theme(plot.title = element_blank())) + 
  wrap_elements(wrap_plots(p_spatial_barplots, ncol=1))

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

# Make two plots
vizgen_meta_plots <- lapply(c("shapeSize", "total_counts"), function(metric){
  titles <- c("shapeSize" = "Total cell size (pixels)",
              "total_counts" = "Total transcripts per cell")
  ggplot(vizgen_metadata_df,
         aes(x=celltype, y=!!sym(metric), fill=celltype, color=celltype)) +
    geom_boxplot(show.legend = FALSE, outlier.shape=NA, linewidth=0.15) +
      scale_fill_manual(values = color_palette, name = "Cell type",
                        labels=add_space_celltype) +
      scale_colour_manual(values = c(rep("black", 2), "white", rep("black", 5))) + 
    #geom_boxplot(color = "black", fill = NA, fatten = NULL, linewidth=0.15, outlier.size=0.1) +
    geom_boxplot(
      color="black",
      fill=NA,
      fatten=NULL,
      #width=0.4,
      linewidth=0.15,
      outlier.size = 0.4,
      outlier.shape = 16,
      outlier.stroke = 0
    ) +
  theme_minimal(base_size=7) +
    scale_y_continuous(labels = scales::comma) +
    guides(fill=guide_legend(nrow=1, title.position = "top")) +
    labs(y=titles[metric]) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth=0.25),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=6),
          plot.title = element_text(size=7, face="bold"),
          legend.position = "bottom",
          #legend.box.margin = margin(t = 0, r= =, b = 0, l = 0),
          legend.key.size = unit(0.4, "cm"),
          legend.spacing.x = unit(0.4, "cm"),
          legend.key.spacing.x = unit(0.3, "cm"),
          legend.title = element_text(size=6, hjust=0.5, margin=margin(b=5)),
          legend.text = element_text(size=6, margin = margin(l=3)))
  
  
})
vizgen_meta_plots


# Wrap element of p_grids
p_grids_wrapped <- wrap_elements(p_grids[[3]] + p_grids[[4]] + p_grids[[5]] + plot_layout(ncol=1) &
                                   theme(plot.title = element_blank()))

p_spatialbarplots_wrapped <- wrap_elements(wrap_plots(p_spatial_barplots) + plot_layout(ncol=1))

p_spatial <- wrap_elements(p_grids_wrapped + p_spatialbarplots_wrapped + plot_layout(ncol=2, widths=c(1,1.2)))

p_vizgen_meta <- wrap_elements(wrap_plots(vizgen_meta_plots) + plot_layout(nrow=1)
                              )

design <- "AABB
           AACD
           AAEE"

combined_plots <- #free(p_grids_wrapped) + (p_cell_hist / p_props_bar) +
  p_spatial +
  (p_props_bar + theme(legend.position = "none")) +
  vizgen_meta_plots[[2]] + vizgen_meta_plots[[1]] +
    p_atlas_counts +
  plot_layout(design = design)
combined_plots
ggsave("visium_hd_liver_combined/plots/vizgen_grid_with_props_UMIs.pdf",
       combined_plots,
       device = cairo_pdf,
       width=8, height=6, dpi=300)



  