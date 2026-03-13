source("visium_hd_liver_combined/0_utils.R")
library(spdep)
library(sf)
library(ggnewscale)

# Figure 5.1
#### Vizgen - plot segmented cells with grid overlay ####
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

crd <- c(10666, 12445, 9776, 11260)
crd_poly <-  st_multipoint(
  matrix(c(crd[1], -crd[3],
           crd[2], -crd[3],
           crd[2], -crd[4],
           crd[1], -crd[4],
           crd[1], -crd[3]), ncol=2, byrow=TRUE)
) %>% st_cast("POLYGON")

bin_size <- 2
grid <- read_sf(paste0("visium_hd_liver_combined/vizgen/grid_", bin_size, "um_roi.shp"))
grid_roi <- st_filter(grid, crd_poly, .predicate = st_within)
lims <- st_bbox(grid_roi)

p_grid <- ggplot() +
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

p_grid

#### VisiumHD - plot segmented nuclei ####
# Read in shape file
segmented_nuclei <- read_sf("visium_hd_liver/rds/segmentation_mask_boundaries_roi.shp")

grid <- read_sf(paste0("visium_hd_liver/rds/grid_2um_roi.shp"))
# grid_roi <- st_filter(grid, crd_poly, .predicate = st_within)
#   lims <- st_bbox(grid_roi)

p_grid_nuclei <- ggplot() +
  geom_sf(fill="#8DD3C7", color = "white", linewidth=0.1, data = segmented_nuclei) +
  geom_sf(fill=NA, color="gray40", linewidth=0.05, data = grid) +
  theme_void() +
  coord_sf(xlim = c(14000, 15000),
           ylim = c(-13500, -12500),
           expand = FALSE)

combined_grids <- p_grids[[1]] + p_grid_nuclei
combined_grids
ggsave(paste0(plot_path, "2um_vizgen_grid_with_segmented_nuclei.pdf"),
       combined_grids,
       device = cairo_pdf,
       width=8, height=4, dpi=300)


#### HISTOGRAMS ####
#### Vizgen - cell and cell types per bin ####
# Get rds from 2_2_plot_vizgen_gridded
ncells_df <- readRDS("visium_hd_liver_combined/vizgen/ncells_ncelltypes_per_grid.rds")

ncells_df_long <- ncells_df %>%
  filter(bin_size == bin_size) %>% 
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
  ggh4x::facet_wrap2(~metric, ncol=1,
                     labeller=labeller(metric=c("n_cells"="Number of cells", "n_celltypes"="Number of cell types"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::label_number(scale_cut=scales::cut_short_scale()),
                     breaks = seq(0, 3e06, by=5e05),
                     name="Frequency") +
  theme_classic(base_size=8) +
  theme(strip.background = element_rect(linewidth = 0.5),
        strip.text.y.right = element_text(size=7),
        strip.text.x.top = element_text(hjust=0.5, size=7),
        panel.spacing.y = unit(5, "mm"),
        panel.grid.major.y = element_line(color="gray90", linewidth=0.2),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, color="black",linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6, margin = margin(r=7)),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(size=5),
        axis.line.y.left = element_blank(),
        axis.line.x.bottom = element_line(linewidth=0.1),
        axis.ticks = element_line(linewidth=0.1))

p_cell_hist

p_grids[[1]] / p_cell_hist +
  plot_layout(heights = c(3,2))

#### VisiumHD - Bin per cells, cells per bin ####
binned_nuclei_df <- bind_rows(read.csv("visium_hd_liver/rds/grids_per_nuclei.csv") %>% 
  rename(n = X) %>% mutate(type = "n_grids"),
  read.csv("visium_hd_liver/rds/nuclei_per_2um_grid.csv") %>% 
  rename(n = X) %>% mutate(type = "n_nuclei")) %>% 
  mutate(type = factor(type, levels=c("n_nuclei", "n_grids")))

# Plot of number of cells and number of cell types per bin size
p_nuclei_hist <- ggplot(binned_nuclei_df, aes(x=n, y=count)) +
  geom_bar(color="black", fill="gray70", linewidth=0.25, stat="identity") +
  ggh4x::facet_wrap2(~type, scales="free", ncol=1,
                     labeller=labeller(type=c("n_grids"="Number of 2µm bins per nucleus",
                                                "n_nuclei"="Number of nuclei per 2µm bin"))) +
  ggh4x::facetted_pos_scales(
    y = list(
      type == "n_nuclei" ~ scale_y_continuous(labels = scales::label_number(scale_cut =scales::cut_short_scale()),
                                             name = "Frequency", expand = expansion(mult=c(0, 0.05))),
      type == "n_grids" ~ scale_y_continuous(labels = scales::comma, expand = expansion(mult=c(0, 0.05)))
    )) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  theme_classic(base_size=8) +
  theme(strip.background = element_rect(linewidth = 0.5),
        strip.text.y.right = element_text(size=7),
        strip.text.x.top = element_text(hjust=0.5, size=7),
        panel.spacing.y = unit(5, "mm"),
        panel.grid.major.y = element_line(color="gray90", linewidth=0.2),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, color="black",linewidth = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6, margin = margin(r=7)),
        axis.text.y = element_text(size=5),
        axis.text.x = element_text(size=5),
        axis.line.y.left = element_blank(),
        axis.line.x.bottom = element_line(linewidth=0.1),
        axis.ticks = element_line(linewidth=0.1))

p_nuclei_hist

binned_nuclei_df %>% filter(type == "n_grids")

expanded <- rep(binned_nuclei_df %>% filter(type == "n_grids") %>% pull(n), binned_nuclei_df %>% filter(type == "n_grids") %>% pull(count))
mean(expanded)
median(expanded)

binned_nuclei_df %>% filter(type == "n_nuclei") %>% 
  mutate(prop = count / sum(count))

combined_plots <- p_cell_hist / p_nuclei_hist

ggsave(paste0(plot_path, "2um_vizgen_grid_with_segmented_nuclei_histogram.pdf"),
       p_cell_hist / p_nuclei_hist,
       device = cairo_pdf,
       width=3, height=7, dpi=300)
