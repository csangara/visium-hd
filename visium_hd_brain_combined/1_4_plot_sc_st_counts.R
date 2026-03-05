library(Seurat)
library(tidyverse)
library(patchwork)
library(ggnewscale)

sc <- readRDS("data/scref_MouseBrain_ABA/WMB-10Xv3_subset.rds")

sc$nCounts <- Matrix::colSums(sc@assays$RNA@counts)
sc$nFeatures <- Matrix::colSums(sc@assays$RNA@counts > 0)

sc_df <- sc@meta.data %>% 
  select(class, nCounts, nFeatures) %>% 
  pivot_longer(cols = c(nCounts, nFeatures),
               names_to = "metric",
               values_to = "value")

color_palette_df <- sc@meta.data %>% 
  distinct(class, class_color) %>% 
  arrange(class)
color_palette <- as.character(color_palette_df$class_color) %>% 
  set_names(as.character(color_palette_df$class))

color_palette

p_sc <- sc_df %>%
  mutate(class = factor(class, levels=rev(color_palette_df$class))) %>% 
  ggplot(aes(y=class, x=value, fill=class)) +
  geom_boxplot(
    linewidth=0.15,
    outlier.size = 0.5,
    outlier.shape = 16,
    outlier.stroke = 0,
    show.legend = FALSE
  ) +
  stat_summary(fun=mean, geom="point", shape=23, size=0.5, stroke=0.25,
               color="black", fill="white") +
  facet_wrap(~metric,
             scales="free_x",
             labeller = labeller(metric = c(nCounts = "Count",
                                         nFeatures = "Feature"))
  ) +
  theme_minimal(base_size=7) +
  scale_fill_manual(values=color_palette,
                    name="Cell type") +
  scale_x_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth=0.25),
        strip.text = element_text(size=7, face="bold"),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

ggsave("visium_hd_brain_combined/plots/boxplot_brain_atlas_counts_per_celltype.pdf",
       device = cairo_pdf,
       plot = p_sc,
       width=7, height=5, dpi=300)


### HISTOGRAM OF VIS HD ###
data_path <- "data/Visium_HD_MouseBrain"

bin_sizes <- c(8)
counts_features_df <- lapply(bin_sizes, function(bin_size) {
  bin_size_str <- sprintf("%03dum", bin_size)
  
  visium_obj_ffpe <- readRDS(paste0(data_path, "/Visium_HD_MouseBrain_", 
                                   bin_size_str, ".rds"))
  visium_obj_ff <- readRDS(paste0(data_path, "_FF/Visium_HD_MouseBrain_FF_",
                                   bin_size_str, ".rds"))
  
  bind_rows(
    data.frame(x = visium_obj_ffpe[[paste0("nCount_Spatial.", bin_size_str), drop=TRUE]], 
               dataset = "ffpe", type = "Count"),
    data.frame(x = visium_obj_ffpe[[paste0("nFeature_Spatial.", bin_size_str), drop=TRUE]], 
               dataset = "ffpe", type = "Feature"),
    data.frame(x = visium_obj_ff[[paste0("nCount_Spatial.", bin_size_str), drop=TRUE]],
               dataset = "fresh_frozen", type = "Count"),
    data.frame(x = visium_obj_ff[[paste0("nFeature_Spatial.", bin_size_str), drop=TRUE]],
               dataset = "fresh_frozen", type = "Feature")
  ) %>% mutate(
    bin_size = paste0(bin_size, "um")
  )
}) %>% bind_rows()

set1_colors <- RColorBrewer::brewer.pal(8, "Set1")[1:2]

# Factor 8, 16, 32
counts_features_df <- counts_features_df %>% 
  mutate(dataset = factor(dataset, levels = c("ffpe", "fresh_frozen"),
                          labels = c("FFPE", "Fresh Frozen")))

means_df <- 
  counts_features_df %>%
  group_by(type, bin_size) %>%
  mutate(max_x = max(x)) %>% 
  group_by(type, bin_size, dataset) %>%
  mutate(text_y = max(table(cut_width(x, width = unique(max_x)/50)))) %>%
  summarise(mean_x = prettyNum(round(mean(x)), big.mark=",", scientific=FALSE),
            n = prettyNum(n(), big.mark=",", scientific=FALSE),
            text_x = mean(x) + unique(max_x)/50*4, text_y = max(text_y))

# Compare counts of the two datasets
p_st <- ggplot(counts_features_df, aes(x = x, fill=dataset)) +
  geom_histogram(color="white", alpha=0.5, bins=50, position="identity", size=0.1) +
  geom_text(aes(x = text_x, y = text_y, color=dataset,
                label = paste0("Mean: ", mean_x, "\nn = ", n)),
            data = means_df, vjust = 1, hjust=0, size=2,
            show.legend = FALSE
  ) +
  facet_wrap(~type, scales = "free",
             labeller = labeller( .multi_line=FALSE)) +
  scale_fill_manual(values = set1_colors, name = "Dataset") +
  scale_color_manual(values = set1_colors) +
  # Use scientific notation for y
  scale_y_continuous(labels = scales::label_number(scale_cut=scales::cut_short_scale())) +
  scale_x_continuous(labels = scales::comma) +
  theme_classic(base_size = 7) +
  theme(legend.position = "bottom",
        # Add right spacing to legend title
        #legend.title = element_text(margin = margin(r = 15), size=8),
        #legend.text = element_text(size = 6),
        #axis.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
  )
ggsave("visium_hd_brain_combined/plots/histogram_brain_counts_features.pdf",
       plot = p_st,
       width=7, height=3, dpi=300)

# Plot count density
data_path <- paste0("data/Visium_HD_MouseBrain")
breaks_list <- list(fresh_frozen = c(100, 1000, 2000, 3000, 4000),
               ffpe = c(100, 250, 500, 750, 1000))

p_spatial_counts <- lapply(c("ffpe", "fresh_frozen"), function(dataset){
  ext <- ifelse(dataset == "fresh_frozen", "_FF", "")
  
  visium_obj <- readRDS(paste0(data_path, ext, "/Visium_HD_MouseBrain", ext, "_008um.rds"))
  
  square_size <- visium_obj@images[[paste0("slice1.008um")]]@scale.factors$spot
  ncounts_features_df <- visium_obj@meta.data %>% 
    rownames_to_column("cell") %>%
    left_join(GetTissueCoordinates(visium_obj), by="cell") %>% 
    mutate(square_size = square_size)
  # Get outliers
  quant <- visium_obj$nCount_Spatial.008um %>% 
    # Get 90th percentile
    quantile(0.99)
  
  ggplot(mapping = aes(x = y, y = x)) +
    geom_tile(mapping = aes(fill = label),
              height=square_size, width=square_size,
              data = ncounts_features_df %>% filter(nCount_Spatial.008um < 100) %>% 
                mutate(label="<100")) +
    scale_fill_manual(values = "lightgrey", name = NULL) +
    new_scale_fill() + 
    geom_tile(mapping = aes(fill = nCount_Spatial.008um),
              height=square_size, width=square_size,
              data = ncounts_features_df %>% filter(nCount_Spatial.008um >= 100)) +
    scale_fill_viridis_c(limits=c(100, quant), oob=scales::squish,
                         breaks = breaks_list[[dataset]],
                         name = "Total counts") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void(base_size=7) +
    theme(#legend.text = element_text(size=6),
          legend.key.size = unit(0.4, "cm"))
})

top_row <- p_spatial_counts[[1]] + p_spatial_counts[[2]] &
  theme(legend.key.size = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title.position = "top")
p_all <- top_row / p_st + plot_layout(heights = c(2, 1))
p_all
ggsave("visium_hd_brain_combined/plots/spatial_overview_brain_counts_hist.pdf",
       p_all, 
       width = 150, height = 150, units = "mm")
