use_color_palette <- "" # Use Karen's or my color palette
source("visium_hd_lung_ila010/0_utils.R")

# Supplementary Figure 4.14
bin_size <- 8 # 8, 16, or 32
bin_size_str <- sprintf("%03dum", bin_size)
visium_obj_ori <- readRDS(paste0(data_path, bin_size_str, ".rds"))
square_size <- visium_obj_ori@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot

ncounts_features_df <- visium_obj_ori@meta.data %>% 
  rownames_to_column("cell") %>%
  left_join(GetTissueCoordinates(visium_obj_ori), by="cell") %>% 
  mutate(square_size = square_size)

# Get outliers
quant <- visium_obj_ori$nCount_Spatial.008um %>% 
  # Get 90th percentile
  quantile(0.99)

p_scatterbar <- ggplot(mapping = aes(x = y, y = x)) +
  geom_tile(mapping = aes(fill = label),
            height=square_size, width=square_size,
            data = ncounts_features_df %>% filter(nCount_Spatial.008um <= 100) %>% 
              mutate(label="\u2264100")) +
  scale_fill_manual(values = "lightgrey", name = NULL) +
  new_scale_fill() + 
  geom_tile(mapping = aes(fill = nCount_Spatial.008um),
            height=square_size, width=square_size,
            data = ncounts_features_df %>% filter(nCount_Spatial.008um > 100)) +
  scale_fill_viridis_c(limits=c(100, quant), oob=scales::squish,
                       breaks = c(seq(0, 1000, 200)),
                       name = "Total counts") +
  scale_y_reverse() +
  coord_fixed() +
  theme_void(base_size=8) +
  theme(legend.text = element_text(size=6),
        legend.key.size = unit(0.4, "cm"))
p_scatterbar

# Histogram
ncounts_features_df_long <- 
  ncounts_features_df %>% select(cell, contains("Spatial.008um")) %>%
  pivot_longer(cols = contains("Spatial.008um"), 
               values_to = "x") %>% 
  mutate(type = gsub("_Spatial.008um", "", name)) %>% 
  mutate(type = gsub("^n", "", type))

means_df <- ncounts_features_df_long %>%
  group_by(type) %>%
  mutate(max_x = max(x)) %>%
  mutate(text_y = max(table(cut_width(x, width = unique(max_x)/50)))) %>%
  summarise(mean_x = prettyNum(round(mean(x)), big.mark=",", scientific=FALSE),
            n = prettyNum(n(), big.mark=",", scientific=FALSE),
            text_x = mean(x) + unique(max_x)/50*4, text_y = max(text_y))

set1_colors <- RColorBrewer::brewer.pal(8, "Set1")[2]

# Plot counts and features
p_histogram <- ggplot(ncounts_features_df_long, aes(x = x)) +
  geom_histogram(color="white", fill=set1_colors,
                 alpha=0.5, bins=50, position="identity", size=0.1) +
  geom_text(aes(x = text_x, y = text_y, color=set1_colors,
                label = paste0("Mean: ", mean_x, "\nn = ", n)),
            data = means_df, vjust = 1, hjust=0, size=2,
            show.legend = FALSE
  ) +
  facet_wrap(~ type, scales = "free",
             labeller = labeller( .multi_line=FALSE)) +
  scale_fill_manual(values = set1_colors) +
  scale_color_manual(values = set1_colors) +
  # Use scientific notation for y
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        # Add right spacing to legend title
        legend.title = element_text(margin = margin(r = 15), size=8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),
        strip.text = element_text(size = 6),
  )

design <- "AAAA
           AAAA
           #BB#"
p <- p_scatterbar + p_histogram + 
  plot_layout(design=design)
p
ggsave(p,
       filename = paste0(plot_path, "supplementary_spatial_deconv_overview_",
                         bin_size_str, "_a.pdf"),
       width = 10, height = 6, bg = "white")

# Old distribution plot - boxplot + violin plot
ori_meta_df <- visium_obj_ori@meta.data %>%
  mutate(
    outlier_lwr = nCount_Spatial.008um < quantile(nCount_Spatial.008um, probs = 0.25) - IQR(nCount_Spatial.008um) * 1.5,
    outlier_upr = nCount_Spatial.008um > quantile(nCount_Spatial.008um, probs = 0.75) + IQR(nCount_Spatial.008um) * 1.5
  )
ori_meta_df %>% 
  ggplot(aes(y=nCount_Spatial.008um, x=0)) +
  geom_violin(linewidth=0.15) +
  geom_boxplot(linewidth=0.15, outlier.shape=NA, width=0.1) +
  geom_point(data = function(x) subset(x, outlier_lwr | outlier_upr),
             position = 'jitter', shape=16, stroke=0, size=0.1, alpha=0.2) +
  # Plot mean as a dot
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 0.75,
               stroke=0.25, fill = "blue", color="blue",
               show.legend = FALSE) +
  theme_minimal(base_size=7) +
  scale_y_continuous(labels = scales::comma) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth=0.25),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title = element_blank())