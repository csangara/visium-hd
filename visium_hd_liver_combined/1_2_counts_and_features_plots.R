# Supplementary Figure 4.3
source("visium_hd_liver_combined/0_utils.R")

data_path_sca <- paste0("data/Visium_HD_Liver/")
data_path_caw <- paste0("data/Visium_HD_Liver_CAW009/")

counts_features_df <- lapply(bin_sizes, function(bin_size) {
  bin_size_str <- sprintf("%03dum", bin_size)
  
  visium_obj_sca <- readRDS(paste0(data_path_sca, "Visium_HD_Liver_",
                                   bin_size_str, ".rds"))
  visium_obj_caw <- readRDS(paste0(data_path_caw, "Visium_HD_Liver_CAW009_",
                                    bin_size_str, ".rds"))
  
  bind_rows(
    data.frame(x = visium_obj_sca[[paste0("nCount_Spatial.", bin_size_str), drop=TRUE]], 
               dataset = "SCA", type = "Count"),
    data.frame(x = visium_obj_sca[[paste0("nFeature_Spatial.", bin_size_str), drop=TRUE]], 
               dataset = "SCA", type = "Feature"),
    data.frame(x = visium_obj_caw[[paste0("nCount_Spatial.", bin_size_str), drop=TRUE]],
               dataset = "CAW", type = "Count"),
    data.frame(x = visium_obj_caw[[paste0("nFeature_Spatial.", bin_size_str), drop=TRUE]],
               dataset = "CAW", type = "Feature")
  ) %>% mutate(
    bin_size = paste0(bin_size, "um")
  )
}) %>% bind_rows()

set1_colors <- RColorBrewer::brewer.pal(8, "Set1")[1:2]

# Factor 8, 16, 32
counts_features_df <- counts_features_df %>% 
  mutate(bin_size = factor(bin_size, levels = c("8um", "16um", "32um")),
         dataset = factor(dataset, levels = c("SCA", "CAW")))

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
ggplot(counts_features_df, aes(x = x, fill=dataset)) +
  geom_histogram(color="white", alpha=0.5, bins=50, position="identity", linewidth=0.1) +
  geom_text(aes(x = text_x, y = text_y, color=dataset,
                label = paste0("Mean: ", mean_x, "\nn = ", n)),
            data = means_df, vjust = 1, hjust=0, size=2,
            show.legend = FALSE
            ) +
  facet_wrap(type ~ bin_size, scales = "free",
             labeller = labeller( .multi_line=FALSE)) +
  scale_fill_manual(values = set1_colors, name = "Dataset") +
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
        strip.text = element_text(size = 6))

ggsave(paste0(plot_path, "histogram_counts_features.pdf"), width = 8, height = 5)


