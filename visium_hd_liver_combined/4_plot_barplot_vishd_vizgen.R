source("visium_hd_liver_combined/0_utils.R")
library(ggpattern)

# Load in deconvolved and ground truth proportions
deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))
deconv_props_8um <- deconv_props_rank %>% filter(bin_size=="008um")

vizgen_props <- lapply(bin_sizes, function(bin_size) {
    read.csv(paste0("visium_hd_liver_combined/vizgen/vizgen_grid_", bin_size, "um_celltype_proportions.csv")) %>%
      pivot_longer(cols=-grid_id, names_to="celltype", values_to="proportion") %>%
      mutate(bin_size = bin_size)
  }) %>% bind_rows()
vizgen_props_8um <- vizgen_props %>% filter(bin_size==8)

deconv_props_summ <- deconv_props_8um %>% 
  mutate(dataset = factor(dataset, levels = c("sca002", "caw009"))) %>% 
  group_by(celltype, dataset, rank, doublet_type) %>% 
  summarise(total_counts = n())

vizgen_props_summ <- vizgen_props_8um %>% 
  group_by(grid_id) %>% 
  filter(proportion > 0) %>% 
  mutate(rank = data.table::frankv(proportion, order=-1),
         doublet_type = case_when(n() == 1 ~ "singlet",
                                  n() > 1 ~ "doublet")) %>% 
  filter(rank <= 2) %>% 
  group_by(celltype, rank, doublet_type) %>%
  summarise(total_counts = n())

# Combine proportions
props_summ <- bind_rows(
  deconv_props_summ %>% mutate(source = "vishd"),
  vizgen_props_summ %>% mutate(dataset = "vizgen", source = "vizgen") %>% 
    mutate(celltype = gsub("\\.", "", celltype)) %>% 
    mutate(celltype = gsub("Kuppfer", "Kupffer", celltype))) %>%
    group_annot_to_vizgen(celltype, "annot_vizgen") %>% 
    # Sum by new annotation
    group_by(annot_vizgen, dataset, rank, doublet_type) %>%
    summarise(total_counts = sum(total_counts)) %>%
    mutate(annot_vizgen = factor(annot_vizgen, levels = names(color_palette_vizgen)),
           dataset = factor(dataset, levels = c("vizgen", "sca002", "caw009"), labels = c("Vizgen", "SCA002", "CAW009")),
           doublet_type = factor(R.utils::capitalize(doublet_type), levels = c("Singlet", "Doublet")))

# % Hepatocytes
props_summ %>% filter(dataset == "Vizgen", annot_vizgen == "Hepatocytes", rank==1) %>% 
  mutate(prop_doublet = total_counts / sum(total_counts))

deconv_props_8um %>% ungroup %>% distinct(dataset, spot, doublet_type) %>% 
  count(dataset, doublet_type) %>% 
  group_by(dataset) %>%
  mutate(prop = n / sum(n))

#### Figure 4.6 ####
# Prepare extra lines for plotting
singlet_threshold <- props_summ %>% filter(rank==1) %>%
  group_by(dataset) %>% summarise(threshold_counts=sum(total_counts)) %>% 
  mutate(text = case_when(dataset == "Vizgen" ~ "Total bins",
                          TRUE ~ ""),
         rank = 1)

extra_lines <- props_summ %>% 
  group_by(dataset, annot_vizgen, rank) %>% 
  filter(n() > 1, doublet_type=="Singlet") %>%
  left_join(singlet_threshold) %>% 
  filter(total_counts > 0.001*threshold_counts) %>% 
  mutate(y_pos = as.numeric(forcats::fct_rev(annot_vizgen))) %>% 
  mutate(y_pos = case_when(dataset != "Vizgen" ~ y_pos -1,
                           TRUE ~ y_pos))

# Pattern fill is the color of color_palette_vizgen and darker by 10%
fill_pattern_color <- sapply(color_palette_vizgen, function(col) {
  col_rgb <- grDevices::col2rgb(col)
  # For unknown, rgb = 25 each
  if (sum(col_rgb) > 100) {
    # Darker color
    new_rgb <- pmax(col_rgb - 30, 0)
  } else {
    # Do a lighter color for unknown
    new_rgb <- pmin(col_rgb + 50, 255)
  }
  grDevices::rgb(new_rgb[1], new_rgb[2], new_rgb[3], maxColorValue = 255)
}) %>% setNames(names(color_palette_vizgen))

p_barplot <- ggplot(props_summ) +
    ggpattern::geom_col_pattern(aes(pattern=doublet_type,
                                    x=total_counts,
                                    y=forcats::fct_rev(annot_vizgen),
                                    fill=annot_vizgen,
                                    pattern_fill=annot_vizgen),
                                    width=1,
                                    pattern_density = 0.25,
                                    pattern_color = NA,
                                    pattern_angle = 40,
                                    position = position_stack(reverse = TRUE)) +
    geom_linerange(data = extra_lines,
                   aes(x=total_counts, ymin=y_pos-0.5,
                       ymax=y_pos+0.5,
                       color=annot_vizgen),
                       linewidth=0.25) +
    geom_vline(data = singlet_threshold,
               aes(xintercept=threshold_counts),
               linetype="dashed", color="navyblue", linewidth=0.25) +
    geom_text(data = singlet_threshold,
              aes(x=threshold_counts, y=(length(color_palette_vizgen)/2)+0.5, label=text),
              vjust=-1, angle=90, size=2, color="navyblue") +
    ggh4x::facet_grid2(dataset~rank, scales = "free", independent = "x",
                       labeller=labeller(rank = c("1" = "Primary cell type",
                                                  "2" = "Secondary cell type"))) +
  scale_fill_manual(values = color_palette_vizgen, name = "Celltype",
                      labels=add_space_celltype, guide="none") +
    scale_pattern_manual(values=c("Singlet" = "stripe",
                                  "Doublet" = "none"),
                         labels = c("Singlet", "Doublet"),
                         name = "Doublet type",
                         guide = guide_legend(override.aes = list(fill="white",
                                                                  pattern_key_scale_factor = 0.5,
                                                                  color="black"), nrow=1)) +
    scale_y_discrete(labels = add_space_celltype) +
    # Pattern fill+ extra line
    scale_pattern_fill_manual(values = fill_pattern_color, guide = "none") +
    scale_color_manual(values = fill_pattern_color, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = scales::label_number(scale_cut=scales::cut_short_scale()),
                       name="Number of spots") +
    theme_bw(base_size=8) +
    theme_barplot_facet +
    theme(strip.text.y.right = element_blank(),
          strip.text.x.top = element_text(margin = margin(t = 0, r = 0, b = 5, l = 0), hjust=0),
          panel.spacing.x = unit(0.5, "cm"),
          axis.text.x = element_text(size=5),
          axis.text.y = element_text(size=6))
          

deconv_props_8um_vizgen <- deconv_props_8um %>% 
  group_annot_to_vizgen(celltype, "annot_vizgen") %>% 
  group_by(annot_vizgen, dataset, spot) %>%
  summarise(sum_props = sum(proportion))

props_combined <- bind_rows(
  deconv_props_8um_vizgen,
  vizgen_props %>% 
    filter(proportion > 0) %>% 
    select(-bin_size) %>% 
    mutate(celltype = gsub("\\.", "", celltype),
           celltype = gsub("Kuppfer", "Kupffer", celltype),
           dataset = "vizgen",
           grid_id = as.character(grid_id)) %>%
    rename(spot = grid_id,
           annot_vizgen = celltype,
           sum_props = proportion) 
) %>%  mutate(annot_vizgen = factor(annot_vizgen, levels = rev(names(color_palette_vizgen))),
              dataset = factor(dataset, levels = c("vizgen", "sca002", "caw009"),
                               labels = c("Vizgen", "SCA002", "CAW009"))) %>% 
  group_by(dataset, annot_vizgen) %>%
  mutate(
    outlier_lwr = sum_props < quantile(sum_props, probs = 0.25) - IQR(sum_props) * 1.5,
    outlier_upr = sum_props > quantile(sum_props, probs = 0.75) + IQR(sum_props) * 1.5
  ) %>%
  ungroup

p_boxplot <- ggplot(props_combined,
              aes(x=sum_props, y=annot_vizgen,
                  fill=annot_vizgen, colour=annot_vizgen)) +
    geom_boxplot(show.legend = FALSE, linewidth=0.15, outlier.shape=NA) +
    ggh4x::facet_wrap2(~dataset, ncol=1, strip.position = "right", axes="x", scales="free_y",
                       )  +
    scale_fill_manual(values = color_palette_vizgen, name = "Celltype",
                      labels = add_space_celltype) +
    scale_colour_manual(values = c(rep("black", length(color_palette_vizgen)-1), "white") %>% 
                          setNames(names(color_palette_vizgen))) +
    geom_boxplot(color = "black", fill = NA, fatten = NULL, linewidth=0.15, outlier.shape=NA) +
                 #outlier.size=0.1) +
    geom_point(data = function(x) subset(x, outlier_lwr | outlier_upr), 
             position = 'jitter', shape=16, stroke=0, size=0.1, show.legend = FALSE) +
    # Plot mean as a dot
    stat_summary(fun = "mean", geom = "point", shape = 21, size = 0.75,
                 stroke=0.25, fill = "white", color="black",
                 show.legend = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, by=0.2),
                     name="Cell type proportion") +
    theme_bw(base_size=8) +
    theme_barplot_facet +
    theme(strip.background = ggh4x::element_part_rect(side = "l", fill = NA, linewidth=0.5),
          strip.text.y.right = element_text(angle=0),
          axis.text.x = element_text(size=5),
          axis.text.y = element_blank(),
          plot.margin = margin(l=5)) 
p_boxplot

p_combined <-  p_barplot + p_boxplot +
  plot_layout(widths = c(2, 1))
ggsave(paste0(plot_path, "barplot_boxplot_vishd_vizgen_8um.pdf"),
       width = 8, height = 6, dpi=300)