library(Seurat)
library(tidyverse)
library(ggpattern)
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

celltype_order <- names(color_palette)
plot_path <- paste0("visium_hd_liver_combined/plots/")
first_run <- FALSE
bin_sizes <- c(8, 16, 32)

#### FIRST RUN - SAVE DECONV PROPS ####

if (first_run){
for (dataset_name in c("sca002", "caw009")){
  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
  proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
  
  for (bin_size in bin_sizes){
    gc()
    bin_size_str <- sprintf("%03dum", bin_size)
    
    visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                                 bin_size_str, ".rds"))
    ext <- "_annot_cd45"
    
    deconv_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                      "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                      bin_size_str, ext),
                               header = TRUE)
    
    removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,paste0("nCount_Spatial.", bin_size_str)] < 100)]
    print(paste(dataset_name, bin_size_str, "removed spots:") )
    print(length(removed_spots))
    # Check if removed rows + leftover rows == total rows (yes)
    length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2]
    
    # Subset visium_obj to only include spots that were not removed
    visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]
    
    # Add rownames to deconv_props
    rownames(deconv_props) <- colnames(visium_obj_subset)
    
    saveRDS(deconv_props, paste0("visium_hd_liver_combined/rds/", "liver", toupper(dataset), "_",
                                     bin_size_str, "_deconv_props.rds"))
  }
}
  
  # Combine all deconv props
  dataset_binsize <- paste0(rep(c("", "_CAW009"), each=length(bin_sizes)), "_", sprintf("%03dum", bin_sizes))
  deconv_props_all <- lapply(dataset_binsize, function(x) {
    deconv_props <- readRDS(paste0("visium_hd_liver_combined/rds/liver", x, "_deconv_props.rds"))
    
    data.frame(deconv_props) %>% 
      rownames_to_column("spot") %>% 
      pivot_longer(cols = -spot, names_to = "celltype",
                   values_to = "proportion") %>% 
      mutate(dataset = ifelse(grepl("CAW009", x), "caw009", "sca002"),
             bin_size = gsub("_CAW009|_", "", x))
  }) %>% bind_rows
  
  deconv_props_rank <- deconv_props_all %>% filter(proportion > 0) %>% 
    group_by(spot, dataset, bin_size) %>% 
    mutate(rank = data.table::frankv(proportion, order=-1),
           doublet_type = case_when(n() == 1 ~ "singlet",
                                    n() == 2 ~ "doublet"))
  
  celltype_order <- deconv_props_rank %>% 
    group_by(celltype, dataset, bin_size, rank, doublet_type) %>% 
    summarise(total_counts = n()) %>% 
    arrange(-total_counts) %>% 
    pull(celltype) %>% unique
  
  deconv_props_rank <- deconv_props_rank %>% 
    mutate(new_celltype = case_when(
      celltype %in% celltype_order[1:9] ~ celltype,
      TRUE ~ "Other")) %>% 
    mutate(new_celltype = factor(new_celltype, levels = c("Other", rev(celltype_order[1:9]))))
  
  saveRDS(deconv_props_rank, paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))
}


#### PRINT TABLE OF MEAN PROPORTIONS ####
deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))

dataset_binsize <- deconv_props_rank %>% ungroup() %>% 
  distinct(dataset, bin_size) %>% unite(dataset_bin, sep = "_") %>% 
  pull(dataset_bin)

mean_props <- lapply(dataset_binsize, function(x) {
  ds <- str_split(x, "_")[[1]][1]
  bin <- str_split(x, "_")[[1]][2]
  deconv_props_rank %>% ungroup() %>% 
    filter(dataset == ds, bin_size == bin) %>% 
    select(spot, celltype, proportion) %>% 
    pivot_wider(names_from = celltype, values_from = proportion, values_fill = 0) %>% 
    pivot_longer(cols = -spot, names_to = "celltype", values_to = "proportion") %>% 
    group_by(celltype) %>%
    summarise(mean_proportion = mean(proportion)) %>% 
    mutate(dataset = ds,
           bin_size = bin)
}) %>% bind_rows() %>% 
  tidyr::complete(dataset, bin_size, celltype,
                  fill = list(mean_proportion = 0))

mean_props %>% arrange(desc(dataset), bin_size) %>%
  pivot_wider(names_from = c(dataset, bin_size), values_from = mean_proportion) %>% 
  # get 2 significant digits in scientific notation
  mutate(across(where(is.numeric), ~ signif(., digits = 2))) %>%
  # Replace Endothelialcells with ECs in celltype
  mutate(celltype = stringr::str_replace_all(celltype, "Endothelialcells", "ECs")) %>% 
  # Add space before the word "cells"
  mutate(celltype = stringr::str_replace_all(celltype, "([A-Za-z]+)(cells)$", "\\1 \\2")) %>%
  write_csv(paste0("visium_hd_liver_combined/tables/deconv_props_mean.csv"))

## Calculate other statistics ##
# Calculate mean per cell type, dataset, bin size
deconv_props_rank %>% 
  group_by(celltype, dataset, bin_size) %>% 
  summarise(mean_proportion = mean(proportion)) %>% 
  arrange(celltype, dataset, bin_size) %>% 
  print(n=Inf)

# How many percent of spots are hepatocytes present in?
deconv_props_rank %>% 
  filter(doublet_type == "doublet") %>% 
  group_by(dataset, bin_size) %>% 
  summarise(total_spots = n_distinct(spot),
            total_hepatocytes = sum(celltype == "Hepatocytes")) %>% 
  mutate(percent_hepatocytes = total_hepatocytes / total_spots * 100)

#### PLOT BARPLOTS ####
deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))

deconv_props_summ <- deconv_props_rank %>% 
  mutate(dataset = factor(dataset, levels = c("sca002", "caw009"))) %>% 
  group_by(celltype, dataset, bin_size, rank, doublet_type) %>% 
  summarise(total_counts = n())

celltype_order <- deconv_props_summ %>% 
  arrange(-total_counts) %>% 
  pull(celltype) %>% unique

deconv_props_agg <- deconv_props_summ %>% 
  mutate(new_celltype = case_when(
    celltype %in% celltype_order[1:9] ~ celltype,
    TRUE ~ "Other")) %>% 
  group_by(new_celltype, dataset, bin_size, rank, doublet_type) %>%
  summarise(total_counts = sum(total_counts)) %>% 
  mutate(new_celltype = factor(new_celltype, levels = c("Other", rev(celltype_order[1:9]))),
         doublet_type = factor(R.utils::capitalize(doublet_type), levels = c("Singlet", "Doublet")))

celltype_format <- rev(levels(deconv_props_agg$new_celltype)) %>% 
  # Replace Endothelialcells with ECs
  stringr::str_replace_all("Endothelialcells", "ECs") %>%
  # Add space before the word "cells"
  stringr::str_replace_all("([A-Za-z]+)(cells)$", "\\1 \\2") %>% 
  # Add space before every capital letter that is not the first one
  stringr::str_replace_all("([a-z])([A-Z])", "\\1 \\2")
  # Wrap long celltype names (15)
  #stringr::str_wrap(width = 20)

new_col <- c(color_palette %>% .[names(.) %in% celltype_order[1:9]],
  "Other" = "#191919")

hep_threshold <- deconv_props_summ %>% filter(celltype=='Hepatocytes', rank==1) %>%
  group_by(dataset, bin_size) %>% summarise(threshold_counts=sum(total_counts)) %>% 
  mutate(text = case_when(
    bin_size == "008um" & dataset == "sca002" ~ "Total bins",
    TRUE ~ ""
  ))
extra_lines <- deconv_props_agg %>% 
  group_by(dataset, bin_size, new_celltype, rank) %>% 
  filter(n() > 1, doublet_type=="Singlet") %>%
  left_join(hep_threshold) %>% 
  filter(total_counts > 0.001*threshold_counts)

datasets <- unique(deconv_props_summ$dataset) %>% sort(decreasing=TRUE)

barplots <- lapply(datasets, function(ds) {
  p <- ggplot(deconv_props_agg %>% filter(dataset==ds)) +
    #geom_bar(stat="identity") +
    ggpattern::geom_col_pattern(aes(pattern = doublet_type,
                                    x=total_counts,
                                    y=new_celltype,
                                    fill=forcats::fct_rev(new_celltype)),
                                width=1,
                                pattern_density = 0.25,
                                pattern_color = NA,
                                pattern_fill = "gray50",
                                pattern_angle = 40,
                                position = position_stack(reverse = TRUE)) +
    geom_linerange(data = extra_lines %>% filter(dataset==ds),
                   aes(x=total_counts, ymin=as.numeric(new_celltype)-0.5,
                       ymax=as.numeric(new_celltype)+0.5),
                   linewidth=0.25, color="gray50") +
    #facet_wrap(bin_size~rank, ncol=2, nrow=3, scales="free_x") +
    ggh4x::facet_grid2(bin_size~rank, scales = "free_x", independent = "x",
                       labeller=labeller(rank=c("1"="Primary cell type", "2"="Secondary cell type"),
                                         bin_size=c("008um" = "8\u00b5m",
                                                    "016um" = "16\u00b5m",
                                                    "032um" = "32\u00b5m")),
                       switch="y") +
    scale_fill_manual(values = new_col, name = "Celltype",
                      labels=celltype_format) +
    scale_pattern_manual(values=c("Singlet" = "stripe",
                                  "Doublet" = "none"),
                         name = "Doublet type",
                         guide = guide_legend(override.aes = list(fill="white",
                                                                pattern_key_scale_factor = 0.5,
                                                                color="black"
                                                                ))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::comma,
                     name="Number of spots") +
    guides(fill = guide_legend(nrow=2)) +
    theme_bw(base_size=8) +
    ggtitle(toupper(ds)) +
    theme(strip.background = element_blank(),
          strip.text.y.left = element_text(angle=0),
          strip.text.x.top = element_text(margin = margin(t = 0, r = 0, b = 5, l = 0), hjust=0),
          panel.grid.major.y = element_blank(),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=6),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=5),
          axis.text.y = element_blank(),
          plot.title = element_text(size=7, face="bold"),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 5),
          legend.key.size = unit(0.5, "cm"),
          legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.position = "bottom",
          legend.direction = "vertical")
  ggsave(paste0(plot_path, "deconv_barplot_", ds, ".pdf"), p,
                width=8, height=6)
  p
}) %>% setNames(datasets)

# Combine plots into one
barplots[["sca002"]] + barplots[["caw009"]] +
  # Title of each plot
  plot_layout(guides="collect") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.key.size = unit(0.4, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.spacing.x = unit(0.4, "cm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6, margin = margin(l=3)))
ggsave(paste0(plot_path, "deconv_barplot_combined.pdf"), width=8, height=5)

## Plot primary barplot by itself ##
p_primary <- ggplot(deconv_props_agg %>% filter(rank==1)) +
    ggpattern::geom_col_pattern(aes(pattern = doublet_type,
       x=total_counts,  y=new_celltype, fill=forcats::fct_rev(new_celltype)),
       width=1, pattern_density = 0.25, pattern_color = NA, pattern_fill = "gray50",
       pattern_angle = 40, position = position_stack(reverse = TRUE)) +
    geom_linerange(data = extra_lines,
                   aes(x=total_counts, ymin=as.numeric(new_celltype)-0.5,
                       ymax=as.numeric(new_celltype)+0.5),
                   linewidth=0.25, color="gray50") +
    geom_vline(data = hep_threshold,
               aes(xintercept=threshold_counts),
               linetype="dashed", color="navyblue", linewidth=0.25) +
    geom_text(data = hep_threshold,
              aes(x=threshold_counts, y=(length(celltype_format)/2)+0.5, label=text),
              vjust=-1, angle=90, size=2, color="navyblue") +
    ggh4x::facet_grid2(bin_size~dataset, scales = "free_x", independent = "x", switch="y",
                       labeller=labeller(bin_size=c("008um" = "8\u00b5m", "016um" = "16\u00b5m", "032um" = "32\u00b5m"),
                                         dataset=toupper)) +
    scale_fill_manual(values = new_col, name = "Celltype", labels=celltype_format) +
    scale_pattern_manual(values=c("Singlet" = "stripe", "Doublet" = "none"), name = "Doublet type",
       guide = guide_legend(override.aes = list(fill="white", pattern_key_scale_factor = 0.5, color="black"))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = scales::comma) +
    guides(fill = guide_legend(nrow=3, override.aes = list(pattern = c(rep("none", length(celltype_format)))))) +
    theme_bw(base_size=8) +
    ggtitle("Number of bins with primary cell type") +
    theme(strip.background = element_blank(),
          strip.text.y.left = element_text(angle=0),
          strip.text.x.top = element_text(margin = margin(t = 0, r = 0, b = 5, l = 0), hjust=0),
          panel.grid.major.y = element_blank(),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=5),
          axis.text.y = element_blank(),
          plot.title = element_text(size=6, face="bold"),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 5),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size=5),
          legend.text = element_text(size=5),
          legend.position = "bottom",
          legend.direction = "vertical")

ggsave(paste0(plot_path, "deconv_primary_barplot_combined.pdf"), 
       p_primary,
       width=4, height=5)


#### PLOT BOXPLOT WITH BARPLOT ####

deconv_props_rank <- readRDS(paste0("visium_hd_liver_combined/rds/deconv_props_all.rds"))

deconv_props_summ <- deconv_props_rank %>% 
  group_by(celltype, dataset, bin_size, rank, doublet_type) %>% 
  summarise(total_counts = n())

celltype_order <- deconv_props_summ %>% 
  arrange(-total_counts) %>% 
  pull(celltype) %>% unique

deconv_props_agg <- deconv_props_summ %>% 
  mutate(new_celltype = case_when(
    celltype %in% celltype_order[1:9] ~ celltype,
    TRUE ~ "Other")) %>% 
  group_by(new_celltype, dataset, bin_size, rank, doublet_type) %>%
  summarise(total_counts = sum(total_counts)) %>% 
  mutate(new_celltype = factor(new_celltype, levels = c("Other", rev(celltype_order[1:9]))),
         doublet_type = factor(R.utils::capitalize(doublet_type), levels = c("Singlet", "Doublet")))

celltype_format <- rev(levels(deconv_props_agg$new_celltype)) %>% 
  # Replace Endothelialcells with ECs
  stringr::str_replace_all("Endothelialcells", "ECs") %>%
  # Add space before the word "cells"
  stringr::str_replace_all("([A-Za-z]+)(cells)$", "\\1 \\2") %>% 
  # Add space before every capital letter that is not the first one
  stringr::str_replace_all("([a-z])([A-Z])", "\\1 \\2")

new_col <- c(color_palette %>% .[names(.) %in% celltype_order[1:9]],
             "Other" = "#191919")
datasets <- unique(deconv_props_summ$dataset) %>% sort(decreasing=TRUE)

barplots <- lapply(datasets, function(ds) {
  p <- ggplot(deconv_props_agg %>% filter(dataset==ds, rank==2)) +
    geom_col(aes(x=total_counts,
                 y=new_celltype,
                 fill=forcats::fct_rev(new_celltype)),
                 position = position_stack(reverse = TRUE)) +
    facet_wrap(~bin_size, ncol=1, scales="free_x") +
    scale_fill_manual(values = new_col, name = "Celltype",
                      labels=celltype_format) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = scales::comma) +
    guides(fill = guide_legend(nrow=2)) +
    theme_bw(base_size=8) +
    theme(strip.background = element_blank(),
          strip.text.x.top = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=5),
          axis.text.y = element_blank(),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 5),
          legend.key.size = unit(0.5, "cm"),
          legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.position = "bottom",
          legend.direction = "vertical")
  # ggsave(paste0(plot_path, "deconv_barplot_", ds, ".pdf"), p,
  #        width=8, height=6)
  p
}) %>% setNames(datasets)


boxplots <- lapply(datasets, function(ds) {
  p <- ggplot(deconv_props_rank %>% filter(dataset == ds),
         aes(x=proportion, y=new_celltype, fill=new_celltype, colour=new_celltype)) +
    geom_boxplot(show.legend = FALSE, linewidth=0.15, outlier.shape=NA) +
    ggh4x::facet_wrap2(~bin_size, ncol=1, strip.position = "left", axes="x",
               labeller = labeller(bin_size=c("008um" = "8\u00b5m",
                                              "016um" = "16\u00b5m",
                                              "032um" = "32\u00b5m")))  +
    scale_fill_manual(values = new_col, name = "Celltype",
                      labels=celltype_format) +
    scale_colour_manual(values = c(rep("black", length(celltype_format)-1), "white") %>%
                          setNames(c(celltype_order[1:(length(celltype_format)-1)], "Other"))) + 
    geom_boxplot(color = "black", fill = NA, fatten = NULL, linewidth=0.15, outlier.size=0.1) +
    # Plot mean as a dot
    stat_summary(fun = "mean", geom = "point", shape = 21, size = 0.75,
                 stroke=0.25, fill = "white", color="black",
                 show.legend = FALSE) +
    theme_bw(base_size=8) +
    theme(strip.background = element_blank(),
          strip.text.y.left = element_text(angle=0),
          panel.grid.major.y = element_blank(),
          panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=5),
          axis.text.y = element_blank())
  if (ds == "caw009") {
    p <- p + theme(strip.text.y.left = element_blank())
  }
  p
}) %>% setNames(datasets)

# BOX - BAR - BOX - BAR
combined_plot <- (boxplots[[1]] + ggtitle("SCA002", subtitle = "Proportion of present cell type")) +
  (barplots[[1]] + ggtitle(" ", subtitle="Number of spots with secondary cell type")) +
  (boxplots[[2]] + ggtitle("CAW009", subtitle = "Proportion of present cell type")) +
  (barplots[[2]] + ggtitle(" ", subtitle="Number of spots with secondary cell type")) +
  plot_layout(ncol=4, guides="collect") &
  theme(panel.grid.minor.x = element_line(linewidth=0.1),
        panel.grid.major.x = element_line(linewidth=0.2),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.x = element_line(linewidth=0.25),
        plot.title = element_text(size=7, face="bold"),
        plot.subtitle = element_text(size=6),
        legend.position = "bottom",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.key.size = unit(0.4, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.spacing.x = unit(0.4, "cm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6, margin = margin(l=3)))

ggsave(paste0(plot_path, "deconv_boxbarboxbarplot_combined.pdf"),
       combined_plot,
       width=8, height=5)

# BOX - BOX - BAR - BAR
p1 <- (boxplots[[1]] + ggtitle("SCA002")) +
  (boxplots[[2]] + ggtitle("CAW009")) &
  theme(plot.title = element_text(size=7),
        panel.grid.minor.x = element_line(linewidth=0.1),
        panel.grid.major.x = element_line(linewidth=0.2),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.x = element_line(linewidth=0.25))

p1 <- p1 + plot_annotation(title = "(a) Proportions of cell type across bins where it is present",
                  theme = theme(plot.title = element_text(size=7.5, face="bold", margin = margin(t=0, r=0, b=0, l=0))))

p2 <- (barplots[[1]] + ggtitle("SCA002")) +
  (barplots[[2]] + ggtitle("CAW009")) &
  theme(legend.position = "none",
        plot.title = element_text(size=7),
        panel.grid.minor.x = element_line(linewidth=0.1),
        panel.grid.major.x = element_line(linewidth=0.2),
        panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
        axis.ticks.x = element_line(linewidth=0.25))
p2 <- p2 + plot_annotation(title = "(b) Number of bins with secondary cell type",
                          theme = theme(plot.title = element_text(size=7.5, face="bold", margin = margin(t=0, r=0, b=0, l=0))))

# Grab legend of barplots[[1]]
legend_bar <- cowplot::get_legend(barplots[[1]] + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.key.size = unit(0.4, "cm"),
        legend.key.spacing.x = unit(0.3, "cm"),
        legend.spacing.x = unit(0.4, "cm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6, margin = margin(l=3))))

p1p2 <- (wrap_elements(p1) + wrap_elements(p2)) +
  plot_layout(ncol=2, widths = c(1.05, 0.95))

p_combined <- p1p2 / wrap_elements(legend_bar) +
  plot_layout(nrow=2, heights=c(2, 0.2))
p_combined
ggsave(paste0(plot_path, "deconv_boxboxbarbarplot_combined.pdf"),
       p_combined,
       width=8, height=5.5)
