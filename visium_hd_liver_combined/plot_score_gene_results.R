library(Seurat)
library(tidyverse)
library(patchwork)

color_palette <- c("Hepatocytes" = "#B4B5B5FF",
                   "Endothelialcells" = "#dca854",
                   "CentralVeinEndothelialcells" = "#FED8B1FF",
                   "LSECs" = "#FBB05FFF",
                   "PortalVeinEndothelialcells" = "#CC7722FF",
                   "LymphaticEndothelialcells" = "#8F4716FF",
                   "Cholangiocytes" = "#C61B84FF",
                   "HsPCs" = "#F19FC3FF",
                   "Stromalcells"="#77161d",
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
                   "Neutrophils" = "#727272",
                   "Unknown"="#fbf4a7")

celltype_order <- names(color_palette)
proper_celltype_names <- celltype_order %>% 
  stringr::str_replace_all(., "MonocytesMonocytederivedcells", "Mono & Mono-derived cells") %>%
  stringr::str_replace_all(., "([A-Za-z]+)(cells)$", "\\1 \\2") %>% 
  setNames(celltype_order)

plot_path <- paste0("visium_hd_liver_combined/plots/")
first_run <- FALSE
bin_sizes <- c(8, 16, 32)
rois <- list("caw009" = c(6600, 7100, 7200, 8200), # roi1
             "caw009" = c(9400, 9900, 11200, 12200), # roi2 (spp1)
             "sca002" = c(14000,15000, 12500, 14500), # roi1
             "sca002" = c(17500, 18500, 29000, 31000) # roi2 (spp1)
)



##### COMBINE PLOTS ####

dataset_name <- "caw009"
bin_size <- 8
bin_size_str <- sprintf("%03dum", bin_size)
gene_filters <- c("all", "q50", "q90")
score_type <- c("raw", "iterative")

for (dataset_name in c("sca002", "caw009")){
  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  rois_oi <- rois[names(rois) == dataset_name]
  roi_names <- sapply(rois_oi, function(roi) paste0("_", paste0(sprintf("%.1f", roi/1000), collapse="-")))
  score_path <- paste0("visium_hd_liver", dataset, "/score_genes")
  
  # Get proportions plot
  scores_df_summ <- lapply(gene_filters, function(gene_filter) {
    lapply(c("", "_iter"), function(score) {
      path <- file.path(score_path, paste0("square_008um_atlas_", gene_filter, "_celltype_scores", score, "_props.csv"))
      read.csv(path) %>% rename(spot=X) %>% 
        as.data.frame() %>%
        mutate(gene_filter = gene_filter,
               score = ifelse(score == "", "raw", "iterative")) %>% 
        group_by(gene_filter, score, annotation) %>%
        summarise(counts = n()) %>%
        ungroup() %>%
        mutate(proportion = counts / sum(counts))
    }) %>% bind_rows()
  }) %>% bind_rows() %>% 
    mutate(annotation = stringr::str_replace_all(annotation, "MonocytesandMonocytesderivedcells", "MonocytesMonocytederivedcells")) %>%
    mutate(annotation = stringr::str_replace_all(annotation, "unknown_celltype", "Unknown")) %>% 
    rename(celltype = annotation) %>% 
    mutate(celltype = droplevels(factor(celltype, levels = celltype_order)))
  
  # Get scores in ROI
  scores_df <- lapply(1:length(roi_names), function(i) {
    visium_obj_roi <- readRDS(paste0("visium_hd_liver_combined/rds/",
                                     "liver", toupper(dataset), "_",
                                     bin_size_str, "_ROI", roi_names[i], ".rds"))
    square_size <- visium_obj_roi@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
    
    lapply(gene_filters, function(gene_filter) {
      lapply(c("", "_iter"), function(score) {
        path <- file.path(score_path, paste0("square_008um_atlas_", gene_filter, "_celltype_scores", score, "_props.csv"))
        read.csv(path) %>% rename(spot=X) %>% 
          as.data.frame() %>%
          mutate(gene_filter = gene_filter,
                 score = ifelse(score == "", "raw", "iterative")) %>% 
          filter(spot %in% colnames(visium_obj_roi)) %>% 
          left_join(GetTissueCoordinates(visium_obj_roi), by = c("spot" = "cell"))
      }) %>% bind_rows()
    }) %>% bind_rows() %>% 
      mutate(roi=i, square_size=square_size, bin_size=bin_size,
             annotation = stringr::str_replace_all(annotation, "MonocytesandMonocytesderivedcells", "MonocytesMonocytederivedcells")) %>%
      mutate(annotation = stringr::str_replace_all(annotation, "unknown_celltype", "Unknown")) %>% 
      rename(celltype = annotation)
  }) %>% bind_rows()
  
  
  for (st in score_type){

    # Make barplots
    p_props_bar <- ggplot(scores_df_summ %>% filter(score==st),
                          aes(y=forcats::fct_rev(celltype), x=proportion, fill=celltype)) +
      geom_bar(stat="identity") +
      theme_minimal(base_size=8) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      scale_fill_manual(values=color_palette, guide="none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,1),
                         breaks = seq(0, 1, 0.2)) +
      facet_wrap(~gene_filter, ncol=1
      ) +
      labs(title = "Mean Cell Type Proportions") +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.border = ggh4x::element_part_rect(side = "lb", fill = NA, linewidth = 0.25),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(linewidth=0.25),
            axis.text.x = element_text(size=5, angle=0, hjust=0.5),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=7, face="bold"))
    # legend.position = "bottom",
    # legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    # legend.key.size = unit(0.4, "cm"),
    # legend.spacing.x = unit(0.4, "cm"),
    # legend.key.spacing.x = unit(0.3, "cm"),
    # legend.title = element_text(size=6, hjust=0.5, margin=margin(b=5)),
    # legend.text = element_text(size=6, margin = margin(l=3)))
    
  
    scores_df_square <- scores_df %>%
      filter(score == st) %>% 
      mutate(group = spot) %>% 
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
    
    scores_df_shapes <- scores_df_square %>%
      mutate(celltype = droplevels(factor(celltype, levels = celltype_order)))
    
    gene_filters <- c("all", "q50", "q90")
    gene_filter_names <- c("All markers", "Top 50%", "Top 10%") %>% setNames(gene_filters)
    gf_roi <- paste0(rep(gene_filters, each=length(roi_names)), "_", 1:length(roi_names))
    # Create the ggplot
    
    p_props_list <- lapply(gf_roi, function(txt) {
      gf <- strsplit(txt, "_")[[1]][1]
      roi_i <- as.numeric(strsplit(txt, "_")[[1]][2])
    
      scores_df_subset <- scores_df %>% 
        filter(score=="raw", gene_filter == gf, roi == roi_i)
      
      square_size <- scores_df_shapes$square_size[1]
      
      p <- ggplot(scores_df_shapes %>% filter(gene_filter == gf, roi == roi_i) ,
                  aes(x = x, y = y)) +
        geom_polygon(aes(fill = celltype, group = group), show.legend = TRUE) +
        # White border
        geom_tile(data = scores_df_subset %>% distinct(x, y),
                  aes(x = y, y = x), height = square_size, width = square_size,
                  fill = NA, color = "white", inherit.aes = FALSE) +
        scale_fill_manual(values = color_palette, labels=proper_celltype_names,
                            drop=FALSE) +
        labs(y = gene_filter_names[gf]) +
        theme_void(base_size=8) +
        scale_y_reverse() +
        coord_fixed() +
        guides(fill = guide_legend(nrow=2)) +
        theme(legend.title = element_blank(),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=6),
              legend.key.size = unit(0.5, "cm"),
              plot.title = element_text(hjust = 0.5, size=7))
      
      if (gf == "all") {
        p <- p + ggtitle(paste0("ROI ", roi_i))
      }
      
      if (roi_i == 1) {
        p <- p + theme(axis.title.y = element_text(size=7, angle=90))
      }
      
      p
    }) %>% setNames(gf_roi)

    design <- "ABG
               CDG
               EFG"
    
    wrapped_plots <- wrap_plots(p_props_list, ncol=2) + p_props_bar +
      plot_layout(guides = "collect", design=design) &
      theme(legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text = element_text(size=5, margin = margin(l=2, r=3)),
            legend.key.size = unit(0.35, "cm"),
            legend.key.spacing.y = unit(0.05, "cm"))
      
    ggsave(paste0(plot_path, "spatialtiles_barplot", dataset, "_score_", st, ".pdf"),
           wrapped_plots, width = 8, height = 6, bg = "white")

  }
}

