library(Seurat)
library(tidyverse)
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
rois <- list("caw009" = c(6600, 7100, 7200, 8200), # roi1
             "caw009" = c(9400, 9900, 11200, 12200), # roi2 (spp1)
             "sca002" = c(14000,15000, 12500, 14500), # roi1
             "sca002" = c(17500, 18500, 29000, 31000) # roi2 (spp1)
)


##### FIRST RUN: SAVE ROI OBJECTS #####

if (first_run){

for (dataset_name in c("sca002", "caw009")){
  cat("Processing dataset:", dataset_name, "\n")
  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
  proportions_path <- paste0("visium_hd_liver", dataset, "/Visium_HD_Liver", toupper(dataset))
  
  for (bin_size in bin_sizes){
    gc()
    bin_size_str <- sprintf("%03dum", bin_size)
    cat("Processing bin size:", bin_size_str, "\n")
    
    visium_obj <- readRDS(paste0(data_path, "Visium_HD_Liver", toupper(dataset), "_",
                                 bin_size_str, ".rds"))
    ext <- "_annot_cd45"
    
    deconv_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                      "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                      bin_size_str, ext),
                               header = TRUE)
    
    removed_spots <- Cells(visium_obj)[which(visium_obj@meta.data[,paste0("nCount_Spatial.", bin_size_str)] < 100)]
    
    # Check if removed rows + leftover rows == total rows (yes)
    length(removed_spots) + dim(deconv_props)[1] == dim(visium_obj)[2]
    
    # Print number of removed spots, total spots, and percentage
    cat("Number of removed spots:", length(removed_spots), "\n")
    cat("Total spots:", dim(visium_obj)[2], "\n")
    cat("Percentage of removed spots:", 
        round(length(removed_spots) / dim(visium_obj)[2] * 100, 2), "%\n")
    
    # Go to next loop
    next
    
    # Subset visium_obj to only include spots that were not removed
    visium_obj_subset <- visium_obj[, !(colnames(visium_obj) %in% removed_spots)]
    
    # Add rownames to deconv_props
    rownames(deconv_props) <- colnames(visium_obj_subset)
    
    # Assign barcode to most abundant cell type per spot
    visium_obj_subset$celltype <- colnames(deconv_props)[max.col(deconv_props)] %>% 
      factor(levels = celltype_order)
    
    rois_oi <- rois[names(rois) == dataset_name]
    
    for (roi in rois_oi) {
      
      roi_name <- paste0("_", paste0(sprintf("%.1f", roi/1000), collapse="-"))
      
      cells_roi <- GetTissueCoordinates(visium_obj_subset) %>%
        filter(x > roi[1] & x < roi[2] & y > roi[3] & y < roi[4]) %>% pull(cell)
      
      visium_obj_roi <- visium_obj_subset %>% .[, colnames(.) %in% cells_roi]
      deconv_props_roi <- deconv_props[rownames(deconv_props) %in% Cells(visium_obj_roi), ]
      
      # Save object
      saveRDS(visium_obj_roi, paste0("visium_hd_liver_combined/rds/", "liver", toupper(dataset), "_",
                                     bin_size_str, "_ROI", roi_name, ".rds"))
      saveRDS(deconv_props_roi, paste0("visium_hd_liver_combined/rds/", "liver", toupper(dataset), "_",
                                     bin_size_str, "_ROI", roi_name, "_deconv_props.rds"))
    }
  }
}

##### FIRST RUN: SAVE HIRES ROI OBJECTS #####

bin_size <- 16
bin_size_str <- sprintf("%03dum", bin_size)

for (dataset_name in c("sca002", "caw009")){
  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  data_path <- paste0("data/Visium_HD_Liver", toupper(dataset), "/")
  rois_oi <- rois[names(rois) == dataset_name] %>% do.call(rbind, .) %>% 
    `colnames<-`(c("ymin", "ymax", "xmin", "xmax")) %>% 
    data.frame

  image_dir <- file.path(data_path, paste0("binned_outputs/square_", bin_size_str, "/spatial"))
  
  image_hires <- Read10X_Image(
    image_dir,
    image.name = "tissue_hires_image.png",
    assay = paste0("Spatial.", bin_size_str),
  )
  
  visium_obj_hires <- Load10X_Spatial(
    data_path,
    bin.size = bin_size,
    image = image_hires,
  )
  
  # For rectangle border: scale and shift y-axis so it starts from the top
  y_range <- range(GetTissueCoordinates(visium_obj_hires)$x)
  hires_scalefactor <- visium_obj_hires@images[[paste0("slice1.", bin_size_str)]]@scale.factors$hires
  
  rois_oi_scaled <- rois_oi %>% 
    mutate(xmin = xmin * hires_scalefactor,
           xmax = xmax * hires_scalefactor,
           ymin = (y_range[2] - ymin + y_range[1])*hires_scalefactor,
           ymax = (y_range[2] - ymax + y_range[1])*hires_scalefactor,
           label = paste0("ROI ", row_number()))

  # Compare this with the plot we get
  # cells_roi1 <- GetTissueCoordinates(visium_obj_hires) %>%
  #   filter(x > rois_oi[1,1] & x < rois_oi[1,2] & y > rois_oi[1,3] & y < rois_oi[1,4]) %>% pull(cell)
  # cells_roi2 <- GetTissueCoordinates(visium_obj_hires) %>%
  #   filter(x > rois_oi[2,1] & x < rois_oi[2,2] & y > rois_oi[2,3] & y < rois_oi[2,4]) %>% pull(cell)
  # cells_roi <- c(cells_roi1, cells_roi2)
  # SpatialDimPlot(visium_obj_hires,
  #                cells.highlight = cells_roi,
  #                image.scale="hires")

  p_roi <- SpatialDimPlot(visium_obj_hires, pt.size.factor = 0,
                          image.scale="hires") +
    geom_rect(data = rois_oi_scaled,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax),
              fill = NA, color = "red", linewidth = 0.5,
              inherit.aes = FALSE) +
    geom_text(data = rois_oi_scaled,
              aes(x = (xmin + xmax) / 2, y = ymin,
                  label = label),
              color = "red", size = 3, vjust = -1,
              inherit.aes = FALSE)
  
  ggsave(paste0(plot_path, "spatialdimplot", toupper(dataset),
                "_ROI", roi_name, "_",
                bin_size_str, "_hires.png"), p_roi,
         width = 8, height = 8, bg = "white")
  
  saveRDS(p_roi, paste0("visium_hd_liver_combined/rds/",
                        "spatialdimplot_ROIs", toupper(dataset), "_",
                 bin_size_str, ".rds"))
}
}  
##### COMBINE PLOTS ####

dataset_name <- "caw009"
bin_size <- 16
bin_size_str <- sprintf("%03dum", bin_size)
features <- c("Glul", "Hal", "Spp1", "Clec4f")

for (dataset_name in c("sca002", "caw009")){

  dataset <- ifelse(dataset_name == "caw009", "_caw009", "")
  rois_oi <- rois[names(rois) == dataset_name]
  roi_names <- sapply(rois_oi, function(roi) paste0("_", paste0(sprintf("%.1f", roi/1000), collapse="-")))
  
  # Load the ROI dimplot
  p_roi <- readRDS(paste0("visium_hd_liver_combined/rds/",
                          "spatialdimplot_ROIs", toupper(dataset), "_",
                          bin_size_str, ".rds"))
  
  square_size <- readRDS(paste0("visium_hd_liver_combined/rds/","liver", toupper(dataset), "_",
                                bin_size_str, "_ROI", roi_names[1], ".rds"))@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
  
  # Load the Visium ROI objects
  features_df <- lapply(1:length(roi_names), function(i) {
    roi_name <- roi_names[i]
    visium_obj_roi <- readRDS(paste0("visium_hd_liver_combined/rds/",
                                     "liver", toupper(dataset), "_",
                                     bin_size_str, "_ROI", roi_name, ".rds"))
    
    visium_obj_roi <- NormalizeData(visium_obj_roi)
    
    # Check if the object is loaded correctly
    print(paste("Loaded ROI object for", roi_name))
    
    features_df <- GetAssayData(visium_obj_roi, layer = "data",
                                assay = paste0("Spatial.", bin_size_str))[features,] %>% 
      t() %>% data.frame() %>% 
      rownames_to_column("cell") %>%
      inner_join(GetTissueCoordinates(visium_obj_roi), by = "cell") %>% 
      pivot_longer(cols = -c(cell, x, y), names_to = "feature", values_to = "value") %>% 
      mutate(roi = i)
    
    features_df
  }) %>% bind_rows()
  
  feature_roi <- paste0(rep(features, each=length(roi_names)), "_", 1:length(roi_names))
  
  # Feature plots
  p_features_list <- lapply (feature_roi, function(txt) {
    feature_name <- strsplit(txt, "_")[[1]][1]
    roi_i <- as.numeric(strsplit(txt, "_")[[1]][2])
  
    p <- ggplot(features_df %>% filter(feature == feature_name, roi == roi_i),
                aes(x = y, y = x)) +
      geom_tile(aes(fill = value), color="gray90", linewidth=0.2,
                height=square_size, width=square_size) +
      scale_fill_viridis_c(breaks=scales::breaks_extended(n=5),
                           limits = range(features_df$value),
                           name = "log1p counts") +
      labs(title = feature_name) +
      scale_y_reverse() +
      coord_fixed() +
      theme_void()
    
    p
  }) %>%  setNames(feature_roi)
  
  
  p_roi <- p_roi + theme(legend.position="none")
  
  p_roi_border <- p_roi +
    coord_fixed(xlim = c(400, 4300)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=5))
  
  # Feature plots
  design <- "
  55
  12
  34
  "
  
  for (roi_i in 1:length(roi_names)){
    p_featureplots <- p_features_list[[paste0("Glul_", roi_i)]] +
      p_features_list[[paste0("Hal_", roi_i)]] +
      p_features_list[[paste0("Spp1_", roi_i)]] +
      p_features_list[[paste0("Clec4f_", roi_i)]] +
      p_roi_border +
      plot_layout(design = design, guides='collect')
    
    ggsave(paste0(plot_path, "featureplots_ROI", toupper(dataset), "_",
                  bin_size_str, "_roi",  roi_i, ".png"),
           p_featureplots, width = 8, height = 11, bg = "white")
  }
  
  
  # Load proportions
  deconv_props_df <- 
    lapply(bin_sizes, function(bin_size){
      bin_size_str <- sprintf("%03dum", bin_size)
      lapply(1:length(roi_names), function(i) {
        
        visium_obj_roi <- readRDS(paste0("visium_hd_liver_combined/rds/",
                                         "liver", toupper(dataset), "_",
                                         bin_size_str, "_ROI", roi_names[i], ".rds"))
        square_size <- visium_obj_roi@images[[paste0("slice1.", bin_size_str)]]@scale.factors$spot
        
        deconv_props_roi <- readRDS(paste0("visium_hd_liver_combined/rds/", "liver", toupper(dataset), "_",
                                           bin_size_str, "_ROI", roi_names[i], "_deconv_props.rds"))
        
        deconv_props_roi %>% 
          rownames_to_column("spot") %>%
          pivot_longer(cols = -spot, names_to = "celltype",
                       values_to = "proportion") %>% 
          filter(proportion > 0) %>% 
          # merge with coordinates
          left_join(GetTissueCoordinates(visium_obj_roi), by = c("spot" = "cell")) %>% 
          # the coords are x and y centroid, so we want x1, y1, x2, y2 as corners of squares
          # get whether there are one or two cell types
          mutate(n_celltypes = n(), .by = "spot",
                 square_size = square_size,
                 bin_size = bin_size,
                 roi = i)
      }) %>% bind_rows()
    }) %>% bind_rows()
  
  
  deconv_props_df_square <- deconv_props_df %>%
    filter(n_celltypes == 1) %>% 
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
  
  # Create barplot in squares
  deconv_props_df_barplot <- deconv_props_df %>% 
    filter(n_celltypes == 2) %>% 
    group_by(spot) %>% arrange(spot, desc(proportion)) %>% 
    mutate(rank = row_number(),
           group = paste0(spot, "_", rank)) %>% 
    mutate(x1 = case_when(rank == 1 ~ y - (square_size / 2) + (proportion*square_size),
                          rank == 2 ~ y + (square_size / 2) - (proportion*square_size)),
           y1 = x - square_size / 2,
           x2 = case_when(rank == 1 ~ y - (square_size / 2) + (proportion*square_size),
                          rank == 2 ~ y + (square_size / 2) - (proportion*square_size)),
           y2 = x + square_size / 2,
           x3 = case_when(rank == 1 ~ y - square_size / 2,
                          rank == 2 ~ y + square_size / 2),
           y3 = x + square_size / 2,
           x4 = case_when(rank == 1 ~ y - square_size / 2,
                          rank == 2 ~ y + square_size / 2),
           y4 = x - square_size / 2
    ) %>% 
    rename(coord_x = x, coord_y = y) %>%
    pivot_longer(cols = c(x1, y1, x2, y2, x3, y3, x4, y4),
                 names_to = c(".value", "corner"),
                 names_pattern = "(x|y)([1-4])")
  
  deconv_props_df_shapes <- bind_rows(deconv_props_df_square, deconv_props_df_barplot) %>% 
    mutate(celltype = factor(celltype, levels = celltype_order))
  
  binsize_roi <- paste0(rep(bin_sizes, each=length(roi_names)), "_", 1:length(roi_names))
  # Create the ggplot
  
  p_props_list <- lapply(binsize_roi, function(txt) {
    bs <- as.numeric(strsplit(txt, "_")[[1]][1])
    roi_i <- as.numeric(strsplit(txt, "_")[[1]][2])
    bin_size_str <- sprintf("%03dum", bs)
    deconv_props_df_subset <- deconv_props_df %>% 
      filter(bin_size == bs, roi == roi_i) 
    square_size <- deconv_props_df_subset$square_size[1]
    
    p <- ggplot(deconv_props_df_shapes %>% filter(roi == roi_i, bin_size == bs),
                aes(x = x, y = y)) +
      geom_polygon(aes(fill = celltype, group = group), show.legend = TRUE) +
      # White border
      geom_tile(data = deconv_props_df_subset %>% distinct(x, y),
                aes(x = y, y = x), height = square_size, width = square_size,
                fill = NA, color = "white", inherit.aes = FALSE) +
      scale_fill_manual(values = color_palette) +
      theme_void() +
      scale_y_reverse() +
      coord_fixed() +
      ggtitle(paste0(bs, "\u00b5m")) +
      guides(fill = guide_legend(ncol=1)) +
      theme(legend.title = element_blank(),
            legend.text = element_text(size=6),
            legend.key.size = unit(0.5, "cm"),
            plot.title = element_text(hjust = 0.5))
    p
  }) %>% setNames(binsize_roi)
  
  # Make custom legend
  legend_square_size <- 0.9
  
  for (roi_i in 1:length(rois_oi)){
    # Check that DCs are only in one resolution
    # dcs <- color_palette %>% .[. == names(which((color_palette %>% table) > 1))] %>% names
    # print(dcs)
    # deconv_props_df %>% filter(roi==roi_i, celltype %in% dcs) %>% distinct(celltype, bin_size)
    # 
    # celltype_order_legend <- celltype_order
    # celltype_order_legend[min(which(celltype_order %in% dcs))] <- "cDC1s, cDC2s, pDCs"
    # 
    spacing <- 0.25
    legend_df <- deconv_props_df %>% filter(roi == roi_i) %>%
      distinct(celltype, bin_size) %>%
      # mutate(celltype = case_when(celltype %in% dcs ~ "cDC1s, cDC2s, pDCs",
      #                             TRUE ~ as.character(celltype))) %>%
      # distinct(celltype, bin_size) %>% 
      mutate(x = case_when(bin_size == 8 ~ legend_square_size/6,
                           bin_size == 16 ~ legend_square_size/2,
                           bin_size == 32 ~ legend_square_size/2+legend_square_size/3),
             celltype = droplevels(factor(celltype, levels = celltype_order)),
             y = as.numeric(celltype) + (as.numeric(celltype)-1)*spacing, 
             celltype_label = stringr::str_replace_all(celltype, "Endothelialcells", "ECs")) %>%
      mutate(celltype_label = stringr::str_replace_all(celltype_label, "MonocytesMonocytederivedcells", "Mono & Mono-derived cells")) %>%
      mutate(celltype_label = stringr::str_replace_all(celltype_label, "([A-Za-z]+)(cells)$", "\\1 \\2")) %>% 
      mutate(celltype_label = stringr::str_replace_all(celltype_label, "([a-z]{2,})([A-Z])", "\\1 \\2"))
    
    legend_presence_df <- data.frame(
      x=c(legend_square_size/6, legend_square_size/2, legend_square_size/2+legend_square_size/3),
      y=c(-4, -3, -2),
      label=c("8\u00b5m", "16\u00b5m", "32\u00b5m")
    ) %>% mutate(y = y +(y+1)*spacing)
    
    p_legend_text_size <- 3
    p_legend <- ggplot(legend_df, aes(x = x, y = y)) +
      # Cell type legend
      geom_tile(aes(fill=celltype), height = legend_square_size, width = legend_square_size/3, show.legend = FALSE) +
      geom_text(data = legend_df %>% distinct(y, celltype_label), 
                aes(label=celltype_label, x=1.5), size=p_legend_text_size, hjust=0) +
      geom_text(data = data.frame(y=0, label="Cell types"), aes(label=label, y=y, x=0), size=p_legend_text_size, hjust=0) +
      # Presence legend
      geom_tile(data = legend_presence_df, aes(y=y), height=legend_square_size) +
      geom_tile(data = data.frame(y=legend_presence_df$y), aes(x=legend_square_size/2, y=y),
                height=legend_square_size, width=legend_square_size, fill=NA, color="black") +
      geom_text(data = legend_presence_df, aes(label=label, y=y, x=1.5), size=p_legend_text_size, hjust=0) +
      geom_text(data = data.frame(y=min(legend_presence_df$y)-1, label="Presence in ROI"),
                aes(label=label, y=y, x=0), size=p_legend_text_size, hjust=0) +
      scale_fill_manual(values = c(color_palette)) +
                                   #"cDC1s, cDC2s, pDCs" = "#893A86FF")) +
      xlim(0, 10) +
      theme_void() +
      scale_y_reverse() +
      coord_fixed()
    
    design <- "
    123
    145
    "
  
    row_1 <- p_roi + p_features_list[[paste0("Glul_", roi_i)]] + 
      p_features_list[[paste0("Hal_", roi_i)]] +
      p_features_list[[paste0("Spp1_", roi_i)]] +
      p_features_list[[paste0("Clec4f_", roi_i)]] +
      plot_layout(design=design, guides="collect") &
      theme(plot.title = element_text(hjust = 0.5, size=7),
            legend.key.size = unit(0.5, "cm"),
            legend.text = element_text(size=6),
            legend.title = element_text(size=6))
    row_1 <- wrap_elements(row_1)
    
    row_3 <- p_props_list[[paste0("16_", roi_i)]] + 
      p_props_list[[paste0("32_", roi_i)]] &
      theme(legend.position="none",
            plot.title = element_text(hjust = 0.5, size=7))
    
    row_23 <- wrap_elements(
      p_props_list[[paste0("8_", roi_i)]] / row_3 +
        plot_layout(heights=c(1.05, 0.5)) &
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5, size=7)))
    
    row_23_with_legend <- wrap_elements(row_23 + p_legend + 
                                          plot_layout(widths = c(0.75, 0.25)))
    
    
    all_rows <- row_1 / row_23_with_legend +
      plot_layout(heights = c(0.3, 0.7)) 
    
    ggsave(paste0(plot_path, "spatialbarplot_allres_featureplot", dataset, "_roi", roi_i, ".pdf"),
           all_rows, width = 8, height = 9, bg = "white")
    
  }
  
}

#### STOPPED HERE ####

ggsave(paste0(plot_path, "spatialbarplot_allres_firstrow", dataset, ".png"),
       first_row,
       width = 15, height = 8, bg = "white")

ggsave(paste0(plot_path, "spatialbarplot_allres_thirdrow", dataset, ".png"),
       third_row,
       width = 15, height = 5, bg = "white")

ggsave(paste0(plot_path, "spatialbarplot_allres_secondandthird", dataset, ".png"),
       second_and_third,
       width = 8, height = 10, bg = "white")

##### STOP HERE

ggplot(deconv_props_df_shapes %>% filter(roi==2, bin_size==16), aes(x = x, y = y)) +
  geom_polygon(aes(fill = celltype, group=group)) +
  # White border
  geom_tile(data = deconv_props_df %>% filter(roi==2, bin_size==16) %>%  distinct(x, y),
            aes(x = y, y = x), height = square_size, width = square_size,
            fill = NA, color = "white") +
  scale_fill_manual(values = color_palette) +
  scale_y_reverse() +
  coord_fixed() +
  theme_void() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, "cm"))

ggsave(paste0(plot_path, "spatialbarplot_doublet_celltype_ROI", roi_name, "_",
              bin_size_str, ".png"),
       width = 10, height = 8, bg = "white")

# Plot region of interest boundary in hires image
# Need to create a new Seurat object with hires image


SpatialDimPlot(visium_obj_roi,
               group.by = "celltype",
               image.alpha = 0, pt.size.factor = 14, shape=22, stroke=0.1) +
  scale_fill_manual(values = color_palette) +
  guides(fill = guide_legend(ncol=1, override.aes = list(size = 3))) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))

#p_celltype
ggsave(paste0(plot_path, "spatialdimplot_celltype_ROI", roi_name, "_", bin_size_str, ".png"),
       width = 10, height = 8, bg = "white") 

deconv_props_df_all <- data.frame(deconv_props) %>% 
  pivot_longer(cols = everything(), names_to = "celltype",
               values_to = "proportion")

# Cell type proportions barplot
deconv_props_summ <- deconv_props_df_all %>% 
  group_by(celltype) %>% 
  summarise(agg_proportion = mean(proportion))

# Barplot ordered by total abundance
p <- ggplot(deconv_props_summ, aes(x=reorder(celltype, agg_proportion), y=agg_proportion)) +
  geom_bar(stat="identity") +
  # Add text of value, rounded to two digits, only if value is > 0
  geom_text(aes(label=ifelse(agg_proportion > 0.001, round(agg_proportion, 3), "")), nudge_y = 0.03, size=2) +
  coord_flip() +
  theme_minimal(base_size = 8) +
  theme(axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  ggtitle("Average cell type proportions across tissue")

p
ggsave(paste0(plot_path, "barplot_avg_across_tissue_", bin_size_str, ".png"), p,
       width = 8, height = 6, bg = "white")


# What is the distribution of cell types across the tissue?
p_boxplot_all <- ggplot(deconv_props_df %>% filter(proportion > 0.0001),
                        aes(y=reorder(celltype, proportion, median), x=proportion)) +
  coord_flip() +
  geom_boxplot(alpha=0.1) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = 8) +
  theme(axis.title.x = element_blank(),
        # angled text
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),) +
  ggtitle("Cell type proportions across tissue (prop > 0.0001)")
p_boxplot_all
ggsave(paste0(plot_path, "boxplot_all_", bin_size_str, ".png"), p_boxplot_all,
       width = 15, height = 8, bg = "white")


# DOUBLET MODE
doublet_props <- read.table(paste0(proportions_path, "_", bin_size_str,
                                   "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                   bin_size_str, ext),
                            header = TRUE)

doublet_info <- read.table(paste0(proportions_path, "_", bin_size_str,
                                  "/proportions_rctd_Visium_HD_Liver", toupper(dataset), "_",
                                  bin_size_str, ext, "_doublet_info.tsv"),
                           header = TRUE)

# Add rownames to deconv_props
rownames(doublet_props) <- colnames(visium_obj_subset)
all(rownames(doublet_props) == doublet_info$spot) # Check

table(doublet_info$spot_class)

# List colors from RColorBrewer paired palette
# RColorBrewer::brewer.pal(n = 6, name = "Paired")

classes_colors <- c("singlet"="#377EB8", "doublet_certain"="#33A02C",
                    "doublet_uncertain"="#B2DF8A", "reject"="#E41A1C")
spotclass_order <- names(classes_colors)

visium_obj_subset$spot_class <- factor(doublet_info$spot_class, levels = spotclass_order)
p_spot_class <- SpatialDimPlot(visium_obj_subset, group.by = "spot_class",
                               image.alpha = 0, stroke = NA) +
  scale_fill_manual(values = classes_colors) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
p_spot_class
ggsave(paste0(plot_path, "spatialdimplot_spotclass_", bin_size_str, ".png"), p_spot_class,
       width = 6, height = 6, bg = "white")

# Let's just plot the ROIs
doublet_props_roi <- doublet_props[rownames(doublet_props) %in% cells_roi, ]
doublet_info_roi <- doublet_info[doublet_info$spot %in% cells_roi, ]

all(colnames(visium_obj_roi) == doublet_info_roi$spot)
visium_obj_roi$spot_class <- factor(doublet_info_roi$spot_class, levels = spotclass_order)

p_spot_class_roi <- SpatialDimPlot(visium_obj_roi, group.by = "spot_class",
                                   image.alpha = 0, pt.size.factor = 14, shape=22, stroke=0.1) +
  scale_fill_manual(values = classes_colors) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right") 
p_spot_class_roi

ggsave(paste0(plot_path, "spatialdimplot_spotclass_ROI_",
              bin_size_str, ".png"), p_spot_class_roi,
       width = 8, height = 6, bg = "white")

# Bar plot of spot classes over the whole tissue
ggplot(doublet_info, aes(y=spot_class)) +
  geom_bar(aes(fill=spot_class), width = 0.6, position="dodge") +
  scale_fill_manual(values = classes_colors) +
  scale_y_discrete(limits = rev(spotclass_order)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic(base_size = 8) +
  theme(axis.title = element_blank(),
        legend.position = "none") +
  ggtitle(paste0("Spot classes across ", bin_size_str, " tissue"))
ggsave(paste0(plot_path, "barplot_spotclasses_", bin_size_str, ".png"),
       width = 5, height = 5, bg = "white")
